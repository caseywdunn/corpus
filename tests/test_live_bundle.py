"""Opt-in all-tool acceptance against a real, read-only mounted bundle.

Requires Linux bubblewrap, Docker with nginx:stable-alpine already present,
and a cached query embedding model. CORPUS_TEST_BUNDLE selects a bundle with
taxonomy, bibliography, lexicon, vector rows and at least one pixel panel ROI.
Nothing modifies the selected bundle or downloads models. This is transport
and immutability acceptance, not an independent accuracy score.
"""
import asyncio
import base64
import io
import json
import os
from pathlib import Path
import secrets
import socket
import subprocess
import sys
import time
import uuid

import httpx
from PIL import Image
import pytest


def _port():
    with socket.socket() as sock:
        sock.bind(("127.0.0.1", 0))
        return sock.getsockname()[1]


def _inventory(root):
    return {str(p.relative_to(root)): (p.stat().st_size, p.stat().st_mtime_ns)
            for p in root.rglob("*") if p.is_file()}


def _discover(root):
    from mcpsrv.indexes import BiblioAuthority, CorpusIndex, TaxonMentionDB
    from pipeline.taxa import TaxonomyDB
    idx = CorpusIndex(root, biblio_db=BiblioAuthority(root / "biblio_authority.sqlite"),
                      taxonomy_db=TaxonomyDB(root / "taxonomy.sqlite"),
                      taxon_mention_db=TaxonMentionDB(root / "taxon_mentions.sqlite"))
    idx.load()
    assert idx.biblio_db and idx.taxonomy_db and idx.taxon_mention_db
    paper = next(p for p in idx.papers.values() if p["n_chunks"] and p["authors"] and p["year"])
    sha = paper["hash"]
    chunks = json.loads((root / "documents" / sha / "chunks.json").read_text())["chunks"]
    taxon = next(t["accepted_name"] for t in idx.taxon_display.values() if t.get("accepted_name"))
    category = next(cat for cat, terms in idx.lexicon_to_papers.items() if terms)
    term = next(iter(idx.lexicon_to_papers[category]))
    work = idx.biblio_db.get_work_by_corpus_hash(sha)
    assert work, "Selected paper must have a bibliographic work"
    panel = None
    for hd in sorted((root / "documents").iterdir()):
        path = hd / "figures.json"
        if not path.is_file():
            continue
        for figure in json.loads(path.read_text()).get("figures", []):
            for roi in figure.get("rois", []):
                if roi.get("roi_px") and roi.get("label"):
                    panel = {"paper_hash": hd.name, "figure_id": figure["figure_id"], "label": roi["label"]}
                    break
            if panel:
                break
        if panel:
            break
    assert panel, "Acceptance requires an actual pixel panel ROI, not a fallback"
    for db in (idx.biblio_db, idx.taxonomy_db, idx.taxon_mention_db):
        db.conn.close()
    return paper, chunks[0], taxon, category, term, work, panel


def _recipes(paper, chunk, taxon, category, term, work, panel):
    sha, author = paper["hash"], paper["authors"][0]["surname"]
    figure = {k: v for k, v in panel.items() if k != "label"}
    return {
        "bundle_info": {}, "corpus_summary": {"top_taxa": 2, "top_terms_per_category": 2},
        "list_papers": {"limit": 2}, "get_papers": {"hashes": [sha]},
        "get_chunks": {"paper_hash": sha, "chunk_ids": [chunk["chunk_id"]]},
        "get_chunks_by_section": {"paper_hash": sha, "limit": 2},
        "get_chunks_for_topic": {"query": paper["title"] or taxon, "k": 2},
        "get_bibliography": {"paper_hash": sha, "limit": 2, "resolved": True},
        "get_intext_citations": {"paper_hash": sha, "limit": 2},
        "get_excerpts_citing": {"work_id": work["work_id"], "limit": 2},
        "get_citation_graph": {"paper_hash": sha, "max_total_edges": 2},
        "resolve_reference": {"query": f"{author} {paper['year']}", "author": author, "year": paper["year"]},
        "format_citations": {"paper_hashes": [sha]},
        "get_missing_references": {"limit": 2},
        "get_original_description": {"taxon_name": taxon},
        "get_works_by_author": {"surname": author, "limit": 2},
        "search_taxon": {"name": taxon, "parent_chain": True},
        "get_papers_for_taxon": {"taxon_name": taxon},
        "get_chunks_for_taxon": {"taxon_name": taxon, "limit": 2},
        "get_taxon_mentions": {"taxon_name": taxon, "limit": 2},
        "get_papers_by_author": {"surname": author},
        "list_valid_species_under": {"parent_taxon_name": taxon},
        "get_taxon_dossier": {"taxon_name": taxon, "max_papers": 2, "max_chunks": 2, "max_figures": 2},
        "get_taxon_lexicon_slice": {"taxon_name": taxon, "category": category, "top_n": 2},
        "get_taxon_subtree_dossier": {"root_taxon_name": taxon, "max_species": 2, "max_aggregate_papers": 2},
        "lexicon_matrix": {"category": category, "terms": [term]},
        "get_lexicon_term_dossier": {"category": category, "term": term, "max_papers": 2},
        "get_figures_for_taxon": {"taxon_name": taxon, "limit": 2},
        "get_figures_for_lexicon_term": {"category": category, "term": term, "limit": 2},
        "get_figure_dossier_for_taxon": {"taxon_name": taxon, "max_figures": 2},
        "get_figure_dossier_for_term": {"category": category, "term": term, "max_figures": 2},
        "get_figure": figure, "list_figure_rois": figure,
        "get_figure_roi_image": {**panel, "profile": "report"},
        "get_figure_image": {**panel, "profile": "report"},
        "get_figure_url": {**figure, "profile": "report"},
        "get_active_profile": {}, "list_output_profiles": {},
    }


async def _exercise(port, proxy_port, token, recipes, panel):
    from mcp import ClientSession
    from mcp.client.sse import sse_client
    from tools.smoke_test_sse import _parse_tool_result
    seen = set()
    async with sse_client(f"http://127.0.0.1:{port}/sse", headers={"Authorization": f"Bearer {token}"}) as streams:
        async with ClientSession(*streams, read_timeout_seconds=300) as session:
            await session.initialize()
            listed = await session.list_tools()
            assert set(recipes) == {t.name for t in listed.tools}
            for name, args in recipes.items():
                result = await session.call_tool(name, args)
                assert not result.is_error, (name, result)
                if name == "get_figure_image":
                    blocks = [c for c in result.content if c.type == "image"]
                    assert blocks, result
                    with Image.open(io.BytesIO(base64.b64decode(blocks[0].data))) as image:
                        assert image.width > 0 and image.height > 0
                else:
                    parsed = _parse_tool_result(result)
                    # Empty list returns may have no text blocks. A nonempty
                    # plain-text reply is not a valid structured tool result.
                    assert parsed is None or isinstance(parsed, (dict, list)), (name, parsed)
                    if name in {"bundle_info", "corpus_summary", "get_papers", "get_chunks", "search_taxon",
                                "format_citations", "get_figure", "list_figure_rois", "get_figure_url"}:
                        assert parsed, (name, parsed)
                    items = parsed if isinstance(parsed, list) else [parsed]
                    assert all(not isinstance(item, dict) or not item.get("error") for item in items), (name, parsed)
                    if name == "get_chunks_for_topic":
                        assert items and all(isinstance(item, dict) and "score" in item for item in items), parsed
                    if name == "get_figure_roi_image":
                        assert parsed["crop"] is True, parsed
                seen.add(name)
                print(f"PASS live read-only MCP: {name}", flush=True)
            sizes = []
            async with httpx.AsyncClient() as client:
                for label in (None, panel["label"]):
                    args = {**panel, "label": label, "profile": "report"}
                    result = _parse_tool_result(await session.call_tool("get_figure_url", args))
                    assert result["auth_header"] is None and token not in json.dumps(result)
                    assert result["url"].startswith(f"http://127.0.0.1:{proxy_port}/figures/")
                    response = await client.get(result["url"])
                    assert response.status_code == 200
                    assert response.headers["cache-control"] == "private, no-store"
                    with Image.open(io.BytesIO(response.content)) as image:
                        sizes.append(image.size)
                    assert (await client.head(result["url"])).status_code == 200
                    assert (await client.get(result["url"] + "&label=tampered")).status_code == 401
            assert sizes[1][0] <= sizes[0][0] and sizes[1][1] <= sizes[0][1] and sizes[0] != sizes[1]
    assert seen == set(recipes)


def test_live_recipes_cover_the_frozen_tool_inventory():
    fixture = Path(__file__).parent / "fixtures/mcp_tool_contract.json"
    recipes = _recipes({"hash": "abc", "authors": [{"surname": "Author"}], "year": 2000, "title": "Title"},
                       {"chunk_id": "c0"}, "Taxon", "category", "term", {"work_id": "work"},
                       {"paper_hash": "abc", "figure_id": "f1", "label": "A"})
    contracts = json.loads(fixture.read_text())
    assert set(recipes) == set(contracts)
    for name, args in recipes.items():
        schema = contracts[name]["input_schema"]
        assert set(schema.get("required", [])) <= set(args), name
        assert set(args) <= set(schema.get("properties", {})), name


@pytest.mark.skipif(not os.environ.get("CORPUS_TEST_BUNDLE"), reason="opt-in real read-only bundle acceptance")
def test_all_tools_with_read_only_bundle_and_live_proxy(tmp_path):
    root = Path(os.environ["CORPUS_TEST_BUNDLE"]).resolve(strict=True)
    assert (root / "bundle_manifest.json").is_file(), "Pass the served bundle, not the build root"
    before = _inventory(root)
    inputs = _discover(root)
    recipes = _recipes(*inputs)
    port, proxy_port = _port(), _port()
    token = secrets.token_hex(32)
    token_path = tmp_path / "token"
    token_path.write_text(token)
    token_path.chmod(0o600)
    deployed = (Path(__file__).parents[1] / "deploy/nginx.conf").read_text()
    route = "location /figures/ {" + deployed.split("location /figures/ {", 1)[1].split("}", 1)[0] + "}"
    route = route.replace("127.0.0.1:8080", f"127.0.0.1:{port}")
    config = tmp_path / "nginx.conf"
    config.write_text("pid /tmp/nginx.pid; error_log /dev/stderr warn; events {} http { access_log off; "
                      "client_body_temp_path /tmp/client; proxy_temp_path /tmp/proxy; "
                      "fastcgi_temp_path /tmp/fastcgi; uwsgi_temp_path /tmp/uwsgi; scgi_temp_path /tmp/scgi; "
                      f"server {{ listen 127.0.0.1:{proxy_port}; {route} }} }}")
    name = "corpus-bundle-test-" + uuid.uuid4().hex
    started = False
    server = None
    env = {**os.environ, "HF_HUB_OFFLINE": "1", "TRANSFORMERS_OFFLINE": "1", "PYTHONDONTWRITEBYTECODE": "1",
           "TMPDIR": str(tmp_path)}
    try:
        with (tmp_path / "server.log").open("w") as log:
            # Everything is read-only except this test's temporary directory.
            # The selected bundle remains on the read-only root mount; private
            # crop caches can only live in TMPDIR, not alongside the evidence.
            server = subprocess.Popen([
                "bwrap", "--die-with-parent", "--ro-bind", "/", "/", "--dev", "/dev", "--proc", "/proc",
                "--bind", str(tmp_path), str(tmp_path), "--", sys.executable, "-m", "mcpsrv.main", str(root),
                "--transport", "sse", "--port", str(port), "--auth-token-file", str(token_path),
                "--public-base-url", f"http://127.0.0.1:{proxy_port}",
            ], env=env, stdout=log, stderr=log)
        deadline = time.monotonic() + 120
        while time.monotonic() < deadline:
            assert server.poll() is None, (tmp_path / "server.log").read_text()[-4000:]
            try:
                with socket.create_connection(("127.0.0.1", port), timeout=0.5):
                    break
            except OSError:
                time.sleep(0.2)
        else:
            pytest.fail("Server startup timed out")
        subprocess.run([
            "docker", "run", "--pull=never", "--rm", "-d", "--name", name, "--network=host", "--read-only",
            "--cap-drop=ALL", "--tmpfs", "/tmp", "--user", f"{os.getuid()}:{os.getgid()}", "--entrypoint", "nginx",
            "-v", f"{config}:/etc/nginx/nginx.conf:ro", "nginx:stable-alpine", "-e", "stderr", "-g", "daemon off;",
        ], check=True, capture_output=True)
        started = True
        deadline = time.monotonic() + 10
        while time.monotonic() < deadline:
            try:
                with socket.create_connection(("127.0.0.1", proxy_port), timeout=0.5):
                    break
            except OSError:
                time.sleep(0.1)
        else:
            pytest.fail("nginx startup timed out")
        from tools.smoke_test_sse import layer1_http
        assert layer1_http("127.0.0.1", port, token) == 0
        asyncio.run(_exercise(port, proxy_port, token, recipes, inputs[-1]))
    finally:
        if server:
            server.terminate()
            try:
                server.wait(timeout=10)
            except subprocess.TimeoutExpired:
                server.kill()
                server.wait()
        if started:
            subprocess.run(["docker", "stop", name], check=True, capture_output=True)
        assert _inventory(root) == before
