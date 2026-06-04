#!/usr/bin/env python3
"""Smoke-test the MCP server's SSE transport + bearer-token auth.

Runs three layers against a locally-launched ``python -m mcpsrv.main``:

  Layer 1 — raw HTTP: unauthenticated → 401, wrong-token → 401,
            correct-token → 200 with ``Content-Type: text/event-stream``.
            Catches middleware / uvicorn / binding issues.

  Layer 2 — MCP handshake: initialize, ``list_tools``,
            ``call_tool(bundle_info)``, ``call_tool(list_papers,
            limit=1)``.  Catches framing, keepalives, and tool-
            dispatch over SSE that stdio masks.

  Layer 3 — broad tool surface: exercises one representative call
            from each tool category (taxonomy, chunks/semantic search,
            bibliography, lexicon, figures, profiles).  Checks that
            tools return plausible non-error results without validating
            corpus-specific content — suitable for both the 4-paper
            demo and a full production corpuscle.

            Layer 3 includes ``get_chunks_for_topic`` (semantic search),
            which loads the ~600 MB BGE-M3 embedding model on first
            call.  Run on a compute node (``salloc`` / batch job) for
            full coverage; on a CPU-only login node the server may OOM.

Starts the server as a subprocess, waits for it to bind, runs the
checks, shuts it down cleanly.  No external deps — stdlib for layer
1, the ``mcp`` Python SDK (already a dep) for layers 2 and 3.

Usage:
    python tools/smoke_test_sse.py output
    python tools/smoke_test_sse.py output --port 18080 --skip-layer2
    python tools/smoke_test_sse.py output --skip-layer3
"""

from __future__ import annotations

import argparse
import asyncio
import json
import os
import secrets
import socket
import subprocess
import sys
import tempfile
import time
import urllib.error
import urllib.request
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
# v0.3 entry point: the old mcp_server.py root script was removed in
# the unified-CLI clean break (#60). The server now lives at
# mcpsrv.main and is run via `python -m mcpsrv.main`.
SERVER_MODULE = "mcpsrv.main"


# ── Tiny TAP-ish reporter ───────────────────────────────────────────

_USE_COLOR = sys.stdout.isatty()


def _color(s: str, code: str) -> str:
    return f"\033[{code}m{s}\033[0m" if _USE_COLOR else s


def _ok(msg: str) -> int:
    print(f"  {_color('PASS', '32')}  {msg}")
    return 0


def _fail(msg: str) -> int:
    print(f"  {_color('FAIL', '31')}  {msg}")
    return 1


def _info(msg: str) -> None:
    print(f"  {_color('•', '36')}     {msg}")


def _section(title: str) -> None:
    print(f"\n{_color('══', '34')} {title}")


# ── Server lifecycle ────────────────────────────────────────────────

def _wait_for_bind(host: str, port: int, timeout: float = 15.0) -> bool:
    deadline = time.time() + timeout
    while time.time() < deadline:
        try:
            with socket.create_connection((host, port), timeout=0.5):
                return True
        except OSError:
            time.sleep(0.2)
    return False


# ── Layer 1: raw HTTP ───────────────────────────────────────────────

def layer1_http(host: str, port: int, token: str) -> int:
    _section("Layer 1: raw HTTP + bearer auth")
    url = f"http://{host}:{port}/sse"
    rc = 0

    # 1a: unauthenticated GET → 401
    try:
        urllib.request.urlopen(url, timeout=3)
        rc |= _fail("unauth GET /sse did not raise")
    except urllib.error.HTTPError as e:
        if e.code == 401:
            _ok(f"unauth GET /sse → 401")
        else:
            rc |= _fail(f"unauth GET /sse → {e.code} (expected 401)")
    except Exception as e:
        rc |= _fail(f"unauth GET /sse raised {type(e).__name__}: {e}")

    # 1b: wrong-token GET → 401
    req = urllib.request.Request(url, headers={"Authorization": "Bearer nope"})
    try:
        urllib.request.urlopen(req, timeout=3)
        rc |= _fail("wrong-token GET did not raise")
    except urllib.error.HTTPError as e:
        if e.code == 401:
            _ok("wrong-token GET /sse → 401")
        else:
            rc |= _fail(f"wrong-token GET /sse → {e.code} (expected 401)")

    # 1c: authenticated GET → 200 + text/event-stream Content-Type
    # SSE is an open-ended stream; we only need to confirm the response
    # line + headers arrive cleanly, not drain the body.
    req = urllib.request.Request(
        url, headers={"Authorization": f"Bearer {token}"}
    )
    try:
        resp = urllib.request.urlopen(req, timeout=3)
    except urllib.error.HTTPError as e:
        rc |= _fail(f"auth GET /sse → {e.code}")
        return rc
    try:
        ctype = resp.headers.get("Content-Type", "")
        if resp.status == 200 and "text/event-stream" in ctype:
            _ok(f"auth GET /sse → 200 {ctype}")
        else:
            rc |= _fail(
                f"auth GET /sse → status={resp.status} ctype={ctype!r}"
            )
    finally:
        resp.close()

    return rc


# ── Layer 2: MCP handshake + core tools ────────────────────────────

async def layer2_mcp_client(host: str, port: int, token: str) -> int:
    _section("Layer 2: MCP client over SSE")
    try:
        from mcp import ClientSession
        from mcp.client.sse import sse_client
    except ImportError as e:
        _info(f"mcp SDK client not importable: {e}")
        _info("skipping layer 2 (install `mcp` or run with --skip-layer2)")
        return 0

    url = f"http://{host}:{port}/sse"
    headers = {"Authorization": f"Bearer {token}"}
    rc = 0

    try:
        async with sse_client(url, headers=headers) as (read, write):
            async with ClientSession(read, write) as session:
                await session.initialize()
                _ok("initialize handshake")

                # list_tools
                try:
                    listing = await session.list_tools()
                    tools = listing.tools
                    if len(tools) < 1:
                        rc |= _fail(f"list_tools returned {len(tools)} tools")
                    else:
                        names = sorted(t.name for t in tools)
                        _ok(f"list_tools → {len(tools)} tools "
                            f"(e.g. {', '.join(names[:3])}, …)")
                except Exception as e:
                    rc |= _fail(f"list_tools raised: {e}")

                # call bundle_info — returns a dict regardless of deploy
                # state.  Parse the first text block as JSON.
                try:
                    result = await session.call_tool("bundle_info", {})
                    parsed = _parse_tool_result(result)
                    if parsed is not None and "bundle_version" in parsed:
                        _ok(
                            f"bundle_info → bundle_version="
                            f"{parsed.get('bundle_version')!r}"
                        )
                    else:
                        rc |= _fail(
                            f"bundle_info result missing bundle_version: "
                            f"{parsed!r}"
                        )
                except Exception as e:
                    rc |= _fail(f"call_tool(bundle_info) raised: {e}")

                # call list_papers(limit=1) — cheap and should succeed
                # on any processed corpus.
                try:
                    result = await session.call_tool(
                        "list_papers", {"limit": 1}
                    )
                    parsed = _parse_tool_result(result)
                    if isinstance(parsed, list):
                        _ok(f"list_papers(limit=1) → list of {len(parsed)}")
                    else:
                        _ok(
                            f"list_papers(limit=1) returned "
                            f"{type(parsed).__name__}"
                        )
                except Exception as e:
                    rc |= _fail(f"call_tool(list_papers) raised: {e}")

    except Exception as e:
        rc |= _fail(f"SSE session setup failed: {type(e).__name__}: {e}")

    return rc


# ── Layer 3: broad tool surface ─────────────────────────────────────

async def layer3_tool_coverage(host: str, port: int, token: str) -> int:
    """Exercise one representative call from each tool category.

    Checks for non-error responses and expected top-level fields;
    does not validate corpus-specific content so the tests work on
    both the 4-paper demo and a full production corpuscle.
    """
    _section("Layer 3: broad tool surface")
    try:
        from mcp import ClientSession
        from mcp.client.sse import sse_client
    except ImportError as e:
        _info(f"mcp SDK not importable: {e} — skipping layer 3")
        return 0

    url = f"http://{host}:{port}/sse"
    headers = {"Authorization": f"Bearer {token}"}
    rc = 0

    async with sse_client(url, headers=headers) as (read, write):
        async with ClientSession(read, write) as session:
            await session.initialize()

            # ── Discover a real paper hash for paper-keyed tests ────
            # limit=1 returns one item → one content block → dict.
            # Handle both dict (single item) and list (multi) cases.
            first_hash = None
            try:
                r = await session.call_tool("list_papers", {"limit": 1})
                papers = _parse_tool_result(r)
                if isinstance(papers, list) and papers:
                    first_hash = papers[0].get("hash")
                elif isinstance(papers, dict):
                    first_hash = papers.get("hash")
            except Exception:
                pass

            # ── Corpus summary ──────────────────────────────────────
            _section("Layer 3a: corpus-level tools")
            try:
                r = await session.call_tool("corpus_summary", {})
                d = _parse_tool_result(r)
                # key is n_papers (not paper_count)
                if isinstance(d, dict) and d.get("n_papers", 0) > 0:
                    _ok(f"corpus_summary → {d['n_papers']} papers, "
                        f"{d.get('n_unique_taxa', '?')} taxa, "
                        f"{d.get('n_figures_total', '?')} figures")
                else:
                    rc |= _fail(f"corpus_summary unexpected: {d!r}")
            except Exception as e:
                rc |= _fail(f"corpus_summary raised: {e}")

            # ── Taxonomy tools ──────────────────────────────────────
            _section("Layer 3b: taxonomy tools")

            # search_taxon — use a broadly-valid name present in any
            # marine corpus; the important check is that the tool runs
            # without error and returns a dict with a "found" key.
            taxon_hit = None
            try:
                r = await session.call_tool(
                    "search_taxon", {"name": "Siphonophorae"}
                )
                d = _parse_tool_result(r)
                # success: {matched_taxon_id, accepted_name, ...}
                # not-found: {not_found: True, queried: name}
                if isinstance(d, dict) and "matched_taxon_id" in d:
                    taxon_hit = d.get("accepted_taxon_id") or "Siphonophorae"
                    _ok(f"search_taxon(Siphonophorae) → found, "
                        f"accepted={d.get('accepted_name')!r}, "
                        f"in_corpus={d.get('in_corpus')}")
                elif isinstance(d, dict) and d.get("not_found"):
                    _ok("search_taxon(Siphonophorae) → not found "
                        "(expected on demo corpus)")
                elif isinstance(d, dict) and "error" in d:
                    _ok(f"search_taxon → error (no taxonomy? {d['error']!r})")
                else:
                    rc |= _fail(f"search_taxon unexpected: {d!r}")
            except Exception as e:
                rc |= _fail(f"search_taxon raised: {e}")

            # list_valid_species_under
            try:
                r = await session.call_tool(
                    "list_valid_species_under",
                    {"parent_taxon_name": "Siphonophorae"},
                )
                d = _parse_tool_result(r)
                if isinstance(d, list):
                    _ok(f"list_valid_species_under(Siphonophorae) → "
                        f"{len(d)} species")
                elif isinstance(d, dict) and "error" in d:
                    _ok(f"list_valid_species_under → error (expected on "
                        f"demo: {d['error']!r})")
                elif isinstance(d, dict) and "accepted_taxon_id" in d:
                    # Single-item list serialised as one block → dict
                    _ok("list_valid_species_under → 1 species (dict form)")
                else:
                    rc |= _fail(f"list_valid_species_under unexpected: {d!r}")
            except Exception as e:
                rc |= _fail(f"list_valid_species_under raised: {e}")

            # get_papers_for_taxon
            try:
                r = await session.call_tool(
                    "get_papers_for_taxon",
                    {"taxon_name": "Physalia", "limit": 3},
                )
                d = _parse_tool_result(r)
                if isinstance(d, list):
                    _ok(f"get_papers_for_taxon(Physalia) → {len(d)} papers")
                elif isinstance(d, dict):
                    # wrapped response or error dict — acceptable
                    _ok(f"get_papers_for_taxon → dict ({list(d)[:3]})")
                else:
                    rc |= _fail(f"get_papers_for_taxon unexpected: {type(d)}")
            except Exception as e:
                rc |= _fail(f"get_papers_for_taxon raised: {e}")

            # get_taxon_dossier — the main workhorse dossier tool
            try:
                r = await session.call_tool(
                    "get_taxon_dossier",
                    {"taxon_name": "Physalia", "max_papers": 3},
                )
                d = _parse_tool_result(r)
                if isinstance(d, dict):
                    keys = set(d.keys())
                    expected = {"taxon", "papers"}
                    if keys & expected:
                        _ok(f"get_taxon_dossier(Physalia) → keys: "
                            f"{sorted(keys)[:6]}")
                    elif "error" in keys:
                        _ok(f"get_taxon_dossier → error dict (demo OK: "
                            f"{d['error']!r})")
                    else:
                        rc |= _fail(f"get_taxon_dossier unexpected keys: {keys}")
                else:
                    rc |= _fail(f"get_taxon_dossier unexpected type: {type(d)}")
            except Exception as e:
                rc |= _fail(f"get_taxon_dossier raised: {e}")

            # ── Chunks / semantic search ────────────────────────────
            _section("Layer 3c: chunk tools")

            # get_chunks_for_topic — semantic search, requires LanceDB
            # and loads the ~600 MB BGE-M3 embedder on first call.
            # Run on a compute node for full coverage.
            try:
                r = await session.call_tool(
                    "get_chunks_for_topic",
                    {"query": "nectophore morphology", "k": 3},
                )
                d = _parse_tool_result(r)
                if isinstance(d, list):
                    _ok(f"get_chunks_for_topic → {len(d)} chunks")
                    if d and isinstance(d[0], dict):
                        has_text = "text" in d[0] or "chunk_id" in d[0]
                        if not has_text:
                            rc |= _fail(
                                f"chunk missing text/chunk_id: {list(d[0])}"
                            )
                elif isinstance(d, dict) and "error" in d:
                    _ok(f"get_chunks_for_topic → error (no embeddings? "
                        f"{d['error']!r})")
                else:
                    rc |= _fail(f"get_chunks_for_topic unexpected: {d!r}")
            except Exception as e:
                rc |= _fail(f"get_chunks_for_topic raised: {e}")

            # get_chunks (paper-keyed, requires a real hash)
            if first_hash:
                try:
                    r = await session.call_tool(
                        "get_chunks",
                        {"paper_hash": first_hash, "with_text": False},
                    )
                    d = _parse_tool_result(r)
                    if isinstance(d, list):
                        _ok(f"get_chunks({first_hash[:8]}…) → "
                            f"{len(d)} chunks")
                    else:
                        rc |= _fail(f"get_chunks unexpected: {type(d)}")
                except Exception as e:
                    rc |= _fail(f"get_chunks raised: {e}")
            else:
                _info("skipping get_chunks (no paper hash available)")

            # ── Bibliography tools ──────────────────────────────────
            _section("Layer 3d: bibliography tools")

            # get_bibliography (paper-keyed)
            if first_hash:
                try:
                    r = await session.call_tool(
                        "get_bibliography",
                        {"paper_hash": first_hash, "limit": 5},
                    )
                    d = _parse_tool_result(r)
                    if isinstance(d, list):
                        _ok(f"get_bibliography({first_hash[:8]}…) → "
                            f"{len(d)} refs")
                    else:
                        rc |= _fail(f"get_bibliography unexpected: {type(d)}")
                except Exception as e:
                    rc |= _fail(f"get_bibliography raised: {e}")
            else:
                _info("skipping get_bibliography (no paper hash available)")

            # get_missing_references — exercises biblio_authority.sqlite
            try:
                r = await session.call_tool(
                    "get_missing_references", {"limit": 5}
                )
                d = _parse_tool_result(r)
                if isinstance(d, list):
                    _ok(f"get_missing_references → {len(d)} missing refs")
                elif isinstance(d, dict) and "error" in d:
                    _ok(f"get_missing_references → error (no biblio DB? "
                        f"{d['error']!r})")
                elif isinstance(d, dict) and "work_id" in d:
                    # Single-item list → one block → dict
                    _ok("get_missing_references → 1 ref (dict form)")
                else:
                    rc |= _fail(f"get_missing_references unexpected: {type(d)}")
            except Exception as e:
                rc |= _fail(f"get_missing_references raised: {e}")

            # get_works_by_author — a well-known siphonophore author
            try:
                r = await session.call_tool(
                    "get_works_by_author", {"surname": "Huxley"}
                )
                d = _parse_tool_result(r)
                if isinstance(d, list):
                    _ok(f"get_works_by_author(Huxley) → {len(d)} works")
                elif isinstance(d, dict):
                    _ok(f"get_works_by_author → dict ({list(d)[:3]})")
                else:
                    rc |= _fail(f"get_works_by_author unexpected: {type(d)}")
            except Exception as e:
                rc |= _fail(f"get_works_by_author raised: {e}")

            # ── Lexicon tools ───────────────────────────────────────
            _section("Layer 3e: lexicon tools")

            try:
                r = await session.call_tool(
                    "lexicon_matrix",
                    {"category": "anatomy"},
                )
                d = _parse_tool_result(r)
                if isinstance(d, dict):
                    if "error" in d:
                        _ok(f"lexicon_matrix → error (no lexicon? "
                            f"{d['error']!r})")
                    elif "term_totals" in d or "terms" in d:
                        # detail=False → {category, detail, paper_count, term_totals}
                        # detail=True  → {category, detail, terms, rows}
                        n = (len(d.get("term_totals") or d.get("terms") or []))
                        _ok(f"lexicon_matrix(anatomy) → {n} terms, "
                            f"paper_count={d.get('paper_count', '?')}")
                    else:
                        rc |= _fail(
                            f"lexicon_matrix unexpected keys: {sorted(d)}"
                        )
                else:
                    rc |= _fail(f"lexicon_matrix unexpected: {type(d)}")
            except Exception as e:
                rc |= _fail(f"lexicon_matrix raised: {e}")

            # ── Figures tools ───────────────────────────────────────
            _section("Layer 3f: figure tools")

            try:
                r = await session.call_tool(
                    "get_figures_for_taxon",
                    {"taxon_name": "Physalia", "limit": 3},
                )
                d = _parse_tool_result(r)
                if isinstance(d, list):
                    _ok(f"get_figures_for_taxon(Physalia) → {len(d)} figures")
                elif isinstance(d, dict) and "error" in d:
                    _ok(f"get_figures_for_taxon → error (demo OK: "
                        f"{d['error']!r})")
                elif isinstance(d, dict) and "figure_id" in d:
                    # Single figure serialised as one block → dict
                    _ok(f"get_figures_for_taxon(Physalia) → 1 figure (dict form)")
                else:
                    rc |= _fail(
                        f"get_figures_for_taxon unexpected: {type(d)}"
                    )
            except Exception as e:
                rc |= _fail(f"get_figures_for_taxon raised: {e}")

            # ── Profiles tools ──────────────────────────────────────
            _section("Layer 3g: profile tools")

            try:
                r = await session.call_tool("list_output_profiles", {})
                d = _parse_tool_result(r)
                if isinstance(d, dict) and "profiles" in d:
                    _ok(f"list_output_profiles → "
                        f"{len(d['profiles'])} profiles: "
                        f"{list(d['profiles'])}")
                elif isinstance(d, dict):
                    _ok(f"list_output_profiles → {sorted(d)[:5]}")
                else:
                    rc |= _fail(f"list_output_profiles unexpected: {type(d)}")
            except Exception as e:
                rc |= _fail(f"list_output_profiles raised: {e}")

    return rc


def _parse_tool_result(result):
    """Best-effort extraction of a JSON payload from an MCP tool
    result's content blocks.  Returns the parsed value or None.

    The MCP framework serialises Python list returns as one content
    block per item (not as a single JSON array block).  When we see
    multiple parseable blocks we re-assemble them into a list so
    callers can do ``isinstance(d, list)`` normally.
    """
    blocks = []
    for block in getattr(result, "content", []) or []:
        text = getattr(block, "text", None)
        if not text:
            continue
        try:
            blocks.append(json.loads(text))
        except Exception:
            # Plain-text block — keep it so we don't lose single
            # plain-text responses.
            blocks.append(text)
    if not blocks:
        return None
    if len(blocks) == 1:
        return blocks[0]
    # Multiple blocks → the server serialised a list one item per
    # block; re-assemble.
    return blocks


# ── Orchestration ───────────────────────────────────────────────────

def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("output_dir", type=Path,
                        help="Processed corpus output dir to serve")
    parser.add_argument("--host", default="127.0.0.1",
                        help="Bind interface (default: 127.0.0.1)")
    parser.add_argument("--port", type=int, default=18080,
                        help="Bind port (default: 18080, to avoid clashing "
                             "with a real 8080 service you may be running)")
    parser.add_argument("--server-python", default=sys.executable,
                        help="Python executable for the mcp_server subprocess")
    parser.add_argument("--skip-layer2", action="store_true",
                        help="Skip the MCP-client layer (useful when the "
                             "mcp SDK isn't installed in the current env)")
    parser.add_argument("--skip-layer3", action="store_true",
                        help="Skip the broad tool-coverage layer")
    parser.add_argument("--server-startup-timeout", type=float, default=15.0,
                        help="Seconds to wait for the server to bind (default: 15)")
    args = parser.parse_args()

    token = secrets.token_hex(16)
    tok_file = tempfile.NamedTemporaryFile(
        mode="w", delete=False, suffix=".token"
    )
    try:
        tok_file.write(token)
        tok_file.close()

        stderr_log = tempfile.NamedTemporaryFile(
            mode="w+", delete=False, suffix=".stderr"
        )
        stderr_log.close()

        cmd = [
            args.server_python, "-m", SERVER_MODULE, str(args.output_dir),
            "--transport", "sse", "--host", args.host, "--port", str(args.port),
            "--auth-token-file", tok_file.name,
        ]
        _section(f"Launching `python -m {SERVER_MODULE}` --transport sse")
        _info(f"cmd: {' '.join(cmd)}")
        _info(f"stderr -> {stderr_log.name}")

        # stderr goes to a file so post-mortem diagnostics are readable
        # without blocking on a pipe.
        with open(stderr_log.name, "w") as errfh:
            server = subprocess.Popen(
                cmd,
                stdout=subprocess.DEVNULL,
                stderr=errfh,
            )

        rc = 0
        try:
            if not _wait_for_bind(args.host, args.port,
                                  args.server_startup_timeout):
                _fail(
                    f"server did not bind on {args.host}:{args.port} "
                    f"within {args.server_startup_timeout}s"
                )
                tail = Path(stderr_log.name).read_text()[-2000:]
                print(f"\nstderr tail:\n{tail}")
                return 1
            _ok(f"server bound on {args.host}:{args.port}")

            rc |= layer1_http(args.host, args.port, token)
            if not args.skip_layer2:
                rc |= asyncio.run(
                    layer2_mcp_client(args.host, args.port, token)
                )
            if not args.skip_layer3:
                rc |= asyncio.run(
                    layer3_tool_coverage(args.host, args.port, token)
                )

        finally:
            server.terminate()
            try:
                server.wait(timeout=5)
            except subprocess.TimeoutExpired:
                server.kill()
                server.wait()

        _section("Result")
        if rc == 0:
            print(f"  {_color('All layers passed.', '32')}")
        else:
            print(f"  {_color('Some layers FAILED.', '31')}")
            print(f"  stderr log: {stderr_log.name}")
        return rc
    finally:
        os.unlink(tok_file.name)


if __name__ == "__main__":
    sys.exit(main())
