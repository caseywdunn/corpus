"""Read-only gold build -> vector rows -> live stdio MCP acceptance.

Prepared only: execute after the coordinator confirms embedding/post/bundle
completion and provides frozen current source-review expectations. No query
embedding, model inference, PDF preparation or extraction is performed here.
"""
from __future__ import annotations

import argparse
import asyncio
import base64
from collections import Counter, defaultdict
from datetime import datetime, timezone
import hashlib
import io
import json
import math
import os
from pathlib import Path
import re
import subprocess
import sys
import traceback


def sha(path):
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def encoded(value):
    return json.dumps(value, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode()


def digest(value):
    return hashlib.sha256(encoded(value)).hexdigest()


def load(path):
    return json.loads(Path(path).read_text())


def dump(path, value):
    Path(path).write_text(json.dumps(value, indent=2, ensure_ascii=False) + "\n")


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def inventory(root):
    return {str(p.relative_to(root)): sha(p) for p in sorted(root.rglob("*")) if p.is_file()}


def git(repo, *args):
    return subprocess.check_output(["git", *args], cwd=repo, text=True).strip()


def decode(response, *, list_result=False):
    require(not response.is_error, f"MCP error: {response}")
    if list_result:
        structured = response.structured_content
        if isinstance(structured, dict) and isinstance(structured.get("result"), list):
            return structured["result"]
    values = [json.loads(block.text) for block in response.content if block.type == "text"]
    if list_result:
        return values[0] if len(values) == 1 and isinstance(values[0], list) else values
    return values[0] if len(values) == 1 else values


def batches(chunks, count, max_bytes):
    """Explicit ID batches; reserve each row's maximum context projection twice
    for SDK text/structured duplication plus JSON envelope overhead.
    """
    batch, size = [], 4096
    for chunk in chunks:
        row_bytes = 2 * (len(encoded({key: chunk.get(key) for key in
                                     ("chunk_id", "text", "headings", "section_class", "figure_refs")}))
                         + 8192 + 1024)
        require(row_bytes + 4096 <= max_bytes,
                f"One chunk exceeds the verifier's response budget: {chunk['chunk_id']}")
        if batch and (len(batch) == count or size + row_bytes > max_bytes
                      or sum(len(c["chunk_id"]) for c in batch) + len(chunk["chunk_id"]) > 65536):
            yield batch
            batch, size = [], 4096
        batch.append(chunk)
        size += row_bytes
    if batch or not chunks:
        yield batch


def preparation(args, receipt):
    from pipeline.embed import ChunkMetadata
    from pipeline.embedding_state import load_document_data, record_payloads

    build, bundle = args.gold / "output", args.gold / "output/corpus_bundle"
    input_receipt = load(args.gold / "input-receipt.json")
    all_inputs = input_receipt["inputs"]
    verified_inputs = []
    for entry in all_inputs:
        actual = sha(args.gold / entry["path"])
        require(actual == entry["sha256"], f"Input changed: {entry['path']}")
        verified_inputs.append({"path": entry["path"], "sha256": actual,
                                "kind": "pdf_document" if Path(entry["path"]).suffix.lower() == ".pdf"
                                else "supporting_input"})
    source_rows = [r for r in all_inputs if Path(r["path"]).suffix.lower() == ".pdf"]
    sources = {r["sha256"][:12]: r for r in source_rows}
    require(len(sources) == len(source_rows), "Duplicate source identities")
    require(len(sources) == input_receipt["gold_document_count"], "Input receipt count mismatch")
    require(set(sources) == {p.name for p in (bundle / "documents").iterdir() if p.is_dir()},
            "Bundle membership differs from the frozen gold inputs")
    require(set(sources) == {p.name for p in (build / "documents").iterdir() if p.is_dir()},
            "Build membership differs from the frozen gold inputs")
    receipt["source_input_receipt_sha256"] = sha(args.gold / "input-receipt.json")
    receipt["source_input_receipt"] = input_receipt
    receipt["verified_inputs"] = verified_inputs
    receipt["input_count"] = len(all_inputs)
    receipt["supporting_input_count"] = len(all_inputs) - len(source_rows)
    receipt["bundle_manifest"] = load(bundle / "bundle_manifest.json")
    receipt["bundle_manifest_sha256"] = sha(bundle / "bundle_manifest.json")
    receipt["embedding_producer"] = load(bundle / "embedding_producer.json")
    receipt["embedding_producer_sha256"] = sha(bundle / "embedding_producer.json")
    controls = args.out / "gold-control-receipts"
    controls.mkdir()
    receipt["gold_control_receipts"] = {}
    for path in sorted(args.gold.glob("*.json")):
        if path.name == "input-receipt.json":
            continue
        copied = controls / path.name
        copied.write_bytes(path.read_bytes())
        receipt["gold_control_receipts"][path.name] = {
            "sha256": sha(copied), "artifact": str(copied.relative_to(args.out))}
    expected, chunksets, documents = {}, {}, {}
    for paper, source in sorted(sources.items()):
        directory = build / "documents" / paper
        served_dir = bundle / "documents" / paper
        doc = load_document_data(directory)
        chunks = load(served_dir / "chunks.json")["chunks"]
        require(chunks == doc["chunks"]["chunks"], f"Bundle chunks differ from build: {paper}")
        require(len({c["chunk_id"] for c in chunks}) == len(chunks), f"Duplicate chunks: {paper}")
        chunksets[paper] = chunks
        for row in record_payloads(paper, doc):
            # Stage 2 validates through its Pydantic metadata model before Arrow
            # storage. Normalize only expected inputs through that same model;
            # do not silently replace nulls/defaults in the observed index rows.
            row["metadata"] = ChunkMetadata.model_validate(row["metadata"]).model_dump()
            key = paper, row["metadata"]["chunk_id"]
            require(key not in expected, f"Duplicate expected vector key: {key}")
            expected[key] = row
        names = ("text.json", "chunks.json", "docling_doc.json", "figures.json", "processed.pdf",
                 "pipeline_state.json", "scan_detection.json", "summary.json", "metadata.json")
        documents[paper] = {"chunk_count": len(chunks), "source_sha256": source["sha256"],
                            "build_artifact_sha256": {n: sha(directory / n) for n in names
                                                      if (directory / n).is_file()}}
    receipt["documents"] = documents
    receipt["paper_count"] = len(sources)
    receipt["expected_chunk_count"] = len(expected)
    require(receipt["bundle_manifest"]["paper_count"] == len(sources), "Manifest paper count mismatch")
    require(receipt["bundle_manifest"]["chunk_count"] == len(expected), "Manifest chunk count mismatch")
    return build, bundle, sources, expected, chunksets


def verify_vectors(args, build, bundle, expected, receipt):
    import lancedb
    from pipeline.embedding_state import marker_problem, read_marker
    from pipeline.model_provenance import same_embedding_space

    database = lancedb.connect(str(bundle / "vector_db/lancedb"))
    table = database.open_table("document_chunks")
    before = table.version
    require(table.schema.field("vector").type.list_size == 1024, "Index dimension is not 1024")
    seen, census, bundled_rows = set(), defaultdict(Counter), {}
    with (args.out / "vector-row-checks.jsonl").open("w") as raw:
        for batch in table.search().select(["text", "metadata", "vector", "embedding_generation"]).limit(None).to_batches():
            for row in batch.to_pylist():
                key = row["metadata"]["pdf_hash"], row["metadata"]["chunk_id"]
                require(key in expected and key not in seen, f"Orphan or duplicate vector row: {key}")
                seen.add(key)
                require(row["text"] == expected[key]["text"], f"Vector text differs: {key}")
                require(row["metadata"] == expected[key]["metadata"], f"Vector metadata/headings differ: {key}")
                vector = row["vector"]
                require(len(vector) == 1024 and all(math.isfinite(v) for v in vector)
                        and any(v != 0 for v in vector), f"Invalid/zero vector: {key}")
                census[key[0]][row["embedding_generation"]] += 1
                bundled_rows[key] = row
                raw.write(json.dumps({"paper_hash": key[0], "chunk_id": key[1],
                    "text_sha256": digest(row["text"]), "metadata_sha256": digest(row["metadata"]),
                    "vector_sha256": digest(vector), "vector_dim": len(vector),
                    "vector_l2_norm": math.sqrt(sum(v * v for v in vector)),
                    "embedding_generation": row["embedding_generation"]}, ensure_ascii=False) + "\n")
    require(seen == set(expected), f"Vector rows missing: {sorted(set(expected) - seen)[:12]}")
    identity = receipt["embedding_producer"]
    require(identity["dimension"] == 1024 and receipt["bundle_manifest"]["embedding_dim"] == 1024,
            "Producer/manifest dimension mismatch")
    markers = {}
    for paper in sorted(receipt["documents"]):
        marker_path = build / "vector_db" / f"{paper}_embedded.done"
        marker = read_marker(marker_path)
        problem = marker_problem(build / "documents" / paper, marker, census,
                                 dim=1024, producer=identity["producer"])
        require(problem is None, f"Embedding completion evidence failed for {paper}: {problem}")
        require(same_embedding_space(marker["embedding_producer"], identity["producer"]),
                f"Embedding space differs: {paper}")
        markers[paper] = {"sha256": sha(marker_path), "marker": marker}
    dump(args.out / "embedding-markers.json", markers)
    after = database.open_table("document_chunks").version
    require(before == after, "Vector index changed during verification")
    # The published bundle precedes the unchanged-resume check. Compare every
    # current build row directly to that immutable snapshot, including vectors;
    # transaction manifests/version numbers themselves need not be equal.
    build_database = lancedb.connect(str(build / "vector_db/lancedb"))
    build_table = build_database.open_table("document_chunks")
    build_before = build_table.version
    compared = set()
    with (args.out / "resume-row-comparison.jsonl").open("w") as raw:
        for batch in build_table.search().select(["text", "metadata", "vector", "embedding_generation"]).limit(None).to_batches():
            for row in batch.to_pylist():
                key = row["metadata"]["pdf_hash"], row["metadata"]["chunk_id"]
                require(key in bundled_rows and key not in compared, f"Resumed build orphan/duplicate row: {key}")
                reference = bundled_rows[key]
                differences = [field for field in reference if row.get(field) != reference[field]]
                require(not differences, f"Resumed build differs from published snapshot: {key}: {differences}")
                compared.add(key)
                raw.write(json.dumps({"paper_hash": key[0], "chunk_id": key[1],
                                      "exact_equal_fields": sorted(reference),
                                      "canonical_row_sha256": digest(row)}, ensure_ascii=False) + "\n")
    require(compared == seen, f"Resumed build missing rows: {sorted(seen - compared)[:12]}")
    build_after = build_database.open_table("document_chunks").version
    require(build_before == build_after, "Build vector table changed during read-only comparison")
    strict_before_path = args.gold / "refresh-20260922-resume-before.json"
    strict_after_path = args.gold / "refresh-20260922-resume-after.json"
    strict_before, strict_after = load(strict_before_path), load(strict_after_path)
    changed = sorted(key for key in set(strict_before["files"]) | set(strict_after["files"])
                     if strict_before["files"].get(key) != strict_after["files"].get(key))
    prefix = "vector_db/lancedb/document_chunks.lance/"
    bookkeeping_only = all(path.startswith(prefix) and (
        path[len(prefix):].startswith("_transactions/") and path.endswith(".txn")
        or path[len(prefix):].startswith("_versions/") and
        (path.endswith(".manifest") or path.endswith("/latest_version_hint.json"))) for path in changed)
    require(bookkeeping_only, f"Unchanged resume modified non-bookkeeping files: {changed}")
    receipt["unchanged_resume_semantic_comparison"] = {
        "status": "pass", "rows_directly_compared": len(compared),
        "all_keys_text_typed_metadata_vectors_and_generations_equal": True,
        "bundle_table_version_before": before, "bundle_table_version_after": after,
        "resumed_build_table_version_before": build_before, "resumed_build_table_version_after": build_after,
        "strict_file_hash_resume_unchanged": strict_after["unchanged"],
        "strict_before_receipt_sha256": sha(strict_before_path),
        "strict_after_receipt_sha256": sha(strict_after_path),
        "strict_receipt_preserved_without_changes": True,
        "changed_file_paths": changed, "changed_files_are_lance_transaction_bookkeeping_only": bookkeeping_only,
        "interpretation": "Semantic row equality is proved independently. Strict file-byte resume equality remains false when Lance commits bookkeeping; no document, marker or vector data-file difference is suppressed.",
        "row_comparison_sha256": sha(args.out / "resume-row-comparison.jsonl")}
    receipt["vector_verification"] = {"rows": len(seen), "dimension": 1024,
        "table_version_before": before, "table_version_after": after,
        "all_rows_unique": True, "exact_build_payloads": True, "all_finite_and_nonzero": True,
        "all_current_generation_and_producer_markers_valid": True,
        "row_checks_sha256": sha(args.out / "vector-row-checks.jsonl"),
        "embedding_markers_sha256": sha(args.out / "embedding-markers.json")}


def frozen_expectations(args, sources):
    data = load(args.expectations)
    require(data.get("status") == "frozen" and data.get("expectations"),
            "Scientific expectations must be frozen and nonempty before serving")
    require(data.get("review_results") and data.get("source_labels"), "Review provenance missing")
    for field in ("review_results", "source_labels"):
        require(sha(data[field]["path"]) == data[field]["sha256"], f"Frozen {field} changed")
    for entry in data["expectations"]:
        paper = entry["paper_hash"]
        require(paper in sources and entry["source_sha256"] == sources[paper]["sha256"],
                f"Expectation source is outside gold inputs: {paper}")
        require(entry.get("source_review_id") and entry.get("source_item_refs"),
                "Expectations require the source-reviewed occurrence, not any matching paper text")
        require(bool(entry.get("literal")) != bool(entry.get("pattern")),
                "Specify exactly one literal or declared pattern per expectation")
    return data


async def verify_mcp(args, bundle, sources, chunksets, expectations, receipt):
    from mcp import ClientSession
    from mcp.client.stdio import StdioServerParameters, stdio_client
    from PIL import Image
    from pipeline.figures import EVIDENCE_FIGURE_TYPES

    params = StdioServerParameters(command=sys.executable,
        args=["-m", "mcpsrv.main", str(bundle), "--transport", "stdio"], cwd=str(args.repo),
        env={**os.environ, "HF_HUB_OFFLINE": "1", "TRANSFORMERS_OFFLINE": "1",
             "PYTHONDONTWRITEBYTECODE": "1"})
    calls, served_by_paper, figures = [], {}, []
    with (args.out / "server.stderr").open("w") as stderr:
        async with stdio_client(params, errlog=stderr) as streams:
            async with ClientSession(*streams, read_timeout_seconds=60) as session:
                initialized = await session.initialize()
                dump(args.out / "initialize.json", initialized.model_dump(mode="json", by_alias=True))
                for paper, chunks in sorted(chunksets.items()):
                    rows = []
                    for number, batch in enumerate(batches(chunks, args.batch_size, args.response_max_bytes)):
                        call = {"paper_hash": paper, "chunk_ids": [c["chunk_id"] for c in batch], "with_text": True}
                        response = await session.call_tool("get_chunks", call)
                        raw = response.model_dump(mode="json", by_alias=True)
                        filename = f"chunks-{paper}-{number:04d}.json"
                        dump(args.out / filename, {"call": call, "response": raw})
                        raw_bytes = len(encoded(raw))
                        require(raw_bytes <= args.response_max_bytes, f"MCP response exceeds byte budget: {filename}")
                        values = decode(response, list_result=True)
                        require(isinstance(values, list), f"Unexpected get_chunks shape: {filename}")
                        require([r.get("chunk_id") for r in values] == call["chunk_ids"],
                                f"MCP missing/duplicate/reordered chunks: {filename}")
                        for value, chunk in zip(values, batch):
                            require("error" not in value and value["text"] == (chunk.get("text") or "")
                                    and value["headings"] == (chunk.get("headings") or []),
                                    f"MCP text/headings differ: {paper}/{chunk['chunk_id']}")
                        rows.extend(values)
                        calls.append({"paper_hash": paper, "batch": number, "rows": len(values),
                                      "response_bytes": raw_bytes, "artifact": filename,
                                      "sha256": sha(args.out / filename)})
                    require(len(rows) == len(chunks), f"Incomplete MCP paper: {paper}")
                    served_by_paper[paper] = rows
                # Deterministic bounded sample: first materialized figure per paper,
                # lexical paper order, skipping absent files/oversized image payloads.
                for paper in sorted(sources):
                    artifact = bundle / "documents" / paper / "figures.json"
                    candidates = load(artifact).get("figures", []) if artifact.is_file() else []
                    candidates = sorted((f for f in candidates if f.get("figure_id") and f.get("filename")
                                         and f.get("figure_type") in EVIDENCE_FIGURE_TYPES
                                         and f.get("caption_status") == "bound"), key=lambda f: f["figure_id"])
                    selected = next((f for f in candidates
                        if (bundle / "documents" / paper / "figures" / f["filename"]).is_file()
                        and (bundle / "documents" / paper / "figures" / f["filename"]).stat().st_size
                        <= args.image_max_bytes), None)
                    if selected is None:
                        continue
                    call = {"paper_hash": paper, "figure_id": selected["figure_id"]}
                    metadata_response = await session.call_tool("get_figure", {**call, "include_licensing": True})
                    require(len(encoded(metadata_response.model_dump(mode="json", by_alias=True)))
                            <= args.response_max_bytes, f"Figure metadata exceeds response budget: {call}")
                    metadata = decode(metadata_response)
                    require("error" not in metadata and metadata["caption_text"] == selected.get("caption_text")
                            and metadata["caption_status"] == selected["caption_status"], f"Figure metadata differs: {call}")
                    name = f"figure-{paper}-{len(figures):02d}"
                    dump(args.out / f"{name}-metadata.json", metadata_response.model_dump(mode="json", by_alias=True))
                    checked = {**call, "publication_clearance": metadata.get("publication_clearance"),
                               "metadata_sha256": sha(args.out / f"{name}-metadata.json"), "profiles": {}}
                    for profile in ("manuscript", "report"):
                        response = await session.call_tool("get_figure_image", {**call, "profile": profile})
                        dump(args.out / f"{name}-{profile}.json", response.model_dump(mode="json", by_alias=True))
                        images = [b for b in response.content if b.type == "image"]
                        strict_permitted = metadata.get("publication_clearance") in {"public_domain", "licensed_open"}
                        if profile == "manuscript" and not strict_permitted:
                            require(response.is_error and not images, f"Strict figure delivered without clearance: {call}")
                            denial = response.structured_content
                            require(denial and denial.get("code") == "forbidden" and denial.get("profile") == profile,
                                    f"Unstructured strict denial: {call}")
                            checked["profiles"][profile] = {"status": "refused", "denial": denial}
                            continue
                        require(not response.is_error and len(images) == 1, f"Permitted figure delivery failed: {call}")
                        image_bytes = base64.b64decode(images[0].data, validate=True)
                        expected = bundle / "documents" / paper / "figures" / selected["filename"]
                        require(hashlib.sha256(image_bytes).hexdigest() == sha(expected), f"Image bytes differ: {call}")
                        image = Image.open(io.BytesIO(image_bytes))
                        image.load()
                        checked["profiles"][profile] = {"status": "delivered", "sha256": sha(expected),
                            "bytes": len(image_bytes), "size": list(image.size), "equal_to_bundle": True}
                    figures.append(checked)
                    if len(figures) >= args.max_figures:
                        break
    checks = []
    for entry in expectations["expectations"]:
        chunks = {c["chunk_id"]: c for c in chunksets[entry["paper_hash"]]}
        matched = []
        for row in served_by_paper[entry["paper_hash"]]:
            refs = {p.get("item_ref") for p in chunks[row["chunk_id"]].get("source_items", [])}
            if not refs.intersection(entry["source_item_refs"]):
                continue
            found = entry["literal"] in row["text"] if entry.get("literal") else re.search(entry["pattern"], row["text"])
            if found:
                matched.append(row["chunk_id"])
        require(matched, f"Reviewed scientific expression absent from its source occurrence: {entry}")
        checks.append({**entry, "served_chunk_ids": matched, "exact_vector_and_bundle_text_verified": True})
    receipt["mcp"] = {"transport": "stdio", "papers": len(served_by_paper), "chunk_calls": len(calls),
        "rows": sum(len(rows) for rows in served_by_paper.values()), "all_rows_exact": True,
        "batch_size_cap": args.batch_size, "response_max_bytes": args.response_max_bytes,
        "largest_response_bytes": max(c["response_bytes"] for c in calls), "calls": calls,
        "figures": figures, "figure_selection": "first bound materialized image per paper, lexical order, size cap",
        "scientific_expression_checks": checks}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gold", type=Path, default=Path("/tmp/corpus-v15-acceptance/gold"))
    parser.add_argument("--repo", type=Path, default=Path("/tmp/corpus-correctness-validation"))
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--expectations", type=Path, required=True)
    parser.add_argument("--batch-size", type=int, default=25)
    parser.add_argument("--response-max-bytes", type=int, default=524288)
    parser.add_argument("--max-figures", type=int, default=6)
    parser.add_argument("--image-max-bytes", type=int, default=8388608)
    args = parser.parse_args()
    args.gold, args.repo, args.out = args.gold.resolve(), args.repo.resolve(), args.out.resolve()
    require(args.gold not in args.out.parents and args.out != args.gold, "Evidence must be outside gold")
    require(1 <= args.batch_size <= 500 and 1 <= args.max_figures <= 20, "Invalid bounded call settings")
    require(not args.out.exists(), "Use a fresh output directory; never overwrite prior acceptance evidence")
    args.out.mkdir(parents=True)
    sys.path.insert(0, str(args.repo))
    bundle = args.gold / "output/corpus_bundle"
    receipt = {"status": "running", "started_at": datetime.now(timezone.utc).isoformat(),
        "script_sha256": sha(__file__), "server_code_revision": git(args.repo, "rev-parse", "HEAD"),
        "server_worktree_changes": git(args.repo, "status", "--short"),
        "scope": "Complete frozen gold document/chunk serving and vector integrity, selected figure deliveries; no retrieval ranking or scientific output precision claim.",
        "expected_limitations": ["No model inference or semantic query", "Finite nonzero vectors plus producer receipts do not independently prove semantic quality",
            "Only source-reviewed present scientific occurrences are checked; omitted/corrupted/unscorable review cases remain in their source review denominator",
            "Figure selection is a bounded contract replay, not exhaustive caption/rights/source-image acceptance",
            "Report-profile delivery does not establish republication permission"]}
    before = inventory(bundle)
    dump(args.out / "bundle-inventory-before.json", before)
    try:
        build, bundle, sources, expected, chunksets = preparation(args, receipt)
        expectations = frozen_expectations(args, sources)
        dump(args.out / "expectations.frozen.json", expectations)
        receipt["scientific_expectations_sha256"] = sha(args.expectations)
        verify_vectors(args, build, bundle, expected, receipt)
        asyncio.run(verify_mcp(args, bundle, sources, chunksets, expectations, receipt))
        receipt["status"] = "pass"
    except BaseException as exc:
        receipt.update(status="fail", error=str(exc), exception_type=type(exc).__name__)
        (args.out / "failure.txt").write_text(traceback.format_exc())
    finally:
        after = inventory(bundle)
        dump(args.out / "bundle-inventory-after.json", after)
        receipt["bundle_inventory_before_sha256"] = digest(before)
        receipt["bundle_inventory_after_sha256"] = digest(after)
        receipt["bundle_unchanged"] = before == after
        if before != after:
            receipt["status"] = "fail"
            receipt["changed_bundle_paths"] = [key for key in sorted(set(before) | set(after)) if before.get(key) != after.get(key)]
        receipt["completed_at"] = datetime.now(timezone.utc).isoformat()
        receipt["server_code_revision_after"] = git(args.repo, "rev-parse", "HEAD")
        if receipt["server_code_revision_after"] != receipt["server_code_revision"]:
            receipt.update(status="fail", code_changed_during_verification=True)
        dump(args.out / "serving-receipt.json", receipt)
    print(json.dumps({key: receipt.get(key) for key in ("status", "paper_count", "expected_chunk_count", "bundle_unchanged", "error")}))
    return 0 if receipt["status"] == "pass" else 1


if __name__ == "__main__":
    raise SystemExit(main())
