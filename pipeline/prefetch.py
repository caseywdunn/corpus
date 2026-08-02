"""Model prefetch + cache inspection (#159).

A first ``corpus run`` on a fresh install has to reach HuggingFace for
four models before it can do any work: docling's page-layout model,
docling's TableFormer (``do_table_structure=True``), docling's
``HybridChunker`` tokenizer, and the embedding model (BGE-M3 by
default). Until this module existed there was no way to
fetch them ahead of time, no way to ask whether they were already
cached, and no documentation of where they land — so an offline or
rate-limited host discovered the problem deep inside the extract or
embed stage, after `corpus check` had already reported green.

Two entry points:

* :func:`model_report` — read-only, **never touches the network**. Feeds
  ``corpus check`` so a missing model is a pre-flight line rather than a
  mid-run traceback.
* :func:`prefetch_all` — does the downloads, with backoff, by exercising
  the real code paths rather than guessing repo ids. Backs
  ``corpus prefetch``.

Why prefetch by doing real work: the exact repo ids a docling release
pulls are a private implementation detail that changes between versions
(2.94.0 fetches ``docling-layout-heron`` + ``docling-models``; earlier
versions differed). Running a genuine conversion downloads whatever the
*pinned* version needs, with no list for us to keep in sync. The repo
ids below are therefore used only for **reporting**, and a repo we don't
recognise is never an error.

Operators who need a specific cache location set ``HF_HOME`` (or the
narrower ``HF_HUB_CACHE``) before running — see INSTALL.md. If your
cluster's compute nodes are network-restricted, run ``corpus prefetch``
on a host that does have outbound access, with ``HF_HOME`` pointed at
shared storage, then run the pipeline with ``HF_HUB_OFFLINE=1`` so an
accidental fetch fails loudly instead of hanging. ``HF_HUB_OFFLINE=1`` is
worth setting even on an unrestricted cluster: it pins the run to the
snapshot you cached rather than whatever ``main`` points at today.
"""
from __future__ import annotations

import logging
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple

logger = logging.getLogger(__name__)


# Repo ids the pinned stack is expected to use. Reporting only — see the
# module docstring. Keep in sync with environment.yaml's docling pin on a
# best-effort basis; being wrong here degrades a `corpus check` line to
# "unknown", it does not break a run.
#
# The chunker tokenizer is listed for a second reason beyond reporting:
# `corpus prefetch` short-circuits ("nothing to do") when every repo here
# is cached, so a model absent from this dict is never fetched at all.
# Leaving it out is what let a "fully prefetched" host still fall back to
# the naive chunker under HF_HUB_OFFLINE=1.
CHUNKER_TOKENIZER_REPO = "sentence-transformers/all-MiniLM-L6-v2"

DOCLING_REPOS: Dict[str, str] = {
    "docling-project/docling-layout-heron": "docling page-layout model",
    "docling-project/docling-models": "docling TableFormer (table structure)",
    CHUNKER_TOKENIZER_REPO: "docling HybridChunker tokenizer (text chunking)",
}

# The embedding default, mirroring pipeline.embeddings._DEFAULT_LOCAL_MODEL.
# Imported lazily in :func:`_embedding_repo` so this module stays cheap to
# import (sentence_transformers pulls torch).
_EMBED_FALLBACK = "BAAI/bge-m3"

# The local VLM is ~16 GB and only needed for `--figure-panels
# vision-local`, so it is opt-in on both entry points.
VISION_REPO = "Qwen/Qwen2.5-VL-7B-Instruct"


def _embedding_repo(model: Optional[str] = None) -> str:
    if model:
        return model
    try:
        from .embeddings import _DEFAULT_LOCAL_MODEL
        return _DEFAULT_LOCAL_MODEL
    except Exception:  # pragma: no cover - defensive
        return _EMBED_FALLBACK


def hf_cache_dir() -> Optional[Path]:
    """Where HuggingFace will read/write models, or None if huggingface_hub
    isn't importable."""
    try:
        from huggingface_hub import constants
        return Path(constants.HF_HUB_CACHE)
    except Exception:
        return None


def cached_repos() -> Dict[str, int]:
    """``{repo_id: size_on_disk_bytes}`` for every repo in the local HF
    cache. Offline and side-effect free; returns ``{}`` when the cache is
    absent or unreadable (a fresh install, typically)."""
    try:
        from huggingface_hub import scan_cache_dir
        info = scan_cache_dir()
    except Exception as e:
        logger.debug("HF cache scan failed (treating as empty): %s", e)
        return {}
    return {r.repo_id: r.size_on_disk for r in info.repos}


def required_repos(
    *, embedding_model: Optional[str] = None, include_vision: bool = False,
) -> Dict[str, str]:
    """``{repo_id: human label}`` for the models a run will need."""
    repos = dict(DOCLING_REPOS)
    repos[_embedding_repo(embedding_model)] = "embedding model (semantic search)"
    if include_vision:
        repos[VISION_REPO] = "local vision-language model (--figure-panels vision-local)"
    return repos


def model_report(
    *, embedding_model: Optional[str] = None, include_vision: bool = False,
) -> Tuple[List[Dict], bool]:
    """``(rows, all_present)`` describing local model-cache state.

    Each row is ``{repo, label, cached, size_bytes}``. Never hits the
    network, so this is safe on an air-gapped host and safe to call from
    ``corpus check``.
    """
    have = cached_repos()
    rows: List[Dict] = []
    for repo, label in required_repos(
        embedding_model=embedding_model, include_vision=include_vision
    ).items():
        rows.append({
            "repo": repo,
            "label": label,
            "cached": repo in have,
            "size_bytes": have.get(repo, 0),
        })
    return rows, all(r["cached"] for r in rows)


def human_size(n: int) -> str:
    f = float(n)
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if f < 1024 or unit == "TB":
            return f"{f:.0f} {unit}" if unit in ("B", "KB") else f"{f:.1f} {unit}"
        f /= 1024
    return f"{f:.1f} TB"  # pragma: no cover


# ---------------------------------------------------------------------------
# Downloads
# ---------------------------------------------------------------------------
#
# HuggingFace throttles anonymous traffic and returns 429. That is not
# hypothetical: it took CI down (#140) when runner IPs got throttled, and
# it is exactly what a shared institutional NAT looks like from the other
# side. So every fetch retries with linear backoff instead of failing on
# the first refusal.

_ATTEMPTS = 6
_BACKOFF_STEP_S = 20


# Failures that retrying cannot fix — they mean the code is wrong, not the
# network. Retrying one of these costs the user ~5 minutes of backoff and
# then reports "failed after 6 attempts", which reads like a flaky network
# and buries the real cause. (Observed: a wrong method name on the
# embedding backend survived six attempts, after the 4.8 GB download had
# already succeeded.)
_PROGRAMMING_ERRORS = (AttributeError, TypeError, NameError, ImportError)


def _with_retry(what: str, fn, attempts: int = _ATTEMPTS) -> None:
    last: Optional[Exception] = None
    for i in range(1, attempts + 1):
        try:
            fn()
            logger.info("%s: ok (attempt %d)", what, i)
            return
        except _PROGRAMMING_ERRORS as e:
            raise RuntimeError(
                f"{what}: {type(e).__name__}: {e}. This is a bug in corpus, "
                "not a network problem — retrying would not help. Please "
                "report it with this message."
            ) from e
        except Exception as e:  # noqa: BLE001 - anything else may be transient
            last = e
            logger.warning("%s: attempt %d/%d failed: %s", what, i, attempts, e)
            if i < attempts:
                time.sleep(i * _BACKOFF_STEP_S)
    raise RuntimeError(f"{what}: failed after {attempts} attempts: {last}")


def prefetch_docling(sample_pdf: Path, attempts: int = _ATTEMPTS) -> None:
    """Warm docling's models by converting *and chunking* one real PDF.

    Mirrors :func:`pipeline.extract` pipeline options so exactly the
    models a run needs are fetched — notably ``do_picture_classification``
    stays off (#140), so we no longer warm a classifier nothing reads.

    Chunking is exercised too, because ``HybridChunker`` pulls its own
    tokenizer (currently ``sentence-transformers/all-MiniLM-L6-v2``) the
    first time it runs — a model that conversion alone never touches.
    Missing it does not fail a run: :mod:`pipeline.chunking` catches the
    error and silently falls back to the naive character chunker, which
    turned a 2-page paper into 1 chunk instead of 16 on a host running
    the documented ``HF_HUB_OFFLINE=1`` recipe. Prefetching by doing the
    real work is the whole point of this module; chunking was the half
    it wasn't doing.
    """
    def _run() -> None:
        from docling.document_converter import DocumentConverter, PdfFormatOption
        from docling.datamodel.base_models import InputFormat
        from docling.datamodel.pipeline_options import PdfPipelineOptions
        opts = PdfPipelineOptions(
            do_ocr=False,
            do_table_structure=True,
            generate_picture_images=True,
            generate_page_images=False,
            do_picture_classification=False,
        )
        conv = DocumentConverter(
            format_options={InputFormat.PDF: PdfFormatOption(pipeline_options=opts)}
        )
        result = conv.convert(str(sample_pdf))

        # Same construction as pipeline.chunking, so whatever tokenizer
        # the pinned docling wants is the one that lands in the cache.
        from docling.chunking import HybridChunker
        list(HybridChunker().chunk(result.document))

    _with_retry("docling models", _run, attempts)


def prefetch_embedding(
    model: Optional[str] = None, attempts: int = _ATTEMPTS,
) -> None:
    """Warm the embedding model by constructing the real backend."""
    repo = _embedding_repo(model)

    def _run() -> None:
        from .embeddings import get_embedder
        emb = get_embedder(repo)
        # Touch the model so weights actually load rather than merely
        # resolving a config file. The method is ``embed`` — the
        # EmbeddingBackend ABC's own name, not sentence-transformers'
        # ``encode``. Calling the wrong one made `corpus prefetch` fail
        # *after* successfully downloading 4.8 GB, then burn six retries
        # on an AttributeError that no amount of retrying could fix.
        emb.embed(["warmup"])

    _with_retry(f"embedding model {repo}", _run, attempts)


def prefetch_vision(attempts: int = _ATTEMPTS) -> None:
    """Warm the local VLM (~16 GB). Opt-in only."""
    def _run() -> None:
        from huggingface_hub import snapshot_download
        snapshot_download(VISION_REPO)

    _with_retry(f"local VLM {VISION_REPO}", _run, attempts)


def prefetch_all(
    sample_pdf: Optional[Path] = None,
    *,
    embedding_model: Optional[str] = None,
    include_vision: bool = False,
    attempts: int = _ATTEMPTS,
) -> None:
    """Fetch every model a run needs. Raises ``RuntimeError`` on failure.

    ``sample_pdf`` is any PDF; docling needs a real document to exercise
    its model loads. The caller (``corpus prefetch``) passes one from the
    configured input directory, falling back to the packaged demo.
    """
    if sample_pdf is not None:
        prefetch_docling(Path(sample_pdf), attempts)
    else:
        logger.warning(
            "no sample PDF available — skipping docling model prefetch "
            "(pass a PDF, or run from a corpuscle with input_pdfs set)"
        )
    prefetch_embedding(embedding_model, attempts)
    if include_vision:
        prefetch_vision(attempts)
