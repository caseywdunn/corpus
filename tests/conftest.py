"""Fixtures for corpus quality tests.

Ground truth lives in tests/ground_truth/*.yaml (one file per paper).
Tests are parametrized over these files automatically.
"""

import hashlib
import json
import os
from pathlib import Path

import pytest
import yaml

REPO_ROOT = Path(__file__).parent.parent
DEMO_DIR = REPO_ROOT / "demo"
GROUND_TRUTH_DIR = Path(__file__).parent / "ground_truth"

# Default output location — override with CORPUS_OUTPUT_DIR env var
_DEFAULT_OUTPUT = Path("/nfs/roberts/project/pi_cwd7/cwd7/output")


def _output_dir():
    override = os.environ.get("CORPUS_OUTPUT_DIR")
    if override:
        return Path(override)
    repo_output = REPO_ROOT / "output"
    if repo_output.is_dir():
        return repo_output
    return _DEFAULT_OUTPUT


_HASH_MAP = None


def _build_hash_map():
    """Map each demo PDF filename stem to its 12-char SHA-256 prefix."""
    global _HASH_MAP
    if _HASH_MAP is not None:
        return _HASH_MAP
    _HASH_MAP = {}
    for pdf in sorted(DEMO_DIR.glob("*.pdf")):
        h = hashlib.sha256(pdf.read_bytes()).hexdigest()[:12]
        _HASH_MAP[pdf.stem] = h
    return _HASH_MAP


def _load_ground_truth_files():
    """Load all YAML ground truth files, return list of (stem, data) tuples."""
    results = []
    for path in sorted(GROUND_TRUTH_DIR.glob("*.yaml")):
        with open(path) as f:
            data = yaml.safe_load(f)
        results.append((path.stem, data))
    return results


def _paper_ids():
    return [stem for stem, _ in _load_ground_truth_files()]


def _load_paper(stem):
    gt_file = GROUND_TRUTH_DIR / f"{stem}.yaml"
    with open(gt_file) as f:
        gt = yaml.safe_load(f)
    pdf_hash = gt.get("hash") or _build_hash_map().get(stem)
    doc_dir = _output_dir() / "documents" / pdf_hash
    return gt, doc_dir


def _load_json(doc_dir, filename):
    path = doc_dir / filename
    if not path.exists():
        pytest.skip(f"{path} not found")
    with open(path) as f:
        return json.load(f)


# ---------------------------------------------------------------------------
# Parametrized fixtures — one test invocation per ground-truth file
# ---------------------------------------------------------------------------

@pytest.fixture(params=_paper_ids())
def paper(request):
    """Yields (ground_truth_dict, doc_dir_path) for each ground-truth file."""
    return _load_paper(request.param)


@pytest.fixture(params=_paper_ids())
def paper_text(request):
    """Yields (ground_truth_dict, text_data_dict)."""
    gt, doc_dir = _load_paper(request.param)
    if "text" not in gt:
        pytest.skip("No text ground truth")
    return gt, _load_json(doc_dir, "text.json")


@pytest.fixture(params=_paper_ids())
def paper_figures(request):
    """Yields (ground_truth_dict, figures_data_dict)."""
    gt, doc_dir = _load_paper(request.param)
    if "figures" not in gt:
        pytest.skip("No figures ground truth")
    return gt, _load_json(doc_dir, "figures.json"), doc_dir


@pytest.fixture(params=_paper_ids())
def paper_references(request):
    """Yields (ground_truth_dict, references_data_dict)."""
    gt, doc_dir = _load_paper(request.param)
    if "references" not in gt:
        pytest.skip("No references ground truth")
    return gt, _load_json(doc_dir, "references.json")


@pytest.fixture(params=_paper_ids())
def paper_metadata(request):
    """Yields (ground_truth_dict, metadata_dict)."""
    gt, doc_dir = _load_paper(request.param)
    if "metadata" not in gt:
        pytest.skip("No metadata ground truth")
    return gt, _load_json(doc_dir, "metadata.json")
