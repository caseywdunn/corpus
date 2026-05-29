"""Pydantic v2 schema for the per-corpuscle config.yaml (#59).

Wraps the existing dict-based loader in :mod:`pipeline.config` with
field-level validation: missing required keys, wrong types, and
unknown enum values surface at startup with errors that point at
the exact key + value (e.g.,
``vision.backend must be one of {local, claude, none}, got 'claud'``)
rather than as ``KeyError`` halfway through Stage 1.

The dict-based ``CONFIG`` singleton in :mod:`pipeline.config` is left
in place so existing call sites (``CONFIG.get("ocr", {}).get(...)``)
keep working unchanged. The unified ``corpus`` CLI (#60) loads the
merged dict, validates it via :func:`validate_config`, and uses the
returned :class:`CorpuscleConfig` model for any new code paths;
legacy stage modules continue reading the dict.
"""
from __future__ import annotations

from pathlib import Path
from typing import List, Literal, Optional, Union

from pydantic import (
    BaseModel,
    ConfigDict,
    Field,
    ValidationError,
    field_validator,
    model_validator,
)


# ---------------------------------------------------------------------------
# Per-corpuscle inputs (new in v0.3 — replace today's CLI flags).
# ---------------------------------------------------------------------------


class TaxonomyConfig(BaseModel):
    """External taxonomy source (#59)."""

    model_config = ConfigDict(extra="forbid")

    source: Optional[Literal["worms", "dwc", "dwca"]] = Field(
        default=None,
        description="Taxonomy source kind. None disables auto-build.",
    )
    root_id: Optional[Union[int, str]] = Field(
        default=None,
        description="Root for the subtree to ingest. For source='worms' "
        "this is an integer WoRMS AphiaID (required). For "
        "source='dwc' or 'dwca' it's a taxonID value as it appears in "
        "the archive — often a bare integer, but DwC-A snapshots from "
        "WoRMS use LSIDs like 'urn:lsid:marinespecies.org:taxname:1371'. "
        "Optional for dwc/dwca; omit to ingest the full archive.",
    )
    path: Optional[Path] = Field(
        default=None,
        description="Local DwC-A archive or DwC directory "
        "(required when source in {'dwc', 'dwca'}).",
    )


class FiguresConfig(BaseModel):
    """Figure panel-ROI detection settings (#102, supersedes the v0.5
    ``vision`` block).

    ``panel_detection`` selects which pass annotates multi-panel figure
    ROIs. Pass 3a (``ocr``) is the cheap CPU-only floor and the default
    — it self-gates to figures whose caption implies multiple panels, so
    it's near-free on single-panel figures. The ``vision-*`` modes run
    Pass 3b instead (open-weights VLM or the Anthropic API); ``off``
    disables panel-ROI detection entirely.
    """

    model_config = ConfigDict(extra="forbid")

    panel_detection: Literal["ocr", "vision-local", "vision-claude", "off"] = Field(
        default="ocr",
        description="ocr = OCR-driven panel ROIs (CPU, default floor); "
        "vision-local = open-weights VLM on CUDA/MPS; "
        "vision-claude = Anthropic API (needs ANTHROPIC_API_KEY); "
        "off = no panel-ROI detection.",
    )
    model: Optional[str] = Field(
        default=None,
        description="Override the per-backend default vision model "
        "(e.g. claude-sonnet-4-6-20251001); used only by the vision-* modes.",
    )


class GrobidConfig(BaseModel):
    """Grobid metadata-extraction service (#59)."""

    model_config = ConfigDict(extra="forbid")

    url: str = Field(
        default="http://localhost:8070",
        description="Grobid service URL.",
    )
    disable: bool = Field(
        default=False,
        description="Skip Grobid even if reachable.",
    )


class BibliographyConfig(BaseModel):
    """Bibliography-build options (#59 / #64)."""

    model_config = ConfigDict(extra="forbid")

    enrich_bhl: bool = Field(
        default=False,
        description="Enrich pre-DOI references against the Biodiversity "
        "Heritage Library. Slow, rate-limited, network-dependent — opt-in.",
    )


class LicensingConfig(BaseModel):
    """Figure-licensing policy (#59 / #51)."""

    model_config = ConfigDict(extra="forbid")

    pd_cutoff_years: int = Field(
        default=95,
        ge=0,
        description="Public-domain cutoff in years before the current year. "
        "Default 95 is the US copyright rule; configurable per corpuscle.",
    )


# ---------------------------------------------------------------------------
# System-wide tuning blocks (carried over from today's _DEFAULT_CONFIG).
# ---------------------------------------------------------------------------


class OcrConfig(BaseModel):
    model_config = ConfigDict(extra="forbid")

    optimize_level: int = Field(default=2, ge=0, le=3)
    ocr_languages_default: List[str] = Field(
        default_factory=lambda: [
            "eng", "deu", "deu_latf", "fra", "rus", "lat",
            "spa", "por", "chi_sim", "chi_tra", "jpn", "ell", "kor",
        ]
    )
    gibberish_threshold: float = Field(default=0.65, ge=0.0, le=1.0)


class ChunkingConfig(BaseModel):
    model_config = ConfigDict(extra="forbid")

    max_tokens: int = Field(default=8191, gt=0)


class StageTimeoutsConfig(BaseModel):
    model_config = ConfigDict(extra="forbid")

    ocr: int = Field(default=1800, gt=0)
    grobid: int = Field(default=300, gt=0)
    docling: int = Field(default=600, gt=0)
    vision_pass: int = Field(default=600, gt=0)


class HugeDocumentConfig(BaseModel):
    model_config = ConfigDict(extra="forbid")

    max_pages: int = Field(default=5000, gt=0)


class QualityGatesConfig(BaseModel):
    model_config = ConfigDict(extra="forbid")

    empty_text_min_chars: int = Field(default=500, ge=0)
    min_chars_per_page: int = Field(default=200, ge=0)
    max_gibberish_score: float = Field(default=0.5, ge=0.0, le=1.0)
    zero_refs_min_pages: int = Field(default=5, ge=0)
    min_median_chunk_chars: int = Field(default=50, ge=0)
    min_figure_mean_intensity: float = Field(default=10, ge=0.0)
    max_figures_sampled: int = Field(default=50, gt=0)


# ---------------------------------------------------------------------------
# Top-level corpuscle config.
# ---------------------------------------------------------------------------


class CorpuscleConfig(BaseModel):
    """The full config.yaml for one corpuscle directory.

    All per-corpuscle inputs (input_pdfs, bib, lexicon, taxonomy, ...)
    fold in here, replacing the long CLI flag lists today's pipeline
    carries. System-wide tuning blocks (ocr, chunking, ...) keep their
    v0.2 shape.
    """

    model_config = ConfigDict(extra="forbid")

    # Per-corpuscle inputs (new in v0.3).
    input_pdfs: Optional[Path] = Field(
        default=None,
        description="Directory of PDFs to process (recursively walked). "
        "Required for `corpus run`; resolved relative to this config "
        "file's directory (#61 path resolution). Glob expansion is "
        "not currently supported — point at a directory.",
    )
    output_dir: Path = Field(
        default=Path("./output"),
        description="Per-corpuscle artifact tree (per-paper directories, "
        "cross-paper SQLite DBs, LanceDB index). Resolved relative to "
        "the config file.",
    )
    bib: Optional[Path] = Field(
        default=None,
        description="Curated BibTeX file (overrides Grobid header parse "
        "for matching PDFs).",
    )
    lexicon: Optional[Path] = Field(
        default=None,
        description="Multi-category lexicon YAML (anatomy, biogeography, "
        "...). Each category emits its own <hash>/<category>.json artifact.",
    )

    taxonomy: TaxonomyConfig = Field(default_factory=TaxonomyConfig)
    figures: FiguresConfig = Field(default_factory=FiguresConfig)
    grobid: GrobidConfig = Field(default_factory=GrobidConfig)
    bibliography: BibliographyConfig = Field(default_factory=BibliographyConfig)
    licensing: LicensingConfig = Field(default_factory=LicensingConfig)

    # System-wide tuning blocks (carried from v0.2 _DEFAULT_CONFIG).
    ocr: OcrConfig = Field(default_factory=OcrConfig)
    chunking: ChunkingConfig = Field(default_factory=ChunkingConfig)
    stage_timeouts: StageTimeoutsConfig = Field(default_factory=StageTimeoutsConfig)
    huge_document: HugeDocumentConfig = Field(default_factory=HugeDocumentConfig)
    quality_gates: QualityGatesConfig = Field(default_factory=QualityGatesConfig)

    @model_validator(mode="before")
    @classmethod
    def _reject_legacy_vision_block(cls, data):
        """The v0.5 ``vision:`` block was renamed to ``figures:`` in v0.6
        (#102). Fail loud with the mapping rather than letting
        ``extra=forbid`` emit a bare "unexpected field" error."""
        if isinstance(data, dict) and "vision" in data:
            raise ValueError(
                "the `vision:` config block was replaced by `figures:` in "
                "v0.6 (#102). Rename it and map `vision.backend` → "
                "`figures.panel_detection`: none → `off` (or `ocr` for the "
                "new CPU OCR-panel floor), local → `vision-local`, claude → "
                "`vision-claude`. Move `vision.model` → `figures.model`."
            )
        return data

    @field_validator("taxonomy")
    @classmethod
    def _taxonomy_consistent(cls, v: TaxonomyConfig) -> TaxonomyConfig:
        if v.source == "worms" and v.root_id is None:
            raise ValueError(
                "taxonomy.root_id is required when taxonomy.source = 'worms'"
            )
        if v.source in {"dwc", "dwca"} and v.path is None:
            raise ValueError(
                f"taxonomy.path is required when taxonomy.source = {v.source!r}"
            )
        return v


def validate_config(merged: dict) -> CorpuscleConfig:
    """Validate a merged config dict against :class:`CorpuscleConfig`.

    Field-level errors point at the exact key + value (e.g.,
    ``vision.backend  Input should be 'local', 'claude' or 'none'``).
    Caller decides how to surface the ValidationError — `corpus check`
    pretty-prints it; library-internal callers can let it propagate.
    """
    return CorpuscleConfig.model_validate(merged)


__all__ = [
    "BibliographyConfig",
    "ChunkingConfig",
    "CorpuscleConfig",
    "GrobidConfig",
    "HugeDocumentConfig",
    "LicensingConfig",
    "OcrConfig",
    "QualityGatesConfig",
    "StageTimeoutsConfig",
    "TaxonomyConfig",
    "ValidationError",
    "VisionConfig",
    "validate_config",
]
