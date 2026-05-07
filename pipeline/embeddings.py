"""Embedding backend abstraction.

:class:`LocalBackend` runs a sentence-transformers model on a local
GPU (Apple Metal / CUDA) or CPU. Default model is **BGE-M3** (1024-dim,
multilingual, 8k context — covers the corpus's German / French /
Russian tail), pulled from HuggingFace and cached under
``~/.cache/huggingface/``. Downloads on first use only.

Failures surface as :class:`EmbeddingError` (no silent zero-vector
fallback). The caller typically catches at the per-document level and
skips a failed doc rather than poisoning the index.

Device selection for LocalBackend, in order of preference:

1. ``cuda`` — for Bouchet (or any NVIDIA GPU).
2. ``mps`` — Apple Silicon (Metal Performance Shaders).
3. ``cpu`` — graceful fallback; embedding the corpus this way is slow
   but functionally correct.

The same code runs on Mac (dev) and Bouchet (production) without
modification, just with different `device` resolved at startup.
"""

from __future__ import annotations

import logging
import os
from abc import ABC, abstractmethod
from typing import List, Optional

logger = logging.getLogger(__name__)


class EmbeddingError(RuntimeError):
    """Raised when the embedding backend fails (network, OOM, etc.)."""


# ---------------------------------------------------------------------------
# Base class
# ---------------------------------------------------------------------------


class EmbeddingBackend(ABC):
    """Abstract embedding backend.

    Subclasses expose:
      * ``dim`` — vector dimension (used to size the LanceDB schema)
      * ``model_name`` — human-readable identifier recorded in the index
      * ``embed(texts)`` — return one float vector per input string
    """

    @property
    @abstractmethod
    def dim(self) -> int: ...

    @property
    @abstractmethod
    def model_name(self) -> str: ...

    @abstractmethod
    def embed(self, texts: List[str]) -> List[List[float]]: ...


# ---------------------------------------------------------------------------
# Local sentence-transformers backend
# ---------------------------------------------------------------------------


_DEFAULT_LOCAL_MODEL = "BAAI/bge-m3"
_DEFAULT_LOCAL_DIM = 1024
# Per-model dim overrides (HuggingFace model card values). Kept as a small
# table so non-default models work without runtime model-introspection.
_LOCAL_MODEL_DIMS = {
    "BAAI/bge-m3": 1024,
    "BAAI/bge-large-en-v1.5": 1024,
    "BAAI/bge-base-en-v1.5": 768,
    "intfloat/multilingual-e5-large": 1024,
    "intfloat/multilingual-e5-base": 768,
    "nomic-ai/nomic-embed-text-v1.5": 768,
    "sentence-transformers/all-MiniLM-L6-v2": 384,
}


def detect_device() -> str:
    """Choose the best torch device available on this machine.

    Order: cuda → mps → cpu. Honours ``CORPUS_DEVICE`` env var as an
    override (set to ``"cpu"`` if you need to force CPU on a machine
    where MPS misbehaves on a particular model).
    """
    override = os.environ.get("CORPUS_DEVICE", "").strip().lower()
    if override in ("cuda", "mps", "cpu"):
        return override

    try:
        import torch
    except ImportError:
        return "cpu"

    if torch.cuda.is_available():
        return "cuda"
    if hasattr(torch.backends, "mps") and torch.backends.mps.is_available():
        return "mps"
    return "cpu"


class LocalBackend(EmbeddingBackend):
    """Sentence-transformers embeddings on a local accelerator.

    Loads the model once per process (it's a 0.5–2 GB allocation
    depending on the model). All ``embed()`` calls go through
    ``SentenceTransformer.encode`` which handles batching internally.
    Outputs are L2-normalized so the LanceDB cosine index uses inner
    product correctly.

    Parameters
    ----------
    model_name:
        HuggingFace model id, default ``BAAI/bge-m3``.
    device:
        Override the autodetected device. Pass ``"mps"``, ``"cuda"``,
        or ``"cpu"``. None → autodetect via :func:`detect_device`.
    batch_size:
        Encoding batch size. 32 is conservative and works on 16 GB
        unified-memory Macs; raise for higher-memory or CUDA GPUs.
    fp16:
        Use float16 weights. Saves ~50% memory; accuracy drop on BGE-M3
        is negligible for retrieval. Default True on CUDA/MPS, False on
        CPU (CPU fp16 is slower than fp32 on most platforms).
    """

    def __init__(
        self,
        model_name: str = _DEFAULT_LOCAL_MODEL,
        *,
        device: Optional[str] = None,
        batch_size: int = 32,
        fp16: Optional[bool] = None,
    ):
        try:
            from sentence_transformers import SentenceTransformer
        except ImportError as e:
            raise EmbeddingError(
                "sentence-transformers is required for the local backend "
                "(pip install sentence-transformers)"
            ) from e

        self._model_name = model_name
        self.device = device or detect_device()
        self.batch_size = batch_size
        if fp16 is None:
            fp16 = self.device in ("cuda", "mps")
        self._dim = _LOCAL_MODEL_DIMS.get(model_name, _DEFAULT_LOCAL_DIM)

        logger.info(
            "Loading local embedding model %s on device=%s fp16=%s",
            model_name, self.device, fp16,
        )
        try:
            kwargs = {}
            if fp16:
                # SentenceTransformer accepts torch_dtype since 3.x; older
                # versions need model_kwargs. Pass through model_kwargs for
                # broadest compatibility.
                kwargs["model_kwargs"] = {"torch_dtype": "float16"}
            self.model = SentenceTransformer(model_name, device=self.device, **kwargs)
        except Exception as e:
            raise EmbeddingError(
                f"Could not load sentence-transformers model {model_name!r} on "
                f"device={self.device}: {e}"
            ) from e

    @property
    def dim(self) -> int:
        return self._dim

    @property
    def model_name(self) -> str:
        return self._model_name

    def embed(self, texts: List[str]) -> List[List[float]]:
        if not texts:
            return []
        try:
            arr = self.model.encode(
                texts,
                batch_size=self.batch_size,
                normalize_embeddings=True,
                show_progress_bar=False,
                convert_to_numpy=True,
            )
        except Exception as e:
            raise EmbeddingError(
                f"Local encoder failed on a batch of {len(texts)} texts: {e}"
            ) from e
        # SentenceTransformer returns numpy float32 (or float16); LanceDB
        # accepts either via its Vector dtype. Convert to plain Python lists
        # for predictable downstream serialization.
        return [v.astype("float32").tolist() for v in arr]


# ---------------------------------------------------------------------------
# Factory
# ---------------------------------------------------------------------------


def get_embedder(model: Optional[str] = None, **kwargs) -> EmbeddingBackend:
    """Construct a :class:`LocalBackend`. ``model`` overrides the default
    HuggingFace model id; any extra ``kwargs`` are forwarded to the
    backend's constructor.
    """
    return LocalBackend(model or _DEFAULT_LOCAL_MODEL, **kwargs)
