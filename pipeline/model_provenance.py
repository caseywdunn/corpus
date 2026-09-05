"""Offline vision producer identities; never download or instantiate a model."""
import hashlib
import importlib.metadata
import json
from pathlib import Path

DEFAULT_VISION_MODELS = {
    "vision-local": "Qwen/Qwen2.5-VL-7B-Instruct",
    "vision-claude": "claude-haiku-4-5-20251001",
}


def _cached_revision(model):
    try:
        from huggingface_hub import try_to_load_from_cache
        cached = try_to_load_from_cache(model, "config.json")
    except (ImportError, ValueError, OSError):
        return None
    if not isinstance(cached, str):
        return None
    snapshot = Path(cached).parent
    if snapshot.parent.name != "snapshots":
        return None
    return snapshot.name


def _local_files(directory):
    hashes = {}
    for path in sorted(directory.rglob("*")):
        relative = path.relative_to(directory)
        if any(part.startswith(".") or part == "__pycache__" for part in relative.parts):
            continue
        if not path.is_file():
            continue
        before = path.stat()
        with path.open("rb") as stream:
            hashes[relative.as_posix()] = hashlib.file_digest(stream, "sha256").hexdigest()
        after = path.stat()
        if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
            raise ValueError(f"Model file changed while fingerprinting: {path}")
    if not hashes:
        raise ValueError(f"No model files in {directory}")
    return {"file_count": len(hashes), "sha256": hashlib.sha256(
        json.dumps(hashes, sort_keys=True).encode()).hexdigest()}


def vision_producer(mode, model=None, *, resolved_revision=None):
    model = model or DEFAULT_VISION_MODELS[mode]
    result = {"model": model, "mode": mode,
              # Includes prompts, generation settings and response handling.
              "implementation_sha256": hashlib.sha256(
                  Path(__file__).with_name("vision.py").read_bytes()).hexdigest()}
    packages = ("anthropic",) if mode == "vision-claude" else ("transformers", "torch", "tokenizers", "accelerate", "qwen-vl-utils")
    result["packages"] = {}
    for name in packages:
        try:
            result["packages"][name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            result["packages"][name] = None
    if mode == "vision-claude":
        result["generation"] = {"max_tokens": 1024}
        return {**result, "verification": "remote-model-id-only"}
    result["generation"] = {"max_new_tokens": 1024, "max_pixels": 1003520, "min_pixels": 3136}
    directory = Path(model).expanduser()
    if directory.is_dir():
        return {**result, "verification": "local-file-content", **_local_files(directory)}
    revision = resolved_revision or _cached_revision(model)
    return {**result, "verification": "repository-revision" if revision else "unverified-cache-miss",
            "revision": revision}
