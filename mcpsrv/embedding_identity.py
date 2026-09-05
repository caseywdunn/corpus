"""Resolve build/query compatibility at startup; never mutate a served bundle."""
import json


def load_identity(root, manifest):
    path = root / "embedding_producer.json"
    if path.exists():
        identity = json.loads(path.read_text())
        if (not isinstance(identity, dict) or identity.get("schema_version") != 1
                or not isinstance(identity.get("model"), str) or not identity["model"]
                or type(identity.get("dimension")) is not int or identity["dimension"] <= 0
                or not isinstance(identity.get("producer"), dict)
                or not identity["producer"].get("verification")):
            raise ValueError("Malformed embedding producer sidecar")
        producer = identity["producer"]
        if producer.get("verification") == "repository-revision" and (
                not isinstance(producer.get("revision"), str) or not producer["revision"]):
            raise ValueError("Malformed embedding repository revision")
        if manifest and (identity["model"] != manifest.get("embedding_model")
                         or identity["dimension"] != manifest.get("embedding_dim")):
            raise ValueError("Embedding producer sidecar disagrees with bundle manifest")
        return identity
    if manifest:
        return {"model": manifest.get("embedding_model"), "dimension": manifest.get("embedding_dim"),
                "producer": None}
    # Build outputs may be served for development. Check all existing document
    # receipts during startup, not a sampled marker in an active query.
    from pipeline.embedding_state import embedding_identity
    return embedding_identity(root, require_all=True)
