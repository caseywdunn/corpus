"""Grobid capability versus persisted evidence, owned by the build plane."""

# Internal request policy, not an operator tuning knob. Increment when the
# persisted evidence contract changes so same-version development builds and
# their subordinate TEI caches cannot silently retain the previous contract.
REFERENCE_EVIDENCE_VERSION = 1


def grobid_input(context, hash_dir=None):
    """Resolve the metadata expectation without downgrading proof on outages.

    A live service expects complete extraction. Fallback receipts deliberately
    disagree with that expectation so recovery retries them. While offline (or
    during a network-free dry run), retain the last enabled receipt's evidence:
    an outage is not a reason to discard a valid result or churn placeholders.
    Configuration changes are compared separately and still invalidate caches.
    """
    if not context["enabled"]:
        return {"status": "disabled", "service_version": None}
    if context.get("available"):
        return {"status": "complete", "service_version": context.get("service_version")}
    if hash_dir is not None:
        from .stages import _load_pipeline_state
        record = _load_pipeline_state(hash_dir)["stages"].get("metadata_extraction") or {}
        prior = (record.get("input_fingerprint") or {}).get("grobid") or {}
        if prior.get("status") in ("complete", "unavailable"):
            return dict(prior)
    return {"status": "unavailable", "service_version": context.get("service_version")}
