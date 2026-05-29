"""Output-type profiles (#101).

Output type is a **client/session property carried per call, never a
corpus or server attribute** — a shared SSE server fields many clients
at once (chat, internal report, manuscript), so a server-start flag or a
per-corpuscle config default cannot distinguish them. The vocabulary is
a fixed set of global built-ins; selection is the per-call ``profile=``
argument on the gated tools, with a server-level *fallback*
(``--default-profile``) only for calls that omit it.

Freeze-critical core (#101): the ``profile=`` parameter, this stable
vocabulary, and the **figure-licensing** gate. The other policy axes are
carried on the profile for forward-compatibility but their enforcement
is additive and may land post-1.0:

* ``require_attribution`` — surfaced in the figure-URL payload today;
  full inline-emission is a follow-up.
* ``citation_provenance`` — ``strict`` enforcement depends on #100.
* ``excerpt_max_words`` — not yet enforced.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Optional


@dataclass(frozen=True)
class OutputProfile:
    name: str
    figure_licensing: str          # "permissive" | "strict"
    require_attribution: bool
    citation_provenance: str       # "warn" | "strict"
    excerpt_max_words: Optional[int]

    def as_dict(self) -> Dict:
        return {
            "name": self.name,
            "figure_licensing": self.figure_licensing,
            "require_attribution": self.require_attribution,
            "citation_provenance": self.citation_provenance,
            "excerpt_max_words": self.excerpt_max_words,
        }


BUILTIN_PROFILES: Dict[str, OutputProfile] = {
    "report": OutputProfile(
        "report", "permissive", False, "warn", None),
    "manuscript": OutputProfile(
        "manuscript", "strict", True, "strict", 90),
    "presentation": OutputProfile(
        "presentation", "strict", True, "warn", None),
}

# Server fallback when a call omits profile=. Permissive by design (#94):
# in-chat figure display is fair use, and refusing it is the common-case
# pain. Publication-bound clients pass profile="manuscript" per call.
DEFAULT_PROFILE = "report"


def get_profile(name: Optional[str]) -> Optional[OutputProfile]:
    """Return the named built-in profile, or None if unknown/None."""
    if name is None:
        return None
    return BUILTIN_PROFILES.get(name)


def resolve_profile(
    per_call: Optional[str], server_default: Optional[str],
) -> OutputProfile:
    """Active profile: per-call wins, then the server default, then the
    built-in ``report``. Unknown names fall through — callers validate
    ``per_call`` up front (via :func:`get_profile`) for a clear error."""
    for candidate in (per_call, server_default, DEFAULT_PROFILE):
        prof = get_profile(candidate)
        if prof is not None:
            return prof
    return BUILTIN_PROFILES[DEFAULT_PROFILE]


def unknown_profile_error(name: str) -> Dict:
    # Uniform error shape (Phase 3 freeze gate): error message + machine
    # code + tool-specific context.
    return {
        "error": f"unknown profile {name!r}",
        "code": "invalid_argument",
        "available": sorted(BUILTIN_PROFILES),
    }
