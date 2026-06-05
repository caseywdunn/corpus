"""Output-profile discovery tools (#101).

The profile vocabulary lives in :mod:`mcpsrv.profiles`; these tools let
a client enumerate it and read the server's fallback so it can decide
which ``profile=`` to pass on the gated figure (and, later, citation)
tools. Selection itself is per-call — these tools are read-only.
"""
from __future__ import annotations

from typing import Dict

from ..app import _need_index, mcp
from ..profiles import BUILTIN_PROFILES, DEFAULT_PROFILE, resolve_profile


@mcp.tool()
def list_output_profiles() -> Dict:
    """List the built-in output profiles (#101) and their policies.

    An output profile is a client/session property you pass per call as
    ``profile=`` to the gated tools (figure image/URL today). It governs
    the figure-licensing gate: ``report`` (permissive — in-chat display)
    vs ``manuscript`` / ``presentation`` (strict — only publishable
    figures). Selection is per call, so concurrent sessions against the
    same server can each use a different profile.

    Returns ``{profiles: [{name, figure_licensing, require_attribution,
    citation_provenance, excerpt_max_words}, ...], server_default}``.
    """
    idx = _need_index()
    server_default = getattr(idx, "default_profile", None) or DEFAULT_PROFILE
    return {
        "profiles": [p.as_dict() for p in BUILTIN_PROFILES.values()],
        "server_default": server_default,
    }


@mcp.tool()
def get_active_profile() -> Dict:
    """Return the server's fallback output profile — the one applied to
    calls that omit ``profile=`` (#101).

    This is informational: the authoritative selection is the per-call
    ``profile=`` argument, so this reports the *default*, not a session
    state. Returns the resolved profile's policy dict plus
    ``is_server_default: true``.
    """
    idx = _need_index()
    active = resolve_profile(None, getattr(idx, "default_profile", None))
    return {**active.as_dict(), "is_server_default": True}
