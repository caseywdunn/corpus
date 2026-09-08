"""Say which dependency is actually missing, not which one we asked for.

The ML backends are reached through one import each — ``transformers``,
``sentence_transformers``, ``anthropic`` — and every one of those pulls a
dependency chain behind it. When a link in that chain is absent the
``ImportError`` surfaces at *our* import statement, so a handler that
names the package it asked for reports the wrong culprit.

Observed on a Mac with the conda env unactivated, where nothing in the
ML stack was installed::

    VisionBackendError: transformers >= 4.45 is required for the local
    VLM backend (pip install transformers>=4.45 qwen-vl-utils torch
    accelerate)

``torch`` was the missing module. The message sent the reader to install
a package that was not the problem, which is worse than a vague error:
it is a confident wrong instruction, and following it leaves the install
just as broken.

``ImportError.name`` carries the module that actually failed, and it
distinguishes the two cases cleanly:

===================================  =====================
situation                            ``exc.name``
===================================  =====================
``transformers`` absent              ``'transformers'``
``transformers`` present, no torch   ``'torch'``
symbol missing from ``transformers`` ``'transformers'``
===================================  =====================

So a missing symbol (an out-of-date package) and an absent package both
report the module we asked for, and a broken dependency chain reports
the link that broke. That is exactly the distinction the message needs.
"""
from __future__ import annotations

from typing import Optional

__all__ = ["missing_dependency_message"]


def missing_dependency_message(
    exc: ImportError,
    *,
    feature: str,
    module: str,
    install: str,
    requirement: Optional[str] = None,
) -> str:
    """Why ``exc`` happened, naming the module that actually failed.

    ``module`` is the import name we asked for (``sentence_transformers``,
    not the ``sentence-transformers`` on PyPI); ``requirement`` is how to
    write it in prose when that differs or carries a version bound;
    ``install`` is the pip argument list that fixes it.
    """
    requirement = requirement or module
    missing = getattr(exc, "name", None) or module
    if missing == module:
        # Either the package is absent or it is too old to carry the
        # symbol. Both are fixed by installing what `install` names, and
        # the version bound in `requirement` covers the second case.
        return f"{requirement} is required for {feature} (pip install {install})"
    return (
        f"{feature} requires {requirement}, and {module} was found — but "
        f"importing it failed because {missing!r} is not installed. That "
        f"is the package to install, not {module}: pip install {install}"
    )
