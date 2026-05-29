"""Single-source version string for the corpus pipeline.

Imported by every CLI and stamped into persistent artifacts
(`bundle_manifest.json`, MCP `bundle_info`) so any output tree carries
the code version that produced it.

Bumped by hand at release time. The release ritual in CONTRIBUTING.md
keeps this in sync with the git tag and CHANGELOG entry. PEP 440
pre-release suffixes (`.devN`, `aN`, `rcN`) on `dev` between releases.
"""

__version__ = "0.6.0.dev0"
