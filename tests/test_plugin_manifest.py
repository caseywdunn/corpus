"""The skills plugin is well-formed, and its two version strings agree (#178).

Skills ship as a Claude Code plugin rather than as project skills, because
almost nothing they do happens inside this repo — a library is assembled in its
own repo, a monograph written in another, and `.claude/skills/` is discovered
only when the cwd is here.

That buys reach and costs a failure mode nothing else in this repo has: the
plugin is consumed by `claude plugin install`, not by pytest or by an import,
so a malformed manifest or a skill with the wrong `name:` fails at a user's
machine rather than in CI. The version appears in *two* files that must agree;
bumping one and forgetting the other is the obvious mistake and is invisible
locally.

So this asserts the contract CONTRIBUTING.md's "Adding a skill" describes:
manifests parse, versions agree, every skill directory has a `SKILL.md`, and
each declares a `name` matching its directory plus a non-empty `description`.

The description check is not cosmetic. It is the only text the model selects a
skill on, and these skills have deliberately adjacent subjects — summarizing a
corpuscle, building one, writing a monograph from one.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest
import yaml

REPO = Path(__file__).resolve().parent.parent
PLUGIN_DIR = REPO / ".claude-plugin"
PLUGIN_JSON = PLUGIN_DIR / "plugin.json"
MARKETPLACE_JSON = PLUGIN_DIR / "marketplace.json"
SKILLS_DIR = REPO / "skills"

pytestmark = pytest.mark.skipif(
    not PLUGIN_DIR.exists(), reason="plugin scaffolding not present (sdist/wheel)"
)


def _load(path: Path) -> dict:
    with path.open(encoding="utf-8") as fh:
        return json.load(fh)


def _skill_dirs() -> list[Path]:
    if not SKILLS_DIR.exists():
        return []
    return sorted(d for d in SKILLS_DIR.iterdir() if d.is_dir())


def _frontmatter(skill_md: Path) -> dict:
    """Parse the leading `---`-delimited YAML block of a SKILL.md."""
    text = skill_md.read_text(encoding="utf-8")
    if not text.startswith("---"):
        raise AssertionError(f"{skill_md.relative_to(REPO)} has no frontmatter block")
    _, _, rest = text.partition("---")
    block, sep, _ = rest.partition("\n---")
    if not sep:
        raise AssertionError(
            f"{skill_md.relative_to(REPO)} frontmatter block is not closed with ---"
        )
    parsed = yaml.safe_load(block)
    if not isinstance(parsed, dict):
        raise AssertionError(
            f"{skill_md.relative_to(REPO)} frontmatter is not a mapping"
        )
    return parsed


def test_manifests_exist_and_parse():
    assert PLUGIN_JSON.exists(), "plugin.json is what `claude plugin install` reads"
    assert MARKETPLACE_JSON.exists(), "marketplace.json is what `marketplace add` reads"
    _load(PLUGIN_JSON)
    _load(MARKETPLACE_JSON)


def test_marketplace_lists_this_plugin_at_the_same_version():
    """The version lives in two files. They drift silently; this is the check."""
    plugin = _load(PLUGIN_JSON)
    market = _load(MARKETPLACE_JSON)

    entries = [p for p in market.get("plugins", []) if p.get("name") == plugin["name"]]
    assert entries, (
        f"marketplace.json lists no plugin named {plugin['name']!r}; "
        f"found {[p.get('name') for p in market.get('plugins', [])]}"
    )
    assert len(entries) == 1, f"{plugin['name']!r} listed more than once"

    assert entries[0].get("version") == plugin.get("version"), (
        f"version drift: plugin.json says {plugin.get('version')!r}, "
        f"marketplace.json says {entries[0].get('version')!r}. "
        "CONTRIBUTING.md's 'Adding a skill' checklist says to bump both."
    )


def test_plugin_declares_a_name_and_version():
    plugin = _load(PLUGIN_JSON)
    for field in ("name", "version", "description"):
        assert plugin.get(field), f"plugin.json is missing {field!r}"


@pytest.mark.parametrize(
    "skill_dir", _skill_dirs(), ids=lambda d: d.name if hasattr(d, "name") else str(d)
)
def test_every_skill_directory_has_a_well_formed_skill_md(skill_dir: Path):
    skill_md = skill_dir / "SKILL.md"
    assert skill_md.exists(), (
        f"{skill_dir.relative_to(REPO)} has no SKILL.md — a skill directory without "
        "one is invisible to Claude Code and silently does nothing"
    )

    fm = _frontmatter(skill_md)

    assert fm.get("name") == skill_dir.name, (
        f"{skill_md.relative_to(REPO)} declares name={fm.get('name')!r} but lives in "
        f"{skill_dir.name!r}; they must match"
    )

    description = (fm.get("description") or "").strip()
    assert description, (
        f"{skill_md.relative_to(REPO)} has no description — that is the only text "
        "the model selects this skill on"
    )
