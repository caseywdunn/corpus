"""Citation-grounding quality tests for #79.

Runs each ``prompts:`` entry in ``tests/ground_truth/*.yaml`` against
the served corpus via a real Claude API round-trip, with the
``format_citations`` MCP tool wired in as a callable. Verifies that
every expected citation is emitted as one of ``format_citations``'s
per-item ``formatted`` outputs (i.e., the model used the tool and
pasted the output verbatim, per ``default_instructions.md``), and that
no ``forbidden_hallucinations`` patterns appear in the response.

This is the *behavioural* check the format_citations tool exists to
support — does an LLM client following the instructions actually
avoid hallucinated citations against this corpus, end to end. The
in-process tests in ``test_format_citation_tool.py`` pin the tool's
own contract; this file pins the client-side discipline.

**Gating** — cost-sensitive (real API calls per prompt). Skipped
unless ``RUN_PROMPT_QUALITY=1``. Within CI, intended for the
release-time T2 lane with ``ANTHROPIC_API_KEY`` in the runner's
secrets.

**Ground-truth shape** (``prompts:`` block in each per-paper YAML):

    prompts:
      - id: marrus_synopsis
        prompt: "Summarize the discovery of Marrus claudanielis with
                 citations from the corpus."
        expected_citations:
          - work_id: "corpus:..."
          - work_id: "10.1234/..."
        forbidden_hallucinations:
          # Real taxonomist-feedback incident: the paper appeared in
          # a generated reference list attributed to the wrong author
          # and the wrong journal. This is the negative trip-wire.
          - author_surname: "Schneider"
            attributed_journal_not: "Bull. Mus. Comp. Zool."
        expected_failure: false  # set true for known gaps; xfails
"""
from __future__ import annotations

import json
import os
import re
import types
from pathlib import Path
from typing import Any, Dict, List, Tuple

import pytest
import yaml

REPO_ROOT = Path(__file__).resolve().parent.parent
GROUND_TRUTH_DIR = REPO_ROOT / "tests" / "ground_truth"

# Cap the round-trip turn count so a runaway tool-use loop can't
# burn the API budget. Even verbose prompts should resolve in well
# under this.
_MAX_TURNS = 20

# Default model for the round-trip. Haiku for cost — these tests
# don't exercise the model's deep reasoning; they exercise its
# discipline at routing citations through the tool.
_DEFAULT_MODEL = "claude-haiku-4-5-20251001"


# ---------------------------------------------------------------------------
# Discovery — runs at pytest collect time, must be cheap and robust
# ---------------------------------------------------------------------------


def _discover_prompts() -> List[Tuple[str, Dict[str, Any]]]:
    """Read every ground-truth YAML; yield ``(yaml_basename, prompt_dict)``.

    Robust to YAMLs without a ``prompts:`` block (today: all three
    of them — this file is the first to consume the block). Returns
    an empty list if the directory is missing entirely.
    """
    if not GROUND_TRUTH_DIR.is_dir():
        return []
    out: List[Tuple[str, Dict[str, Any]]] = []
    for path in sorted(GROUND_TRUTH_DIR.glob("*.yaml")):
        try:
            data = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
        except yaml.YAMLError:
            continue  # not our problem here — other tests will flag it
        for prompt in data.get("prompts") or []:
            if "id" in prompt and "prompt" in prompt:
                out.append((path.name, prompt))
    return out


_PROMPTS = _discover_prompts()


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


def _gate_skip_reason() -> str | None:
    """Return a skip reason string, or None if the test should run.

    Env var set + API key missing is an actionable misconfiguration,
    not a quiet skip — those should fail loud."""
    if os.environ.get("RUN_PROMPT_QUALITY") != "1":
        return "Set RUN_PROMPT_QUALITY=1 to opt into the citation-quality lane"
    return None


@pytest.fixture(scope="session")
def biblio_db_path() -> Path:
    """Locate biblio_authority.sqlite from CORPUS_OUTPUT_DIR (or the
    demo's output as the default). Skips if the file is missing —
    can't run these tests against a corpus that hasn't been built.
    """
    default = REPO_ROOT / "demo" / "output"
    output_dir = Path(os.environ.get("CORPUS_OUTPUT_DIR", default))
    biblio = output_dir / "biblio_authority.sqlite"
    if not biblio.is_file():
        pytest.skip(
            f"biblio_authority.sqlite not found at {biblio}; run "
            f"`corpus run` against the demo (or set CORPUS_OUTPUT_DIR) first."
        )
    return biblio


@pytest.fixture(scope="session")
def mcp_index_injected(biblio_db_path: Path):
    """Inject a BiblioAuthority into the MCP module-global so the
    in-process format_citations tool call path resolves against the
    real demo build."""
    from mcpsrv import app as mcp_app
    from mcpsrv.indexes import BiblioAuthority

    fake_index = types.SimpleNamespace(biblio_db=BiblioAuthority(biblio_db_path))
    original = mcp_app._INDEX
    mcp_app.set_index(fake_index)
    yield fake_index
    mcp_app.set_index(original)


@pytest.fixture(scope="session")
def anthropic_client():
    """Lazy client construction so the import + env-var check don't
    fire on collection of skipped tests."""
    if not os.environ.get("ANTHROPIC_API_KEY"):
        pytest.fail(
            "RUN_PROMPT_QUALITY=1 is set but ANTHROPIC_API_KEY is not — "
            "the citation-quality lane needs both. Either set the key or "
            "unset RUN_PROMPT_QUALITY to skip."
        )
    import anthropic
    return anthropic.Anthropic()


# ---------------------------------------------------------------------------
# Tool definition + round-trip
# ---------------------------------------------------------------------------


def _format_citations_tool_def() -> Dict[str, Any]:
    """Anthropic-API tool definition mirroring the batched
    format_citations MCP tool. The description nudges the model to
    actually use it — the system prompt (from default_instructions.md)
    carries the verbatim-paste rule.
    """
    return {
        "name": "format_citations",
        "description": (
            "Return fully-assembled citation strings for works in the "
            "corpus bibliographic authority DB. Use for EVERY citation "
            "you emit — never recombine author/year/journal/title in your "
            "own context. Pass exactly one of: queries (free-text 'Author "
            "Year [Title]' strings), work_ids (canonical DOI / corpus: / "
            "bhl:), paper_hashes (12-hex SHA-256 prefixes of corpus "
            "papers) — each a list. Returns citations[] in input order. "
            "Batch a whole reference list into one call."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "queries": {"type": "array", "items": {"type": "string"}},
                "work_ids": {"type": "array", "items": {"type": "string"}},
                "paper_hashes": {"type": "array", "items": {"type": "string"}},
            },
        },
    }


def _system_prompt() -> str:
    """Use the same default_instructions.md the served corpus ships."""
    p = REPO_ROOT / "mcpsrv" / "default_instructions.md"
    return p.read_text(encoding="utf-8")


def _run_round_trip(
    client, prompt_text: str, system: str, model: str,
) -> Tuple[str, List[Dict[str, Any]]]:
    """Loop tool-use turns until the model emits a final response.

    Returns (final_text, tool_calls) where tool_calls is the list of
    ``{"input": ..., "output": ...}`` records — one per *resolved
    citation*. A single batched ``format_citations`` call that returns
    N citations expands into N records here, so the per-item scorer
    downstream stays unchanged.
    """
    from mcpsrv.tools.bibliography import format_citations

    messages: List[Dict[str, Any]] = [{"role": "user", "content": prompt_text}]
    tool_calls: List[Dict[str, Any]] = []
    tools = [_format_citations_tool_def()]

    for _ in range(_MAX_TURNS):
        resp = client.messages.create(
            model=model,
            max_tokens=4000,
            system=system,
            tools=tools,
            messages=messages,
        )

        # Always append the assistant message; tool_result content
        # must reference its tool_use ids.
        messages.append({"role": "assistant", "content": resp.content})

        tool_results: List[Dict[str, Any]] = []
        for block in resp.content:
            if block.type == "tool_use" and block.name == "format_citations":
                try:
                    result = format_citations(**(block.input or {}))
                except Exception as e:  # surface as tool error, not crash
                    result = {"error": "tool_exception", "detail": str(e)}
                # Flatten the batch into one record per citation so the
                # scorer can match each work_id / formatted / inline.
                for cite in (result.get("citations") or []):
                    tool_calls.append(
                        {"input": dict(block.input or {}), "output": cite}
                    )
                tool_results.append({
                    "type": "tool_result",
                    "tool_use_id": block.id,
                    "content": json.dumps(result),
                })

        if resp.stop_reason != "tool_use":
            final = "".join(b.text for b in resp.content if b.type == "text")
            return final, tool_calls

        messages.append({"role": "user", "content": tool_results})

    raise RuntimeError(
        f"model exceeded {_MAX_TURNS} tool-use turns — runaway loop"
    )


# ---------------------------------------------------------------------------
# Scoring
# ---------------------------------------------------------------------------


# Parenthetical citation regex: matches "(Surname, 2010)",
# "(Surname et al., 2010)", "(A & B, 2010)". Used to count
# candidate-citation tokens in the final response so we can detect
# tokens that DIDN'T come from a logged tool call (likely
# hand-written ⇒ hallucinated).
_PAREN_CITE_RE = re.compile(
    r"\(([A-Z][A-Za-zÀ-ÿ\-]+(?:\s+et\s+al\.)?"
    r"(?:\s+&\s+[A-Z][A-Za-zÀ-ÿ\-]+)?,?\s+\d{4}[a-z]?)\)"
)


def _score(
    prompt_spec: Dict[str, Any],
    final_text: str,
    tool_calls: List[Dict[str, Any]],
) -> Dict[str, Any]:
    """Compare emitted citations against the expected panel.

    A work counts as "emitted" iff format_citations returned an entry with
    a work_id that resolved correctly AND the returned ``formatted``
    string appears in ``final_text`` (i.e., the model pasted it
    verbatim).

    Hallucination detection scans for parenthetical-citation patterns
    in ``final_text`` that don't match any logged tool call's
    ``inline`` output. Crude but catches the recombination class the
    tool exists to prevent.
    """
    expected = {
        c["work_id"] for c in prompt_spec.get("expected_citations") or []
        if c.get("work_id")
    }

    emitted_work_ids: set[str] = set()
    inline_outputs: set[str] = set()
    for tc in tool_calls:
        out = tc.get("output") or {}
        if not isinstance(out, dict):
            continue
        work_id = out.get("work_id")
        formatted = out.get("formatted") or ""
        inline = out.get("inline") or ""
        if inline:
            inline_outputs.add(inline)
        if work_id and formatted and formatted in final_text:
            emitted_work_ids.add(work_id)

    expected_emitted = expected & emitted_work_ids
    missing = expected - emitted_work_ids

    # Hallucination scan: parenthetical citations in final_text that
    # weren't surfaced by a tool call. False-positive risk on quoted
    # excerpts that contain parentheticals, but those are rare in
    # synthesis outputs.
    hallucinated: List[str] = []
    for m in _PAREN_CITE_RE.finditer(final_text):
        token = m.group(0)
        if token not in inline_outputs:
            hallucinated.append(token)

    # Forbidden hallucinations: a name appears anywhere alongside a
    # forbidden journal attribution. Conservative pattern: surname
    # and forbidden-journal substring co-occur within a ~200-char
    # window (matches the typical reference-list entry length).
    forbidden_hits: List[Dict[str, str]] = []
    for fh in prompt_spec.get("forbidden_hallucinations") or []:
        surname = (fh.get("author_surname") or "").strip()
        journal_not = (fh.get("attributed_journal_not") or "").strip()
        if not surname or not journal_not:
            continue
        window = 200
        for m in re.finditer(re.escape(surname), final_text):
            chunk = final_text[m.start(): m.start() + window]
            if journal_not in chunk:
                forbidden_hits.append({
                    "author_surname": surname,
                    "attributed_journal_not": journal_not,
                    "context": chunk,
                })

    return {
        "expected": sorted(expected),
        "emitted_work_ids": sorted(emitted_work_ids),
        "expected_emitted": sorted(expected_emitted),
        "missing": sorted(missing),
        "hallucinated": hallucinated,
        "forbidden_hits": forbidden_hits,
    }


# ---------------------------------------------------------------------------
# Parametrized test
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "yaml_name, prompt_spec",
    _PROMPTS,
    ids=[f"{name}::{p['id']}" for name, p in _PROMPTS] or None,
)
def test_prompt_grounds_citations(
    yaml_name: str,
    prompt_spec: Dict[str, Any],
    mcp_index_injected,
    anthropic_client,
):
    """Per-prompt verdict.

    Pass conditions:
    - No ``forbidden_hallucinations`` patterns appear in the response.
    - Every ``expected_citations`` work_id was emitted via the tool.
    - No parenthetical-citation tokens appear that don't trace to a
      tool call.

    Prompts with ``expected_failure: true`` xfail (used for known
    gaps so a regression elsewhere doesn't get masked).
    """
    skip = _gate_skip_reason()
    if skip:
        pytest.skip(skip)

    model = os.environ.get("PROMPT_QUALITY_MODEL", _DEFAULT_MODEL)
    final_text, tool_calls = _run_round_trip(
        anthropic_client, prompt_spec["prompt"], _system_prompt(), model,
    )
    score = _score(prompt_spec, final_text, tool_calls)

    # Forbidden hallucinations: always a hard fail, regardless of
    # expected_failure. The whole point of those trip-wires is to
    # block known regressions on the way in.
    assert not score["forbidden_hits"], (
        f"{yaml_name}::{prompt_spec['id']}: forbidden hallucination(s) "
        f"in response — {score['forbidden_hits']}"
    )

    if prompt_spec.get("expected_failure"):
        # xfail: at least one assertion below should fail. If the
        # prompt passes cleanly, the maintainer should flip the flag
        # off.
        if not score["missing"] and not score["hallucinated"]:
            pytest.fail(
                f"{yaml_name}::{prompt_spec['id']} unexpectedly passed; "
                f"remove expected_failure: true"
            )
        pytest.xfail(
            f"known gap: missing={score['missing']} "
            f"hallucinated={score['hallucinated']}"
        )

    assert not score["missing"], (
        f"{yaml_name}::{prompt_spec['id']}: expected citations not emitted "
        f"via format_citations — {score['missing']}. Final response excerpt: "
        f"{final_text[:300]!r}"
    )
    assert not score["hallucinated"], (
        f"{yaml_name}::{prompt_spec['id']}: parenthetical citations in "
        f"response don't trace to a format_citations call — "
        f"{score['hallucinated']}. The model hand-wrote a citation; "
        f"this is the hallucination class #79 exists to prevent."
    )


# ---------------------------------------------------------------------------
# Unit tests on the harness machinery itself — always run
# ---------------------------------------------------------------------------


def test_paren_cite_regex_catches_common_shapes():
    cases = [
        "(Totton, 1965)",
        "(Smith & Jones, 2010)",
        "(Smith et al., 2010)",
        "(Smith et al., 2010a)",
    ]
    for c in cases:
        assert _PAREN_CITE_RE.search(c), c


def test_score_marks_missing_when_expected_not_emitted():
    spec = {"expected_citations": [{"work_id": "10.1/foo"}]}
    score = _score(spec, final_text="", tool_calls=[])
    assert score["missing"] == ["10.1/foo"]


def test_score_marks_emitted_only_when_formatted_in_final_text():
    """The model has to PASTE the tool output, not just call the tool.
    A logged tool call whose formatted string never appears in the
    response shouldn't count as emitted."""
    spec = {"expected_citations": [{"work_id": "10.1/foo"}]}
    score = _score(spec, final_text="some other text", tool_calls=[
        {"input": {"work_id": "10.1/foo"},
         "output": {"work_id": "10.1/foo", "formatted": "Smith (2010)."}},
    ])
    assert score["missing"] == ["10.1/foo"]
    assert score["expected_emitted"] == []


def test_score_emitted_when_formatted_appears():
    spec = {"expected_citations": [{"work_id": "10.1/foo"}]}
    score = _score(spec, final_text="As shown by Smith (2010).", tool_calls=[
        {"input": {"work_id": "10.1/foo"},
         "output": {"work_id": "10.1/foo", "formatted": "Smith (2010)."}},
    ])
    assert score["missing"] == []
    assert score["expected_emitted"] == ["10.1/foo"]


def test_score_detects_unsourced_parenthetical():
    """Parenthetical citation in response that doesn't match any tool
    call's `inline` field is the hallucination class to catch."""
    score = _score(
        prompt_spec={},
        final_text="As shown elsewhere (Hallucinated, 2020).",
        tool_calls=[],
    )
    assert "(Hallucinated, 2020)" in score["hallucinated"]


def test_score_forbidden_hallucination_co_occurrence():
    spec = {"forbidden_hallucinations": [{
        "author_surname": "Schneider",
        "attributed_journal_not": "Bull. Mus. Comp. Zool.",
    }]}
    score = _score(
        prompt_spec=spec,
        final_text="Schneider (1891). Some title. Bull. Mus. Comp. Zool., 2, 1-50.",
        tool_calls=[],
    )
    assert len(score["forbidden_hits"]) == 1
    assert score["forbidden_hits"][0]["author_surname"] == "Schneider"


def test_score_forbidden_hallucination_no_hit_when_separated():
    """Surname + journal far apart in the text shouldn't trip the
    co-occurrence window."""
    spec = {"forbidden_hallucinations": [{
        "author_surname": "Schneider",
        "attributed_journal_not": "Bull. Mus. Comp. Zool.",
    }]}
    score = _score(
        prompt_spec=spec,
        final_text="Schneider published in Zeitschrift. " + ("X" * 300) +
                   " Some other paper in Bull. Mus. Comp. Zool.",
        tool_calls=[],
    )
    assert score["forbidden_hits"] == []
