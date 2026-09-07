"""Inflected surface forms for non-English lexicon translations (#165).

An enumerated surface-form set is the wrong shape for an inflecting
language. The lexicon lists `Luftblase`; Eschscholtz 1825 prints
`Luftblasen`. Vanhöffen 1906 prints `Schwimmglocken` 41 times against 7
of `Schwimmglocke`, so base-form-only matching found a minority of that
document's own mentions — and a German paper could report anatomy
coverage of exactly zero, which reads as "nothing here" rather than
"not indexed".

The endings are curated per language rather than derived by a stemmer,
and every generated form is added to the variant map explicitly, so
matching stays whole-word exact and the map stays greppable.
"""
from __future__ import annotations

import pytest

from pipeline.taxa import (
    _build_lexicon_matcher,
    _inflected_forms,
    extract_lexicon_mentions,
)


def _lex(**translations):
    return {
        "pneumatophore": {
            "synonyms": ["pneumatophores", "float"],
            "translations": translations,
            "description": "Apical gas-filled float.",
        }
    }


def _chunks(*texts):
    return [{"chunk_id": f"c{i}", "text": t} for i, t in enumerate(texts)]


# --- the reported case ----------------------------------------------------


def test_a_german_plural_is_matched():
    """The exact case in #165: Eschscholtz prints `Luftblasen`, the
    lexicon has `Luftblase`."""
    out = extract_lexicon_mentions(
        _chunks("Die Luftblasen sind gross."), _lex(de=["Luftblase"]),
    )
    assert out["total_mentions"] == 1
    assert out["mentions"][0]["canonical"] == "pneumatophore"
    assert out["mentions"][0]["matched_text"] == "Luftblasen"


@pytest.mark.parametrize("printed", [
    "Schwimmglocken",      # Vanhöffen 1906, 41 occurrences
    "Deckstücke", "Deckstücken", "Deckstückes", "Deckstücks",
    "Tentakeln", "Nesselkapseln", "Larven", "Tastern",
])
def test_the_measured_german_endings_all_match(printed):
    base = {"Schwimmglocken": "Schwimmglocke", "Deckstücke": "Deckstück",
            "Deckstücken": "Deckstück", "Deckstückes": "Deckstück",
            "Deckstücks": "Deckstück", "Tentakeln": "Tentakel",
            "Nesselkapseln": "Nesselkapsel", "Larven": "Larve",
            "Tastern": "Taster"}[printed]
    out = extract_lexicon_mentions(_chunks(f"und {printed} hier"),
                                   _lex(de=[base]))
    assert out["total_mentions"] == 1, f"{printed} from {base}"


@pytest.mark.parametrize("printed,base", [
    ("tentacules", "tentacule"), ("nématocystes", "nématocyste"),
    ("bractées", "bractée"), ("palpones", "palpone"),
])
def test_french_plurals_match(printed, base):
    out = extract_lexicon_mentions(_chunks(f"les {printed} sont"),
                                   _lex(fr=[base]))
    assert out["total_mentions"] == 1


@pytest.mark.parametrize("printed", [
    # Consonant-final stem: the ending stacks on.
    "нектофора", "нектофоры", "нектофоров", "нектофорами",
    "нектофором", "нектофорам", "нектофорах", "нектофору", "нектофоре",
])
def test_russian_case_endings_match_a_consonant_stem(printed):
    out = extract_lexicon_mentions(_chunks(f"два {printed} тут"),
                                   _lex(ru=["нектофор"]))
    assert out["total_mentions"] == 1, printed


@pytest.mark.parametrize("printed", ["личинки", "личинками", "личинкам"])
def test_russian_declension_replaces_a_stem_final_vowel(printed):
    """`личинка` declines to `личинки`, not `личинкаи` — so the ending
    has to reach the stem without its final vowel too."""
    out = extract_lexicon_mentions(_chunks(f"эти {printed} здесь"),
                                   _lex(ru=["личинка"]))
    assert out["total_mentions"] == 1, printed


@pytest.mark.parametrize("printed,base,lang", [
    # Russian genitive plural inserts a fill vowel into the stem:
    # личинка -> личин|о|к, which no ending appended to any prefix of
    # the stem produces.
    ("личинок", "личинка", "ru"),
    # German umlaut plurals change a stem vowel: Magen -> Mägen,
    # Fangfaden -> Fangfäden. Eschscholtz prints both.
    ("Saugmägen", "Saugmagen", "de"),
    ("Fangfäden", "Fangfaden", "de"),
])
def test_stem_changing_inflections_are_a_known_gap(printed, base, lang):
    """Recorded rather than silently absent. These need a per-language
    morphology table or the curator listing the form, and suffixing
    cannot reach them however long the ending list gets — so a lexicon
    that needs them should list them under `synonyms`."""
    out = extract_lexicon_mentions(_chunks(f"... {printed} ..."),
                                   _lex(**{lang: [base]}))
    assert out["total_mentions"] == 0
    # And the documented workaround does work.
    entry = _lex(**{lang: [base]})
    entry["pneumatophore"]["synonyms"].append(printed)
    assert extract_lexicon_mentions(
        _chunks(f"... {printed} ..."), entry)["total_mentions"] == 1


# --- what must NOT happen -------------------------------------------------


def test_english_is_not_inflected():
    """`cnida` + an ending matches `Cnidaria` 4,444 times and
    `cnidarian(s)` 3,740 more across the reference library — the phylum,
    not the nematocyst. English variants are hand-listed on purpose."""
    lex = {"cnida": {"synonyms": [], "translations": {}, "description": ""}}
    out = extract_lexicon_mentions(
        _chunks("The phylum Cnidaria contains cnidarians and cnidarian polyps"),
        lex,
    )
    assert out["total_mentions"] == 0
    assert _inflected_forms("cnida", "en") == []
    # Nor via a translation list that happens to be tagged English.
    assert _inflected_forms("float", "en") == []


def test_a_curated_form_is_never_displaced_by_a_generated_one():
    """Two entries where one's generated form collides with the other's
    hand-written one: the curated term wins."""
    lex = {
        "bract": {"synonyms": [], "translations": {"de": ["Brakteen"]},
                  "description": ""},
        "brakte_thing": {"synonyms": [], "translations": {"de": ["Brakte"]},
                         "description": ""},
    }
    _, variants = _build_lexicon_matcher(lex)
    # "Brakteen" is curated under `bract`; `Brakte` + "en" would generate
    # the same string for `brakte_thing`.
    assert variants["brakteen"] == "bract"


def test_short_stems_are_left_alone():
    """Below five characters a stem plus an ending is as likely to be an
    unrelated word. The shortest real translation form is `Larve`."""
    assert _inflected_forms("os", "de") == []
    assert _inflected_forms("aile", "fr") == []
    assert _inflected_forms("Larve", "de") != []


def test_multiword_and_hyphenated_translations_are_left_alone():
    """A phrase inflects on its head, not by taking an ending on the end
    of the phrase."""
    assert _inflected_forms("Schwimm- und Athmungshöhle", "de") == []
    assert _inflected_forms("cloche natatoire", "fr") == []


def test_an_unknown_language_tag_generates_nothing():
    assert _inflected_forms("Schwimmglocke", "zz") == []
    assert _inflected_forms("Schwimmglocke", "") == []
    assert _inflected_forms("Schwimmglocke", None) == []


def test_matching_stays_whole_word():
    """Generated forms go through the same whole-word alternation, so a
    generated form inside a longer word still does not match."""
    out = extract_lexicon_mentions(
        _chunks("Schwimmglockenzone und Luftblasenwand"),
        _lex(de=["Schwimmglocke", "Luftblase"]),
    )
    assert out["total_mentions"] == 0


def test_the_base_form_still_matches():
    out = extract_lexicon_mentions(
        _chunks("Die Luftblase und die Schwimmglocke"),
        _lex(de=["Luftblase", "Schwimmglocke"]),
    )
    assert out["total_mentions"] == 2
