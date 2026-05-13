# round2_paper

A single siphonophore PDF (`Siebert_etal2011.pdf`) held outside the
demo's `input_pdfs` scope so the operator (or CI) can exercise the
implicit-resume "add a paper and re-run" flow on top of the
4-paper [`demo/`](../../../demo/) corpuscle.

Usage:

```bash
# Round 1: build the demo from scratch on the 4 bundled papers.
cd demo && corpus run --no-vision

# Round 2: drop the held-back paper into demo/ and re-run. Implicit
# resume only re-processes the new paper; the other 4 cache-hit.
cp ../tests/fixtures/round2_paper/Siebert_etal2011.pdf .
corpus run --no-vision

# Cleanup (so the working tree matches the committed shape):
rm Siebert_etal2011.pdf
rm -rf output     # demo/output is gitignored regardless
```

Why this lives here instead of `demo/_round2/`: the demo's
`config.yaml` sets `input_pdfs: .` which rglobs the whole `demo/`
tree. A subdirectory inside `demo/` would be picked up on round 1,
defeating the held-back purpose.

The paper has a corresponding entry in
[`demo/siphonophores.bib`](../../../demo/siphonophores.bib) so its
metadata enrichment via `--bib` works on round 2 without any extra
configuration.

See [issue #74](https://github.com/caseywdunn/corpus/issues/74) for
the rationale behind the 4 + 1 fixture shape.
