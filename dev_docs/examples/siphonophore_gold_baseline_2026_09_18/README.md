# Siphonophore gold: intermediate CPU extraction baseline

This is a source-fidelity measurement of the existing 35-document gold set,
not final v1.5 acceptance. Extraction ran at `6fbf4e0` with OCR panel detection;
scorers ran at `a3e638d`. This predates surname v3, spacing v2 and regional PDF
CMap recovery. It includes neither Qwen inference nor embeddings/post-build
materialization. Do not compare it directly with a differently configured
historical release or treat token coverage as scientific-symbol accuracy.

All 35 document summaries report success, no stage failures, and completed
extraction, metadata, chunking, figure materialization/cross-reference and
annotation records. The application logged completion after 7056.6 seconds.
The retained tool session reported exit 143; its cause is unknown. That wrapper
result is retained separately rather than rewritten as exit zero. No extraction
was restarted. Each of the three subsequent read-only scorers exited zero.

Of the 761 transcribed pages, curator `keeppages` directives exclude 86; all
675 included pages were scored. Median prose-token coverage is 0.9448; taxon
coverage averages 0.8902 over 229 eligible pages. The scorer classifies 48 pages
as extraction-empty and one as script-missing. These categories require source
review; successful stage execution does not establish correct content.

Physical figure counts give 0.8803 recall / 0.8642 precision; the default MCP
types give 0.8245 / 1.0000. Count agreement cannot establish that the right
objects were detected. Caption identity binding gives 0.5757 recall / 0.9817
precision against 839 typed identities. Panel exact-set rate is 0.6735. The
receipt and summaries preserve strata and denominators, including weak results.

`receipt.json` retains input/code/report identities, per-document artifact
hashes and compact score summaries. Text files are the actual scorer outputs
with scratch paths normalized. Full per-page reports remain with the build;
original PDFs and transcriptions are not duplicated here. Reproduce using the
receipt commands and the pinned library source manifest; read its
`CROSSCHECK_REPORT.md` before interpreting scores. Final producers, independent
retrieval coverage and full reference-corpus release validation remain pending.
