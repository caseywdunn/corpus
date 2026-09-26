# September 25 Qwen source review

**Acceptance fails: one of nine selected target crops passes; eight fail.**
This is a purposive pilot result, not a corpus-wide error rate. Both #305 and
#342 remain open. No production model, prompt or cropping behavior was changed
to obtain these results.

[Machine-readable judgments and artifact hashes](review.json) record direct
inspection of each original page, fresh ROI overlay and target crop. Inference
ran in SLURM job `27479041`; this review uses its retained extraction rasters,
not a fresh full extraction. The three named issue cases and six additional
deterministic controls are distinct from the unavailable original six-control
manifest. Preservation of that original known-good population is unestablished.

| Source and target | Visual result | Evidence |
| --- | --- | --- |
| Siebert 2013, physical p. 5, B | Pass | Complete colony and B label; no neighboring-panel substitution. |
| Deevey–Brooks 1971, p. 12, A | Fail | Most of the time series, axis values and units are outside the crop; it cuts through the peak. |
| Hays 2018, p. 6, A | Fail | Lower birds and penguin are clipped by the A/B boundary. |
| Rees 1966, p. 5, A | Fail | Adult/arrow/context clipped and A label lost; extra whole-image figure ROI. |
| Deevey–Brooks 1971, p. 9, A | Fail | Narrow strip loses A label and most data, and includes B content; D/E absent from response. |
| Park–Lee 2022, p. 6, A | Fail | Specimen margins, A label and scale bar clipped. |
| Castriota 2017, p. 3, C | Fail | Float and ruler omitted; only right-hand tentacles retained. Input raster already contains only C although the caption requests A/B/C. |
| MacBride 1914, p. 5, A | Fail | Oval embryo clipped and A label omitted. |
| Gamero 2015, p. 6, A | Fail | Upper sampling symbols clipped. Full overlay also splits B's legend and drops shared scale/coordinates. |

## Coordinate evidence versus content

All nine outputs have in-bounds boxes. Their recorded frames match the actual
generation-call processor dimensions, and conversion of recorded model boxes
to raster boxes agrees within one pixel. For the named Deevey–Brooks and Hays
failures, the actual generation frame is byte-identical to the inspected frame:
the complete content was available to the model. The residual crop failures
therefore persist after the coordinate conversion correction; they are not
evidence that the old double-resize arithmetic remains in these captures.

The pilot's original `numeric_checks_passed` field is preserved for provenance,
but its name was too broad: it required `pass3_status == completed`. Two of its
four failures returned `partial_vision`; Rees and Park returned
`completed_compound` and also had valid bounds/frame provenance. Conversely,
four of the five original numeric passes fail visual review. Neither a complete
status nor a bounded rectangle proves the right scientific content was cropped.

Raw pages, processor frames, responses, overlays and crops are retained under
[`scratch/v15-bouchet-20260925/vision-pilot-v2`](../../../scratch/v15-bouchet-20260925/vision-pilot-v2).
Each case's `review.json` now holds its completed judgment. The original capture
summary remains preserved separately when its review rollup is updated.

## Remaining acceptance

Correct the residual target-box/content behavior and validate against this
same frozen selection, including the passing Siebert control. Keep upstream
raster/caption mismatches explicit, especially Castriota. Recovered original
controls are needed before claiming reproduction or preservation of the earlier
pilot. Broader inference or a full-library vision run is not established here.
