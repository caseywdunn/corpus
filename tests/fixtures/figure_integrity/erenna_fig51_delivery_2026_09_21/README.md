Source Figure 51 in Pugh and Haddock (2016), physical page 44, passes the
bounded #329 crop-to-delivery acceptance at `eb50d75`. The complete printed
`E. insidiator` label, all five specimen images and all five scale bars survive;
the neighboring caption and prose are excluded. The source page, legacy clipped
image and returned image were visually compared by a Codex team agent. A second
team agent separately inspected the returned image and source edge; neither
review is claimed as human or independent domain review.

`receipt.json` retains the actual legacy figure record, a later saved Docling
picture observation with matching bounds, original PDF/artifact hashes, repaired
geometry and returned PNG hash. Production `render_figures` completed the narrowly
clipped bounds using the PDF's embedded image extent. Production `package`,
`CorpusIndex.load` and `MCPServer.call_tool("get_figure_image", profile="report")`
then delivered identical bytes: 2871 × 1133 pixels, SHA-256
`5ca44a19da9633bdabb928a5719934178e5e85213c55d271e39cb6f56dfe1108`.
The source inputs and bundle were unchanged by delivery. No PDFs or image corpus
are committed here.

To repeat from the repository root in the installed corpus environment, supply
the original library PDF and a new scratch directory:

```sh
PYTHONPATH=. python - SOURCE_PDF NEW_SCRATCH_DIR <<'PY'
import asyncio, base64, hashlib, json, sys
from pathlib import Path
from pipeline.figures import render_figures
from mcpsrv.bundle import package
from mcpsrv.indexes import CorpusIndex
from mcpsrv import app
import mcpsrv.tools.figures
receipt = json.loads(Path('tests/fixtures/figure_integrity/erenna_fig51_delivery_2026_09_21/receipt.json').read_text())
source, scratch = map(Path, sys.argv[1:])
assert not scratch.exists()
with source.open('rb') as f:
    assert hashlib.file_digest(f, 'sha256').hexdigest() == receipt['source_pdf']['sha256']
paper = scratch / 'build/documents/da370f1e0434'
figures = paper / 'figures'
figures.mkdir(parents=True)
figure = receipt['original_figure']
figure['file_path'] = str((figures / figure['filename']).resolve())
assert render_figures(source, [figure], figures, native=True)['rendered'] == 1
(paper / 'figures.json').write_text(json.dumps({'figures': [figure], 'total_figures': 1}))
package(scratch / 'build', scratch / 'bundle', 'issue-329-replay', False, False)
index = CorpusIndex(scratch / 'bundle'); assert index.load() == 1
app.set_index(index)
result = asyncio.run(app.mcp.call_tool('get_figure_image', receipt['served']['arguments']))
assert not result.is_error and result.content[0].type == 'image'
png = base64.b64decode(result.content[0].data)
assert png == (figures / 'fig_51.png').read_bytes()
(scratch / 'served-whole.png').write_bytes(png)
print(hashlib.sha256(png).hexdigest())
PY
python -m pytest -q tests/test_figure_raster_edges.py
```

The existing 14 controls passed. Four rotations, two nonzero-cropbox cases and
the nearby-prose veto also received source/crop visual comparison. The veto
deliberately leaves its artificial clipped edge unchanged because expansion
would include additional prose. These controls remain separate from the real
source case; panel ROI conversion remains separately tested.

The scratch bundle is only an image-delivery probe: missing bibliography and
taxon databases produce expected warnings, and no embeddings/models are loaded.
This replay does not replace fresh layout extraction, a full candidate build,
external deployment acceptance, or corpus-wide clipping/recall measurement.
Original local evidence and the executed replay script were retained under
`/tmp/corpus-issue329-delivery/`; those paths are historical provenance, not
dependencies of the reproduction above. PNG bytes can vary with PyMuPDF versions;
the recorded run used PyMuPDF 1.28.0, and source geometry/visual contents are the
primary acceptance evidence.
