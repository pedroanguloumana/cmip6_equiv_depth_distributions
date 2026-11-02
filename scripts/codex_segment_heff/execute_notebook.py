#!/usr/bin/env python3
"""Execute the segment-depth comparison from audited compact products."""

import nbformat
from nbclient import NotebookClient

from src.paths import DATA_DIR, PROJECT_ROOT


source = PROJECT_ROOT / "notebooks/codex_segment_heff.ipynb"
destination = DATA_DIR / "wave_significance/codex_segment_heff/codex_segment_heff.executed.ipynb"
if destination.exists():
    raise FileExistsError(f"Preserving the executed comparison: {destination}")
notebook = nbformat.read(source, as_version=4)
client = NotebookClient(notebook, timeout=900, kernel_name="python3",
                        resources={"metadata": {"path": str(PROJECT_ROOT)}})
client.execute()
images = sum("image/png" in output.get("data", {})
             for cell in notebook.cells for output in cell.get("outputs", []))
if images != 17:
    raise RuntimeError(f"Expected 17 standalone figures, found {images}.")
nbformat.write(notebook, destination)
print(f"Executed {len(notebook.cells)} cells and {images} figures: {destination}")
