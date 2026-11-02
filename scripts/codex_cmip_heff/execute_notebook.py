#!/usr/bin/env python3
"""Execute the agent-authored CMIP6 comparison notebook with its Jupyter kernel."""

import json

from jupyter_client import KernelManager

from src.paths import DATA_DIR, PROJECT_ROOT


source = PROJECT_ROOT / "notebooks/codex_heff_comparison_cmip6.ipynb"
OUTPUT = DATA_DIR / "wave_significance/codex_cmip_heff_comparison"
destination = OUTPUT / "codex_heff_comparison_cmip6.individual_models.executed.ipynb"
if destination.exists():
    raise FileExistsError(f"Preserving existing executed notebook: {destination}")
notebook = json.loads(source.read_text())
manager = KernelManager(kernel_name=notebook["metadata"]["kernelspec"]["name"])
manager.start_kernel(cwd=str(PROJECT_ROOT))
client = manager.client()
client.start_channels()
try:
    client.wait_for_ready(timeout=60)
    for number, cell in enumerate(notebook["cells"], start=1):
        if cell["cell_type"] != "code":
            continue
        print(f"Executing cell {number}/{len(notebook['cells'])}", flush=True)
        source_code = cell["source"]
        if isinstance(source_code, list):
            source_code = "".join(source_code)
        message_id = client.execute(source_code, allow_stdin=False)
        outputs = []
        count = None
        error = None
        while True:
            message = client.get_iopub_msg(timeout=900)
            if message["parent_header"].get("msg_id") != message_id:
                continue
            kind = message["msg_type"]
            content = message["content"]
            if kind == "status" and content["execution_state"] == "idle":
                break
            if kind == "execute_input":
                count = content["execution_count"]
            elif kind == "stream":
                outputs.append(dict(output_type="stream", name=content["name"],
                                    text=content["text"]))
            elif kind in ("display_data", "execute_result"):
                output = dict(output_type=kind, data=content["data"],
                              metadata=content.get("metadata", {}))
                if kind == "execute_result":
                    output["execution_count"] = content["execution_count"]
                outputs.append(output)
            elif kind == "clear_output":
                outputs.clear()
            elif kind == "error":
                error = content
                outputs.append(dict(output_type="error", ename=content["ename"],
                                    evalue=content["evalue"], traceback=content["traceback"]))
        cell["execution_count"] = count
        cell["outputs"] = outputs
        print(f"Cell {number} complete; {len(outputs)} outputs", flush=True)
        if error is not None:
            raise RuntimeError(f"Cell {number} failed: {error['ename']}: {error['evalue']}")
finally:
    client.stop_channels()
    manager.shutdown_kernel(now=True)

images = sum("image/png" in output.get("data", {})
             for cell in notebook["cells"] for output in cell.get("outputs", []))
if images != 35:
    raise RuntimeError(f"Expected 35 standalone figures, got {images}")
destination.write_text(json.dumps(notebook, indent=1, ensure_ascii=False) + "\n")
print(f"Executed {len(notebook['cells'])} cells; {images} figures in {destination}", flush=True)
