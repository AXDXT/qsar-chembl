"""Execute the trusted Gate 2 notebook offline without a kernel/server listener.

Usage (repo root): .venv/bin/python scripts/run_gate2_notebook.py [--output PATH]
Runs every code cell in one IPython namespace; fails on the first cell error and
only writes a complete, validated notebook. No output notebook is written on failure.
"""
from __future__ import annotations

import argparse
import os
from pathlib import Path

import nbformat
from IPython.core.interactiveshell import InteractiveShell
from IPython.utils.capture import capture_output


ROOT = Path(__file__).resolve().parents[1]
NOTEBOOK = ROOT / "notebooks/gate2_dpp4_curation_and_leakage.ipynb"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=NOTEBOOK)
    args = parser.parse_args()
    output = args.output.resolve()
    notebook = nbformat.read(NOTEBOOK, as_version=4)
    nbformat.validate(notebook)
    os.chdir(ROOT)
    shell = InteractiveShell.instance()
    count = 0
    for index, cell in enumerate(notebook.cells, start=1):
        if cell.cell_type != "code":
            continue
        count += 1
        print(f"Executing code cell {count} (notebook cell {index})", flush=True)
        cell.outputs = []
        with capture_output() as captured:
            result = shell.run_cell(cell.source, store_history=True)
        # Raise outside capture so failures remain visible and return a nonzero exit.
        result.raise_error()
        cell.execution_count = count
        for channel in ("stdout", "stderr"):
            value = getattr(captured, channel)
            if value:
                cell.outputs.append(nbformat.v4.new_output("stream", name=channel, text=value))
        for rich in captured.outputs:
            cell.outputs.append(nbformat.v4.new_output(
                "display_data", data=rich.data, metadata=rich.metadata))
    nbformat.validate(notebook)
    nbformat.write(notebook, output)
    print(f"PASS: executed all {count} code cells; saved {output.relative_to(ROOT) if output.is_relative_to(ROOT) else output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
