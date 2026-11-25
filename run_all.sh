#!/usr/bin/env bash
# One-command pipeline for the steps that exist so far (Gate 1).
# Usage: ./run_all.sh [TARGET_ID]      (default CHEMBL284, provisional recommendation)
# A target that already has data/raw/<TARGET_ID>/manifest.json is NOT re-pulled:
# data/raw is immutable. Use `python -m src.retrieve <ID> --force` to replace a pull deliberately.
set -euo pipefail
cd "$(dirname "$0")"

PY=".venv/bin/python"
TARGET="${1:-CHEMBL284}"

# -x skips macOS "._*" AppleDouble files that the exFAT drive creates next to every source file.
"$PY" -m compileall -q -x '/\._' src
"$PY" -m pytest -q

if [[ -f "results/licence/manifest.json" ]]; then
    echo "Gate 1 licence snapshot exists; skipping live fetch."
else
    "$PY" -m src.provenance
fi

if [[ -f "results/candidate_counts/candidate_counts.csv" ]]; then
    echo "Gate 1 candidate-count snapshot exists; skipping live recount."
else
    "$PY" -m src.candidates
fi

if [[ -f "data/raw/${TARGET}/manifest.json" ]]; then
    echo "data/raw/${TARGET} already pulled; skipping retrieval."
else
    "$PY" -m src.retrieve "$TARGET"
fi

"$PY" -m src.profile "$TARGET"
