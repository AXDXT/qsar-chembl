"""Single source of truth for paths and the global random seed.

Every module that does anything stochastic must call ``set_global_seed()``
(or pass ``SEED`` explicitly as ``random_state``). Nothing else in the repo
should hard-code a seed.
"""
from __future__ import annotations

import os
import random
from pathlib import Path

import numpy as np

SEED = 42

ROOT = Path(__file__).resolve().parent.parent
DATA_RAW = ROOT / "data" / "raw"          # immutable: unmodified API responses only
DATA_CURATED = ROOT / "data" / "curated"  # derived from data/raw by code in src/
DATA_FEATURES = ROOT / "data" / "features"  # Gate 3 feature cache, derived from data/curated
RESULTS = ROOT / "results"

CHEMBL_API_BASE = "https://www.ebi.ac.uk/chembl/api/data"


def set_global_seed(seed: int = SEED) -> int:
    """Seed Python, NumPy and the hash seed for child processes."""
    os.environ["PYTHONHASHSEED"] = str(seed)
    random.seed(seed)
    np.random.seed(seed)
    return seed
