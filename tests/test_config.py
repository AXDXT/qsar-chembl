import os
import random

import numpy as np

from src import config


def test_seed_is_single_fixed_value():
    assert config.SEED == 42


def test_set_global_seed_makes_draws_reproducible():
    config.set_global_seed()
    first = (random.random(), np.random.rand())
    config.set_global_seed()
    assert (random.random(), np.random.rand()) == first
    assert os.environ["PYTHONHASHSEED"] == str(config.SEED)
