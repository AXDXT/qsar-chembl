"""Final training uses the full-training selection and serialises reproducibly."""
import json

import joblib
import pytest

from src import gate4_refit as refit
from test_gate4_combined import _inputs


@pytest.mark.parametrize('block', ['combined', 'descriptors'])
def test_freeze_and_resume(tmp_path, block, monkeypatch):
    rows, descriptors = _inputs()
    report = {'identity_sha256': 'search', 'final_fpSize': 1024,
              'final_selection': {'selected_parameters': {'min_samples_leaf': 3}}}
    monkeypatch.setattr(refit.cv, 'RF_SPEC', {**refit.cv.RF_SPEC, 'n_estimators': 3})
    directory = tmp_path / 'model'
    args = (directory, rows, descriptors, report, block, {'split': 'random'}, 1)
    record = refit._freeze(*args)
    assert record == json.loads((directory / 'manifest.json').read_text())
    assert record['reload_verified_rows'] == len(rows)
    model = joblib.load(directory / 'forest.joblib')
    assert model.min_samples_leaf == 3
    assert model.n_features_in_ == descriptors.shape[1] + (1024 if block == 'combined' else 0)
    if block == 'combined':
        assert record['training']['fingerprint_spec']['fpSize'] == 1024
    monkeypatch.setattr(refit.RandomForestRegressor, 'fit', lambda *a: pytest.fail('refitted'))
    assert refit._freeze(*args) == record
    report['identity_sha256'] = 'different'
    with pytest.raises(ValueError, match='identity'):
        refit._freeze(*args)
    report['identity_sha256'] = 'search'
    (directory / 'forest.joblib').write_bytes(b'incomplete')
    with pytest.raises(ValueError, match='forest'):
        refit._freeze(*args)


def test_refuses_test_rows_and_incomplete_search(tmp_path, monkeypatch):
    rows, descriptors = _inputs()
    rows[-1]['partition'] = 'test'
    with pytest.raises(ValueError, match='training-only'):
        refit._freeze(tmp_path / 'bundle', rows, descriptors, {}, 'combined', {}, 1)
    monkeypatch.setattr(refit, 'RESULTS', tmp_path)
    with pytest.raises(FileNotFoundError, match='incomplete'):
        refit.refit('CHEMBL284', 'random', 'combined')
    assert not (tmp_path / 'gate4_models').exists()
