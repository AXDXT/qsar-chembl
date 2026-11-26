"""Frozen models are linked to their training search, without fitting or scoring."""
import json

import joblib
import pytest

from src import gate4_combined as combined, gate4_model_io as model_io
from src import gate4_refit as refit, nested_cv as cv
from test_gate4_combined import _inputs


@pytest.fixture(params=['descriptors', 'combined'])
def bundle(tmp_path, monkeypatch, request):
    monkeypatch.setattr(cv, 'RF_SPEC', {**cv.RF_SPEC, 'n_estimators': 3})
    rows, descriptors = _inputs()
    block = request.param
    search = tmp_path / 'search'
    settings = {'source_checksums': {'data': 'toy'}, 'folds': 2,
                'grid': {'max_features': [1.0], 'min_samples_leaf': [1]},
                'bootstrap_repeats': 2}
    if block == 'combined':
        report = combined.run_combined_cv(rows, descriptors, search, tmp_path / 'plan', **settings)
    else:
        report = cv.run_nested_cv(rows, descriptors, search,
                                 feature_spec={'block': block}, **settings)
    directory = tmp_path / 'bundle'
    refit._freeze(directory, rows, descriptors, report, block,
                  {'target': 'CHEMBL284', 'split': 'scaffold'}, 1)
    # Loading must never train, generate features or predict on either partition.
    def forbidden(*args, **kwargs):
        pytest.fail('Loader fitted, featurised or predicted')
    monkeypatch.setattr(cv.RandomForestRegressor, 'fit', forbidden)
    monkeypatch.setattr(cv.RandomForestRegressor, 'predict', forbidden)
    monkeypatch.setattr(combined, '_combined_values', forbidden)
    return directory, search, rows, descriptors, block


def load(bundle, **kwargs):
    directory, search, rows, descriptors, block = bundle
    return model_io._load_bundle(directory, search, rows, descriptors,
        kwargs.pop('source_checksums', {'data': 'toy'}), target='CHEMBL284',
        split=kwargs.pop('split', 'scaffold'), block=block, **kwargs)


def test_load_exposes_search_metadata_without_retraining_or_scoring(bundle):
    directory, search, _, _, _ = bundle
    before = {path: path.read_bytes() for parent in (directory, search)
              for path in parent.iterdir() if path.is_file()}
    model, metadata = load(bundle)
    manifest = json.loads((search / 'manifest.json').read_text())
    assert metadata['versions'] == manifest['versions']
    assert metadata['input_checksums'] == {'data': 'toy'}
    assert metadata['seed'] == model.random_state == cv.SEED
    assert metadata['training']['heldout_evaluations'] == 0
    assert metadata['search_code_sha256'] == manifest['code_sha256']
    assert before == {path: path.read_bytes() for path in before}


@pytest.mark.parametrize('changed', ['binary', 'model_manifest', 'report', 'search_manifest',
                                     'source', 'code', 'version', 'rows', 'features', 'split',
                                     'estimator'])
def test_drift_is_rejected(bundle, monkeypatch, changed):
    directory, search, rows, descriptors, block = bundle
    kwargs = {}
    if changed == 'binary':
        (directory / 'forest.joblib').write_bytes(b'incomplete')
    elif changed == 'model_manifest':
        path = directory / 'manifest.json'
        record = json.loads(path.read_text())
        record['training']['training_rows'] += 1
        path.write_text(json.dumps(record))
    elif changed in ('report', 'search_manifest'):
        path = search / ('report.json' if changed == 'report' else 'manifest.json')
        record = json.loads(path.read_text())
        if changed == 'report':
            record['result']['final_selection']['selected_parameters']['min_samples_leaf'] = 99
        else:
            record['seed'] += 1
        path.write_text(json.dumps(record))
    elif changed == 'source':
        kwargs['source_checksums'] = {'data': 'changed'}
    elif changed == 'code':
        monkeypatch.setattr(model_io, '_sha256', lambda path: 'changed')
    elif changed == 'version':
        monkeypatch.setattr(model_io.joblib, '__version__', 'changed')
    elif changed == 'rows':
        rows[0]['partition'] = 'test'
    elif changed == 'features':
        descriptors[0, 0] += 1
    elif changed == 'split':
        kwargs['split'] = 'random'
    else:
        # Even a binary with an updated checksum must match declared parameters.
        path = directory / 'forest.joblib'
        model = joblib.load(path)
        model.min_samples_leaf = 99
        joblib.dump(model, path)
        record_path = directory / 'manifest.json'
        record = json.loads(record_path.read_text())
        record['model_sha256'] = refit._sha256(path)
        record_path.write_text(json.dumps(record))
    if changed != 'estimator':
        monkeypatch.setattr(model_io.joblib, 'load', lambda *a, **k: pytest.fail('Loaded invalid binary'))
    with pytest.raises(ValueError):
        load(bundle, **kwargs)


def test_public_loader_requires_complete_artifacts_and_valid_names(tmp_path, monkeypatch):
    monkeypatch.setattr(model_io, 'RESULTS', tmp_path)
    monkeypatch.setattr(model_io.training, 'load_training_inputs',
                        lambda *a: pytest.fail('Prepared missing model'))
    with pytest.raises(FileNotFoundError, match='missing'):
        model_io.load_model('CHEMBL284', 'scaffold', 'combined')
    for args in [('CHEMBL284/../../', 'scaffold', 'combined'),
                 ('CHEMBL284', 'other', 'combined'), ('CHEMBL284', 'random', 'other')]:
        with pytest.raises(ValueError, match='Choose'):
            model_io.load_model(*args)
    assert list(tmp_path.iterdir()) == []


@pytest.mark.parametrize('split', ['scaffold', 'random'])
def test_saved_descriptor_models_load_read_only_with_current_verified_inputs(split, monkeypatch):
    def forbidden(*args, **kwargs):
        pytest.fail('Loading a production model trained or scored it')
    monkeypatch.setattr(cv.RandomForestRegressor, 'fit', forbidden)
    monkeypatch.setattr(cv.RandomForestRegressor, 'predict', forbidden)
    model, metadata = model_io.load_model('CHEMBL284', split, 'descriptors')
    assert model.n_features_in_ == len(model_io.DESCRIPTOR_NAMES)
    assert metadata['training']['split'] == split
    assert metadata['training']['training_rows'] == metadata['reload_verified_rows']
    assert metadata['input_checksums']['curated_compounds']
