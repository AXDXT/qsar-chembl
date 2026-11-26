"""Audit fold-local representation, grouped tuning and resumable combined fits."""
import copy
import json

import numpy as np
import pytest
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
from sklearn.model_selection import KFold

from src import gate4_combined as combined, gate4_fingerprints as fp, nested_cv as cv
from test_baselines import _cached_inputs
from test_gate4_fingerprints import _toy


@pytest.fixture(autouse=True)
def small_forest(monkeypatch):
    monkeypatch.setattr(cv, 'RF_SPEC', {**cv.RF_SPEC, 'n_estimators': 3})


def _inputs():
    rows = [{**r, 'pIC50_median': str(4 + i / 5)} for i, r in enumerate(_toy())]
    x = np.column_stack([np.arange(len(rows)) + j for j in range(12)]).astype(float)
    return rows, x


def _run(path, rows=None, x=None, **kwargs):
    default_rows, default_x = _inputs()
    settings = {'source_checksums': {'data': 'toy'}, 'bootstrap_repeats': 20,
                'grid': {'max_features': ['sqrt', 1.], 'min_samples_leaf': [1]}}
    settings.update(kwargs)
    return combined.run_combined_cv(default_rows if rows is None else rows,
        default_x if x is None else x, path / 'combined', path / 'plan', **settings)


def test_all_fits_use_their_own_representation_and_strictly_grouped_fitting_rows(tmp_path, monkeypatch):
    rows, x = _inputs()
    # Create differing local decisions: the final-fit size must not be reused in CV.
    def evidence(rows):
        return {'measurement_seconds': 0., 'rows': [
            {'inchikey': r['inchikey'], 'identifiers': 2,
             'lost_by_bit_length': {'1024': int(i == 0), '2048': 0, '4096': 0}}
            for i, r in enumerate(rows)]}
    monkeypatch.setattr(fp, '_collision_evidence', evidence)
    actual, calls = cv._fit, []
    def audited(rows, values, labels, fold, parameters, **kwargs):
        assert set(fold['train']).isdisjoint(fold['validation'])
        assert {rows[i]['scaffold'] for i in fold['train']}.isdisjoint(
            rows[i]['scaffold'] for i in fold['validation'])
        expected_bits = 2048 if 0 in fold['train'] else 1024
        assert values.shape == (len(rows), 12 + expected_bits)
        np.testing.assert_array_equal(values[:, :12], x)
        generator = rdFingerprintGenerator.GetMorganGenerator(
            radius=2, includeChirality=False, fpSize=expected_bits)
        for i in fold['train'] + fold['validation']:
            np.testing.assert_array_equal(values[i, 12:], generator.GetFingerprintAsNumPy(
                Chem.MolFromSmiles(rows[i]['canonical_smiles'])))
        calls.append(copy.deepcopy(fold))
        return actual(rows, values, labels, fold, parameters, **kwargs)
    monkeypatch.setattr(cv, '_fit', audited)
    report = _run(tmp_path, rows, x)
    manifest = json.loads((tmp_path / 'combined/manifest.json').read_text())
    assert len(calls) == 65
    assert report['final_fpSize'] == 2048
    for number, fold in enumerate(manifest['outer_plan']):
        for inner in calls[number * 11:number * 11 + 10]:
            assert set(inner['train'] + inner['validation']) == set(fold['train'])
        assert calls[number * 11 + 10] == fold
    assert report['heldout_evaluations'] == 0
    assert report['oof_metrics'] == cv.point_metrics(report['oof_predictions'])
    assert report['oof_bootstrap'] == cv.clustered_intervals(report['oof_predictions'], repeats=20)
    assert len(report['oof_predictions']) == len(rows)
    for prediction in report['oof_predictions']:
        i = next(i for i, r in enumerate(rows) if r['inchikey'] == prediction['inchikey'])
        assert i in manifest['outer_plan'][prediction['outer_fold']]['validation']
    for path in (tmp_path / 'combined').glob('*-fold*.json'):
        record = json.loads(path.read_text())
        assert record['sha256'] == cv._digest(record['result'])
        assert record['result']['fpSize'] in {1024, 2048}


def test_outer_validation_labels_cannot_change_inner_scores_or_selected_model(tmp_path):
    rows, x = _inputs()
    before = _run(tmp_path / 'before', rows, x)
    outer = cv._folds(rows, list(range(len(rows))), 5)
    changed = copy.deepcopy(rows)
    for i in outer[0]['validation']:
        changed[i]['pIC50_median'] = '10000'
    after = _run(tmp_path / 'after', changed, x)
    a, b = before['outer_folds'][0], after['outer_folds'][0]
    assert a['selection']['selected_parameters'] == b['selection']['selected_parameters']
    assert [c['fold_metrics'] for c in a['selection']['candidates']] == [
        c['fold_metrics'] for c in b['selection']['candidates']]
    assert a['metrics'] != b['metrics']


def test_interrupt_resume_preserves_fits_and_completed_rerun_never_featurises(tmp_path, monkeypatch):
    actual, calls = cv._fit, []
    def interrupt(*args, **kwargs):
        calls.append(1)
        if len(calls) == 4:
            raise KeyboardInterrupt()
        return actual(*args, **kwargs)
    monkeypatch.setattr(cv, '_fit', interrupt)
    with pytest.raises(KeyboardInterrupt):
        _run(tmp_path)
    assert len(list((tmp_path / 'combined').glob('outer*-fold*.json'))) == 3
    result = _run(tmp_path, n_jobs=2)
    assert len(calls) == 66
    monkeypatch.setattr(cv, '_fit', lambda *a, **k: pytest.fail('Refitted'))
    monkeypatch.setattr(combined, '_combined_values', lambda *a: pytest.fail('Refeaturised'))
    assert result == _run(tmp_path, n_jobs=3)


@pytest.mark.parametrize('changed', ['labels', 'descriptors', 'source', 'structure', 'grid',
                                     'code', 'trees', 'bootstrap', 'version'])
def test_cache_rejects_changed_inputs(tmp_path, monkeypatch, changed):
    _run(tmp_path)
    rows, x = _inputs()
    kwargs = {}
    if changed == 'labels':
        rows[0]['pIC50_median'] = '99'
    elif changed == 'descriptors':
        x[0, 0] += 1
    elif changed == 'source':
        kwargs['source_checksums'] = {'data': 'different'}
    elif changed == 'structure':
        rows[0]['canonical_smiles'] = 'N'
    elif changed == 'grid':
        kwargs['grid'] = {'max_features': [.1], 'min_samples_leaf': [1]}
    elif changed == 'code':
        monkeypatch.setattr(combined, '_sha256', lambda p: 'changed')
    elif changed == 'trees':
        monkeypatch.setattr(cv, 'RF_SPEC', {**cv.RF_SPEC, 'n_estimators': 4})
    elif changed == 'bootstrap':
        kwargs['bootstrap_repeats'] = 21
    else:
        monkeypatch.setattr(combined.joblib, '__version__', 'changed')
    with pytest.raises(ValueError, match='inputs/settings/code changed'):
        _run(tmp_path, rows, x, **kwargs)


def test_cached_report_tamper_is_rejected(tmp_path):
    _run(tmp_path)
    path = tmp_path / 'combined/report.json'
    record = json.loads(path.read_text())
    record['result']['oof_metrics']['rmse'] = 0
    path.write_text(json.dumps(record))
    with pytest.raises(ValueError, match='checksum/identity'):
        _run(tmp_path)


@pytest.mark.parametrize('invalid', ['test', 'shape', 'nan', 'workers', 'grid', 'source'])
def test_invalid_inputs_never_fit(tmp_path, monkeypatch, invalid):
    rows, x = _inputs()
    kwargs = {}
    if invalid == 'test':
        rows[0]['partition'] = 'test'
    elif invalid == 'shape':
        x = x[:, :-1]
    elif invalid == 'nan':
        x[0, 0] = np.nan
    elif invalid == 'workers':
        kwargs['n_jobs'] = 0
    elif invalid == 'grid':
        kwargs['grid'] = {'n_estimators': [2]}
    else:
        kwargs['source_checksums'] = {}
    monkeypatch.setattr(cv, '_fit', lambda *a, **k: pytest.fail('Fitted invalid input'))
    with pytest.raises(ValueError):
        _run(tmp_path, rows, x, **kwargs)
    assert not (tmp_path / 'combined').exists()


def test_non_grouped_inner_folds_are_rejected(tmp_path, monkeypatch):
    actual, calls = cv.GroupKFold, []
    class Ungrouped:
        def __init__(self, n_splits):
            self.cv = KFold(n_splits=n_splits, shuffle=True, random_state=42)
        def split(self, indices, groups):
            return self.cv.split(indices)
    def splitter(n_splits):
        calls.append(1)
        return actual(n_splits=n_splits) if len(calls) == 1 else Ungrouped(n_splits)
    monkeypatch.setattr(cv, 'GroupKFold', splitter)
    with pytest.raises(AssertionError, match='Scaffold leakage'):
        _run(tmp_path)
    assert not (tmp_path / 'combined').exists()


def test_adapter_only_exposes_saved_training_keys_and_separate_caches(tmp_path, monkeypatch):
    inputs, seen = _cached_inputs(tmp_path), []
    def spy(rows, x, cache_dir, plan_dir, **kwargs):
        assert all(r['partition'] == 'train' for r in rows)
        assert all(Chem.MolToInchiKey(Chem.MolFromSmiles(r['canonical_smiles'])) == r['inchikey']
                   for r in rows)
        assert cache_dir.name == 'combined' and plan_dir.name == 'fingerprint_plan'
        assert cache_dir.parent == plan_dir.parent
        seen.append(rows)
        return {'identity_sha256': 'toy'}
    monkeypatch.setattr(combined, 'run_combined_cv', spy)
    assert set(combined.run_combined_search(*inputs, tmp_path / 'cache')) == {'scaffold', 'random'}
    assert len(seen) == 2
    for dest in (*inputs, tmp_path):
        with pytest.raises(ValueError, match='separate'):
            combined.run_combined_search(*inputs, dest)
    actual = combined._read_csv
    def changed(path):
        result = actual(path)
        with path.open('a') as stream:
            stream.write('\n')
        return result
    monkeypatch.setattr(combined, '_read_csv', changed)
    with pytest.raises(ValueError, match='changed during'):
        combined.run_combined_search(*inputs, tmp_path / 'cache')
    assert len(seen) == 2


def test_misaligned_fingerprint_decision_is_refused_before_fitting(tmp_path, monkeypatch):
    actual = combined.plan_fingerprints
    def misaligned(*args, **kwargs):
        report = copy.deepcopy(actual(*args, **kwargs))
        decision = report['fit_decisions']['outer0-inner-fold0']
        decision['train_indices'] = decision['train_indices'][1:]
        return report
    monkeypatch.setattr(combined, 'plan_fingerprints', misaligned)
    monkeypatch.setattr(cv, '_fit', lambda *a, **k: pytest.fail('Fitted a misaligned plan'))
    with pytest.raises(ValueError, match='does not match fitting rows'):
        _run(tmp_path)
