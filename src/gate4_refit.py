"""Refit a selected Gate 4 forest on its saved training pool and freeze it.

Run: python -m src.gate4_refit CHEMBL284 --split scaffold --block descriptors
Requires the completed training search. Outputs go to results/gate4_models/.
"""
import argparse
import json
from pathlib import Path
import tempfile

import joblib
import numpy as np
from sklearn.ensemble import RandomForestRegressor

from src import gate4_combined as combined
from src import gate4_training as training
from src import nested_cv as cv
from src.config import DATA_CURATED, DATA_FEATURES, DATA_RAW, RESULTS
from src.features import DESCRIPTOR_NAMES
from src.retrieve import _sha256


def _check_saved(directory, expected):
    record = json.loads((directory / 'manifest.json').read_text())
    if record['training'] != expected or record['model_sha256'] != _sha256(directory / 'forest.joblib'):
        raise ValueError('Saved forest or its training identity differs')
    return record


def _freeze(directory, rows, descriptors, report, block, training_id, n_jobs):
    if not rows or any(row['partition'] != 'train' for row in rows):
        raise ValueError('Final refit requires a nonempty training-only pool')
    parameters = {**cv.RF_SPEC, **report['final_selection']['selected_parameters']}
    bits = report['final_fpSize'] if block == 'combined' else None
    training_id = {**training_id, 'feature_block': block, 'search_identity': report['identity_sha256'],
                   'report_sha256': cv._digest(report), 'estimator': parameters,
                   'descriptor_columns': list(DESCRIPTOR_NAMES), 'fpSize': bits,
                   'fingerprint_spec': {**combined.FINGERPRINT_SPEC, 'fpSize': bits} if bits else None,
                   'training_rows': len(rows), 'heldout_evaluations': 0}
    if directory.exists():
        return _check_saved(directory, training_id)
    directory.parent.mkdir(parents=True, exist_ok=True)
    values = combined._combined_values(rows, descriptors, bits) if bits else descriptors
    labels = np.array([float(row['pIC50_median']) for row in rows])
    # Publish the manifest and model together, only after binary round-trip validation.
    with tempfile.TemporaryDirectory(dir=directory.parent) as scratch:
        staged = Path(scratch) / 'bundle'
        staged.mkdir()
        model = RandomForestRegressor(**{**parameters, 'n_jobs': n_jobs}).fit(values, labels)
        model.n_jobs = 1  # Deterministic prediction accumulation and portable reload.
        path = staged / 'forest.joblib'
        joblib.dump(model, path, compress=3)
        predicted = model.predict(values)
        if not np.array_equal(predicted, joblib.load(path).predict(values)):
            raise ValueError('Reloaded forest failed training prediction parity')
        record = {'training': training_id, 'model_sha256': _sha256(path),
                  'reload_verified_rows': len(predicted), 'fit_n_jobs': n_jobs}
        cv._atomic_json(staged / 'manifest.json', record)
        try:
            staged.rename(directory)
        except OSError:
            if not directory.is_dir():
                raise
            return _check_saved(directory, training_id)
    return record


def refit(target, split, block, n_jobs=1):
    if split not in ('scaffold', 'random') or block not in ('combined', 'descriptors') or n_jobs < 1:
        raise ValueError('Choose a saved split, feature block and positive worker count')
    target = target.strip().upper()
    if not target.startswith('CHEMBL') or not target[6:].isdigit():
        raise ValueError('Expected a ChEMBL target ID')
    cache = RESULTS / 'gate4_training' / target
    if not (cache / split / block / 'report.json').is_file():
        raise FileNotFoundError(f'Training search still incomplete for {split}/{block}')
    sources = (DATA_RAW / target, DATA_CURATED / target, RESULTS / 'splits' / target,
               DATA_FEATURES / target)
    inputs = training.load_training_inputs(*sources)
    if inputs['target_chembl_id'] != target:
        raise ValueError('Requested target differs from verified inputs')
    pool = inputs['splits'][split]
    rows, descriptors = pool['rows'], pool['descriptors']
    source_checksums = pool['source_checksums']
    if block == 'descriptors':
        report = cv.run_nested_cv(rows, descriptors, cache / split / block,
            source_checksums=source_checksums, feature_spec={'block': block,
            'columns': list(DESCRIPTOR_NAMES), 'saved_split': split,
            'preprocessing': 'None; stateless descriptors'})
    else:
        path = sources[1] / 'compounds.csv'
        structures = {r['inchikey']: r['canonical_smiles'] for r in training._read_csv(path)}
        if _sha256(path) != source_checksums['curated_compounds']:
            raise ValueError('Curated structures changed during final refit preparation')
        rows = [{**r, 'canonical_smiles': structures[r['inchikey']]} for r in rows]
        report = combined.run_combined_cv(rows, descriptors, cache / split / block,
            cache / split / 'fingerprint_plan', source_checksums=source_checksums)
    return _freeze(RESULTS / 'gate4_models' / target / split / block, rows, descriptors,
                   report, block, {'target': target, 'split': split}, n_jobs)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('target', nargs='?', default='CHEMBL284')
    parser.add_argument('--split', choices=('scaffold', 'random'), required=True)
    parser.add_argument('--block', choices=('descriptors', 'combined'), required=True)
    parser.add_argument('--n-jobs', type=int, default=1)
    args = parser.parse_args()
    record = refit(args.target, args.split, args.block, args.n_jobs)
    print(f"{args.split}/{args.block}: verified {record['reload_verified_rows']} training rows")


if __name__ == '__main__':
    main()
