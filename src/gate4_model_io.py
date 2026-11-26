"""Read-only loading of frozen Gate 4 forests with verified search provenance.

Only load this repository's trusted joblib artifacts. Checksums detect drift;
they do not make pickle-based files from an untrusted source safe to open.
No fitting, prediction, parameter selection or held-out evaluation occurs here.
"""
from __future__ import annotations

import hashlib
import io
import json
from pathlib import Path
import platform

import joblib
import numpy as np
from rdkit import rdBase
import sklearn
from sklearn.ensemble import RandomForestRegressor

from src import gate4_training as training
from src import nested_cv as cv
from src.config import DATA_CURATED, DATA_FEATURES, DATA_RAW, RESULTS
from src.features import DESCRIPTOR_NAMES, _read_csv
from src.retrieve import _sha256
from src.splits import FINGERPRINT_SPEC


def _load_bundle(directory, search_dir, rows, descriptors, source_checksums, *,
                 target, split, block):
    """Validate links to the completed search before deserialising model bytes."""
    directory, search_dir = Path(directory), Path(search_dir)
    manifest = json.loads((search_dir / 'manifest.json').read_text())
    wrapper = json.loads((search_dir / 'report.json').read_text())
    report = wrapper['result']
    identity = cv._digest(manifest)
    if (wrapper.get('name') != 'report' or wrapper.get('identity') != identity
            or wrapper.get('sha256') != cv._digest(report)
            or report.get('identity_sha256') != identity):
        raise ValueError('Frozen model search checksum/identity mismatch')
    versions = {'python': platform.python_version(), 'numpy': np.__version__,
                'sklearn': sklearn.__version__, 'rdkit': rdBase.rdkitVersion,
                'joblib': joblib.__version__}
    if manifest['versions'] != versions:
        raise ValueError('Frozen model library versions changed')
    if (manifest['source_checksums'] != source_checksums
            or any(_sha256(Path(training.__file__).with_name(name)) != digest
                   for name, digest in manifest['code_sha256'].items())):
        raise ValueError('Frozen model input or search code changed')
    columns = ('inchikey', 'scaffold', 'pIC50_median', 'partition')
    expected_rows = ([{key: row[key] for key in columns} for row in rows]
                     if block == 'descriptors' else rows)
    descriptors = np.asarray(descriptors, dtype=np.float64)
    feature_hash = 'features_sha256' if block == 'descriptors' else 'descriptor_sha256'
    if (not rows or any(row['partition'] != 'train' for row in rows)
            or manifest['rows'] != expected_rows
            or descriptors.shape != (len(rows), len(DESCRIPTOR_NAMES))
            or not np.isfinite(descriptors).all()
            or manifest[feature_hash] != hashlib.sha256(descriptors.tobytes(order='C')).hexdigest()):
        raise ValueError('Frozen model training rows or descriptors changed')
    bits = report['final_fpSize'] if block == 'combined' else None
    parameters = {**cv.RF_SPEC, **report['final_selection']['selected_parameters']}
    expected = {'target': target, 'split': split, 'feature_block': block,
                'search_identity': identity, 'report_sha256': cv._digest(report),
                'estimator': parameters, 'descriptor_columns': list(DESCRIPTOR_NAMES),
                'fpSize': bits,
                'fingerprint_spec': {**FINGERPRINT_SPEC, 'fpSize': bits} if bits else None,
                'training_rows': len(rows), 'heldout_evaluations': 0}
    record = json.loads((directory / 'manifest.json').read_text())
    if (record['training'] != expected or report['heldout_evaluations'] != 0
            or record['reload_verified_rows'] != len(rows)):
        raise ValueError('Frozen model metadata differs from its completed training search')
    # Hash and load the same byte snapshot, avoiding a second path read at unpickle.
    binary = (directory / 'forest.joblib').read_bytes()
    if hashlib.sha256(binary).hexdigest() != record['model_sha256']:
        raise ValueError('Frozen model binary checksum mismatch')
    model = joblib.load(io.BytesIO(binary))
    if (not isinstance(model, RandomForestRegressor)
            or model.get_params() != RandomForestRegressor(**parameters).get_params()
            or model.n_features_in_ != len(DESCRIPTOR_NAMES) + (bits or 0)
            or len(model.estimators_) != parameters['n_estimators']):
        raise ValueError('Frozen forest does not match its declared estimator/features')
    metadata = {**record, 'versions': manifest['versions'],
                'input_checksums': manifest['source_checksums'],
                'search_code_sha256': manifest['code_sha256'], 'seed': manifest['seed']}
    return model, metadata


def load_model(target, split, block):
    """Return (forest, provenance) after checking current saved training inputs.

    Raises on missing/incomplete artifacts or drift; never resumes a search or
    refits a missing model. Train/test file checks verify provenance only. Labels
    returned by the adapter and passed to the loader are training labels only.
    """
    target = target.strip().upper()
    if (not target.startswith('CHEMBL') or not target[6:].isdigit()
            or split not in ('scaffold', 'random') or block not in ('descriptors', 'combined')):
        raise ValueError('Choose a ChEMBL target, saved split and feature block')
    directory = RESULTS / 'gate4_models' / target / split / block
    search_dir = RESULTS / 'gate4_training' / target / split / block
    for path in (directory / 'manifest.json', directory / 'forest.joblib',
                 search_dir / 'manifest.json', search_dir / 'report.json'):
        if not path.is_file():
            raise FileNotFoundError(f'Frozen model or completed search missing: {path}')
    inputs = training.load_training_inputs(DATA_RAW / target, DATA_CURATED / target,
        RESULTS / 'splits' / target, DATA_FEATURES / target)
    if inputs['target_chembl_id'] != target:
        raise ValueError('Requested target differs from verified inputs')
    pool = inputs['splits'][split]
    rows = pool['rows']
    if block == 'combined':
        path = DATA_CURATED / target / 'compounds.csv'
        structures = {row['inchikey']: row['canonical_smiles'] for row in _read_csv(path)}
        if _sha256(path) != pool['source_checksums']['curated_compounds']:
            raise ValueError('Curated structures changed during model loading')
        rows = [{**row, 'canonical_smiles': structures[row['inchikey']]} for row in rows]
    return _load_bundle(directory, search_dir, rows, pool['descriptors'],
                        pool['source_checksums'], target=target, split=split, block=block)
