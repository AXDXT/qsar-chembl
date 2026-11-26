"""Gate 4 combined descriptors+Morgan search; training rows only.

Run: .venv/bin/python -m src.gate4_combined CHEMBL284 --n-jobs 4
Each forest uses its own fitting-row collision decision. No final model is fitted
and no saved test row is featurised, predicted or scored by this module.

The orchestration is separate from nested_cv so completed descriptor-control and
fingerprint-plan cache identities remain valid. Grouping, fitting, atomic caches,
metrics and uncertainty use the existing Gate 2/3/4 helpers unchanged.
"""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import platform
from time import perf_counter

import joblib
import numpy as np
from rdkit import Chem, rdBase
from rdkit.Chem import rdFingerprintGenerator
import sklearn
from sklearn.model_selection import ParameterGrid

from src import nested_cv as cv
from src.config import DATA_CURATED, DATA_FEATURES, DATA_RAW, RESULTS
from src.features import DESCRIPTOR_NAMES, _read_csv
from src.gate4_fingerprints import plan_fingerprints
from src.gate4_training import load_training_inputs
from src.retrieve import _sha256
from src.splits import FINGERPRINT_SPEC


def _combined_values(rows, descriptors, bits):
    """Stateless features; the caller supplies the fit-local selected bit length."""
    generator = rdFingerprintGenerator.GetMorganGenerator(
        radius=FINGERPRINT_SPEC['radius'], includeChirality=FINGERPRINT_SPEC['includeChirality'],
        fpSize=bits)
    fingerprints = []
    for row in rows:
        mol = Chem.MolFromSmiles(row['canonical_smiles'])
        if mol is None or Chem.MolToInchiKey(mol) != row['inchikey']:
            raise ValueError('Curated SMILES/InChIKey mismatch')
        fingerprints.append(generator.GetFingerprintAsNumPy(mol))
    return np.column_stack((descriptors, np.asarray(fingerprints, dtype=np.uint8)))


def run_combined_cv(rows, descriptors, cache_dir: Path, plan_dir: Path, *,
                    source_checksums: dict, folds: int = 5, grid: dict | None = None,
                    bootstrap_repeats: int = cv.BOOTSTRAP_REPEATS, n_jobs: int = 1) -> dict:
    """Resume fold-local combined search, using the descriptor control's grid/rule.

    Bit-length planning independently verifies its cached inputs/code and groups.
    Its report digest pins every decision here. Inner validation and outer test
    rows cannot enter a fit's representation decision or hyperparameter selection.
    Computing stateless features for training-pool validation rows is permitted;
    no fitted preprocessing or column selection uses those rows.
    """
    if not rows or any(r.get('partition') != 'train' for r in rows):
        raise ValueError('Combined CV accepts saved training rows only')
    descriptors = np.asarray(descriptors, dtype=np.float64)
    y = np.array([float(r['pIC50_median']) for r in rows])
    if (descriptors.shape != (len(rows), len(DESCRIPTOR_NAMES))
            or not np.isfinite(descriptors).all() or not np.isfinite(y).all()):
        raise ValueError('Require aligned finite descriptor columns and labels')
    if not isinstance(n_jobs, int) or isinstance(n_jobs, bool) or n_jobs < 1:
        raise ValueError('n_jobs must be a positive worker count')
    if not isinstance(bootstrap_repeats, int) or bootstrap_repeats < 2:
        raise ValueError('Need at least two bootstrap draws')
    grid = cv.DEFAULT_GRID if grid is None else grid
    if set(grid) != {'max_features', 'min_samples_leaf'}:
        raise ValueError('Search only max_features and min_samples_leaf')
    candidates = list(ParameterGrid(grid))
    if not candidates:
        raise ValueError('Require a nonempty grid')
    for parameters in candidates:
        cv.RandomForestRegressor(**{**cv.RF_SPEC, **parameters})._validate_params()
    indices = list(range(len(rows)))
    outer = cv._folds(rows, indices, folds)
    inner = [cv._folds(rows, f['train'], folds) for f in outer]
    final = cv._folds(rows, indices, folds)
    if Path(cache_dir).resolve() == Path(plan_dir).resolve():
        raise ValueError('Search and fingerprint plan need separate cache directories')
    plan = plan_fingerprints(rows, plan_dir, source_checksums=source_checksums, folds=folds)
    identity = json.loads(cv._json({
        'schema': 1, 'rows': rows,
        'descriptor_sha256': hashlib.sha256(descriptors.tobytes(order='C')).hexdigest(),
        'descriptor_columns': DESCRIPTOR_NAMES, 'source_checksums': source_checksums,
        'fingerprint_spec': FINGERPRINT_SPEC, 'fingerprint_plan_sha256': cv._digest(plan),
        'candidates': candidates, 'estimator': cv.RandomForestRegressor(**cv.RF_SPEC).get_params(),
        'versions': {'python': platform.python_version(), 'numpy': np.__version__,
                     'sklearn': sklearn.__version__, 'rdkit': rdBase.rdkitVersion,
                     'joblib': joblib.__version__},
        'code_sha256': {name: _sha256(Path(__file__).with_name(name)) for name in
                        ('gate4_combined.py', 'nested_cv.py', 'gate4_fingerprints.py',
                         'gate4_training.py', 'features.py', 'splits.py', 'leakage_probe.py',
                         'baselines.py', 'retrieve.py', 'config.py')},
        'folds': folds, 'outer_plan': outer, 'inner_plan': inner, 'final_plan': final,
        'bootstrap_repeats': bootstrap_repeats, 'seed': cv.SEED,
    }))
    cache = cv._Cache(cache_dir, identity)
    matrices = {}

    def fit(name, decision_name, fold, parameters):
        decision = plan['fit_decisions'][decision_name]
        if decision['train_indices'] != fold['train']:
            raise ValueError('Fingerprint decision does not match fitting rows')
        bits = decision['selected_fpSize']

        def compute():
            if bits not in matrices:
                matrices[bits] = _combined_values(rows, descriptors, bits)
            result = cv._fit(rows, matrices[bits], y, fold, parameters, n_jobs=n_jobs)
            return {**result, 'fingerprint_decision': decision_name, 'fpSize': bits}
        return cache.get(name, compute)

    def search(prefix, fold_plan):
        scores = []
        for candidate, parameters in enumerate(candidates):
            results = [fit(f'{prefix}-candidate{candidate}-fold{number}',
                           f'{prefix}-fold{number}', fold, parameters)
                       for number, fold in enumerate(fold_plan)]
            scores.append({'parameters': parameters,
                'mean_validation_rmse': float(np.mean([r['metrics']['rmse'] for r in results])),
                'fold_metrics': [r['metrics'] for r in results],
                'fit_predict_seconds': [r['fit_predict_seconds'] for r in results],
                'fold_fpSizes': [r['fpSize'] for r in results]})
        winner = min(range(len(scores)), key=lambda i: scores[i]['mean_validation_rmse'])
        return {'selected_parameters': candidates[winner], 'selected_candidate': winner,
                'candidates': scores}

    def complete():
        started = perf_counter()
        results, predictions = [], []
        for number, fold in enumerate(outer):
            selected = search(f'outer{number}-inner', inner[number])
            evaluation = fit(f'outer{number}-evaluation', f'outer{number}-evaluation',
                             fold, selected['selected_parameters'])
            predictions.extend({**r, 'outer_fold': number} for r in evaluation['predictions'])
            results.append({'outer_fold': number, 'selection': selected,
                            'metrics': evaluation['metrics'], 'fpSize': evaluation['fpSize'],
                            'fit_predict_seconds': evaluation['fit_predict_seconds']})
        final_selection = search('final-inner', final)
        counts = Counter(cv._json(r['selection']['selected_parameters']) for r in results)
        spread = {}
        for metric in ('rmse', 'mae', 'r2'):
            valid = [r['metrics'][metric] for r in results if r['metrics'][metric] is not None]
            spread[metric] = {'values': [r['metrics'][metric] for r in results],
                             'mean': float(np.mean(valid)) if valid else None,
                             'sample_sd': float(np.std(valid, ddof=1)) if len(valid) > 1 else None}
        searches = [r['selection'] for r in results] + [final_selection]
        seconds = sum(sum(c['fit_predict_seconds']) for s in searches for c in s['candidates'])
        seconds += sum(r['fit_predict_seconds'] for r in results)
        return {'identity_sha256': cache.identity, 'feature_block': 'descriptors+Morgan',
            'grid': candidates, 'grid_rationale': cv.GRID_RATIONALE,
            'selection_rule': 'Minimum unweighted mean inner-fold RMSE; ties use logged grid order.',
            'fold_grouping': 'GroupKFold in both loops, using saved Gate 2 exact scaffold labels.',
            'fingerprint_plan_identity': plan['identity_sha256'],
            'fingerprint_plan_sha256': cv._digest(plan),
            'final_fpSize': plan['fit_decisions']['final-fit']['selected_fpSize'],
            'outer_folds': results, 'outer_fold_spread': spread,
            'oof_predictions': predictions, 'oof_metrics': cv.point_metrics(predictions),
            'oof_bootstrap': cv.clustered_intervals(predictions, repeats=bootstrap_repeats),
            'uncertainty': 'Whole-scaffold percentile bootstrap conditions on fixed OOF predictions '
                'and folds; excludes selection/refit/split-seed uncertainty. Fold SD is descriptive, '
                'not a confidence interval: training folds overlap.',
            'selection_stability': {'all_outer_folds_agree': len(counts) == 1,
                'distinct_selections': len(counts), 'counts': [
                    {'parameters': json.loads(k), 'outer_folds': v} for k, v in sorted(counts.items())]},
            'final_selection': final_selection,
            'timing': {'cumulative_fit_predict_seconds': seconds,
                'completion_invocation_seconds': perf_counter() - started,
                'n_jobs': n_jobs, 'search_workers': 1,
                'note': 'Fit times include reused work; report retains original completion timing. '
                        'Each fit records its worker count. Feature generation is in invocation '
                        'time, outside fit/predict time; fingerprint-plan timings are separate.'},
            'heldout_evaluations': 0,
            'disclaimer': 'Literature/patent-reported IC50; no efficacy, safety or suitability claim.'}
    return cache.get('report', complete)


def run_combined_search(raw_dir, curated_dir, split_dir, feature_dir, cache_dir, *,
                        splits=('scaffold', 'random'), n_jobs=1):
    """Verify sources, then join structures onto saved training rows only."""
    if not splits or len(set(splits)) != len(splits) or set(splits) - {'scaffold', 'random'}:
        raise ValueError('Choose unique saved split names: scaffold and/or random')
    dest = Path(cache_dir).resolve()
    for source in (raw_dir, curated_dir, split_dir, feature_dir, DATA_RAW):
        source = Path(source).resolve()
        if dest == source or source in dest.parents or dest in source.parents:
            raise ValueError('Combined caches must be separate from all input trees')
    inputs = load_training_inputs(raw_dir, curated_dir, split_dir, feature_dir)
    path = Path(curated_dir) / 'compounds.csv'
    structures = {r['inchikey']: r['canonical_smiles'] for r in _read_csv(path)}
    if any(_sha256(path) != inputs['splits'][s]['source_checksums']['curated_compounds'] for s in splits):
        raise ValueError('Curated structures changed during combined preparation')
    reports = {}
    for split in splits:
        training = inputs['splits'][split]
        rows = [{**r, 'canonical_smiles': structures[r['inchikey']]} for r in training['rows']]
        reports[split] = run_combined_cv(
            rows, training['descriptors'], dest / split / 'combined', dest / split / 'fingerprint_plan',
            source_checksums=training['source_checksums'], n_jobs=n_jobs)
    return reports


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('target_id', nargs='?', default='CHEMBL284')
    parser.add_argument('--split', choices=('scaffold', 'random', 'both'), default='both')
    parser.add_argument('--n-jobs', type=int, default=1)
    parser.add_argument('--cache-dir', type=Path)
    args = parser.parse_args(argv)
    target = args.target_id.strip().upper()
    reports = run_combined_search(DATA_RAW / target, DATA_CURATED / target,
        RESULTS / 'splits' / target, DATA_FEATURES / target,
        args.cache_dir or RESULTS / 'gate4_training' / target,
        splits=('scaffold', 'random') if args.split == 'both' else (args.split,), n_jobs=args.n_jobs)
    for split, report in reports.items():
        print(f"{split} combined search complete: {report['identity_sha256']}; training-only")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
