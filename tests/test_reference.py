"""Check immutable reference hashes, method rosters and true-null rejection rates."""
import csv
import hashlib
import json
from pathlib import Path

root = Path(__file__).resolve().parents[1]
ref = root / 'reference'
provenance = json.loads((ref / 'provenance.json').read_text())
for name, digest in provenance['reference_sha256'].items():
    assert hashlib.sha256((ref / name).read_bytes()).hexdigest() == digest, name
for name, digest in provenance.get('figures_sha256', {}).items():
    assert hashlib.sha256((ref / 'figures' / name).read_bytes()).hexdigest() == digest, name

def read(name):
    with (ref / name).open() as stream:
        return list(csv.DictReader(stream))

summary = read('simulation_summary.csv')
assert len(summary) == 108
keys = {(int(row['scenario_id']), row['method']) for row in summary}
assert len(keys) == 108 and {id for id, _ in keys} == set(range(1, 13))
for row in summary:
    assert int(row['B_expected']) == int(row['B_success']) == int(row['B_valid']) == 1000
    assert all(int(row[key]) == 0 for key in ('B_failed', 'B_missing', 'B_fallback'))
    assert 0 <= float(row['coverage']) <= 1 and float(row['mean_volume']) > 0
    if row['method'] == 'MPPI':
        weights = [float(row[f'mean_w{i}']) for i in range(4)]
        assert min(weights) >= -1e-10 and abs(sum(weights) - 1) < 1e-8
lookup = {(int(row['scenario_id']), row['method']): float(row['coverage']) for row in summary}
for row in read('mean_test_curves.csv'):
    if abs(float(row['distance'])) < 1e-10:
        assert abs(float(row['rejection_frequency']) - (1 - lookup[int(row['scenario_id']), row['method']])) < 1e-10
for row in read('regression_test_curves.csv'):
    if float(row['radius']) == 0:
        assert abs(float(row['rejection_frequency']) - (1 - lookup[int(row['scenario_id']), row['method']])) < 1e-10
# Freeze the formally reported domain-shift rows to catch an accidental old-reference mix.
volumes = {(int(row['scenario_id']), row['method']): float(row['mean_volume']) for row in summary}
for scenario, expected in [(3, 3.7460), (6, 0.9953), (9, 1.7087), (12, 31.7663)]:
    assert abs(volumes[scenario, 'MPPI'] - expected) < 0.0001, (scenario, volumes[scenario, 'MPPI'])
normal = read('normal_diagnostics.csv')
assert len(normal) == 6000
for id in (1, 2, 3, 7, 8, 9):
    assert sum(int(row['scenario_id']) == id for row in normal) == 1000
assert all(row['method'] == 'MPPI' for row in normal)
assert len(read('dxa_results.csv')) == 54
print('REFERENCE_TESTS_PASSED')
