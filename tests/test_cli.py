"""Exercise entry-point option validation and root discovery from another directory."""
from pathlib import Path
import subprocess
import shutil
import tempfile

root = Path(__file__).resolve().parents[1]
with tempfile.TemporaryDirectory(prefix='mppi-cli-') as temp:
    for arguments in (['--B=0'], ['--B=2.5'], ['--B=2', '--B=3'], ['--unknown=1']):
        proc = subprocess.run(['Rscript', '--vanilla', str(root / 'scripts/run_simulations.R'), *arguments],
                              cwd=temp, capture_output=True, text=True)
        assert proc.returncode != 0 and 'Error' in proc.stderr, (arguments, proc.stdout, proc.stderr)
    proc = subprocess.run(['Rscript', '--vanilla', str(root / 'scripts/reproduce.R'), '--source=invalid'],
                          cwd=temp, capture_output=True, text=True)
    assert proc.returncode != 0 and 'source must be reference or run' in proc.stderr
    proc = subprocess.run(['Rscript', '--vanilla', str(root / 'scripts/preflight.R'), '--B=2', '--results=' + temp],
                          cwd=temp, capture_output=True, text=True)
    assert proc.returncode == 0 and 'PREFLIGHT_PASSED' in proc.stdout, proc.stderr
    alias = Path(temp) / 'checkout with spaces (test)'
    alias.mkdir()
    for name in ('R', 'scripts', 'config'):
        shutil.copytree(root / name, alias / name)
    shutil.copyfile(root / 'README.md', alias / 'README.md')
    proc = subprocess.run(['Rscript', '--vanilla', str(alias / 'scripts/preflight.R'), '--B=2',
                           '--results=' + str(Path(temp) / 'path-check')], cwd=temp,
                          capture_output=True, text=True)
    assert proc.returncode == 0 and 'PREFLIGHT_PASSED' in proc.stdout, proc.stderr
print('CLI_TESTS_PASSED')
