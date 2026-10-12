"""Local control-flow tests; fake Rscript/nproc, no SSH and no real server jobs."""
import os, pathlib, subprocess, tempfile, time
root = pathlib.Path(__file__).resolve().parents[1]
with tempfile.TemporaryDirectory(prefix='mppi-launcher-') as temp:
    temp = pathlib.Path(temp); bin_dir = temp/'bin'; bin_dir.mkdir()
    (bin_dir/'nproc').write_text('#!/usr/bin/env bash\necho 16\n')
    (bin_dir/'Rscript').write_text('''#!/usr/bin/env bash
printf '%s\\n' "$*" >> "$MPPI_TEST_CALLS"
case "$MPPI_TEST_MODE:$*" in
 preflight:*preflight.R*) exit 1;;
 partial:*--scenario=2\\ *) exit 2;;
 infra:*--scenario=2\\ *) exit 3;;
 hold:*preflight.R*) touch "$MPPI_TEST_READY"; sleep 1;;
esac
exit 0
''')
    for f in bin_dir.iterdir(): f.chmod(0o755)
    def environment(mode, calls):
        return dict(os.environ, PATH=str(bin_dir)+os.pathsep+os.environ['PATH'],
                    MPPI_TEST_MODE=mode, MPPI_TEST_CALLS=str(calls), MPPI_B='1',
                    MPPI_AVAILABLE_CORES='8', SLURM_CPUS_PER_TASK='8')
    def run(mode, directory, **changes):
        calls=temp/(directory.name+'-'+mode+'.calls'); env=environment(mode,calls); env.update(changes)
        p=subprocess.run(['bash',str(root/'scripts/run_all_server.sh'),str(root),str(directory)],env=env,capture_output=True,text=True)
        return p, calls.read_text().splitlines() if calls.exists() else []
    out=temp/'partial'
    p,calls=run('partial',out)
    assert p.returncode==2 and len(calls)==15, (p.returncode,p.stderr)
    assert all('--workers=3' in c for c in calls if 'run_simulations.R' in c)
    assert all(c.startswith('--vanilla ') for c in calls)
    assert not (out/'.launcher_lock').exists()
    before=list((out/'logs').iterdir());p,_=run('success',out)
    assert p.returncode==0 and len(list((out/'logs').iterdir()))==len(before)+1
    p,calls=run('infra',temp/'infra');assert p.returncode==3 and len(calls)==3
    p,calls=run('preflight',temp/'preflight');assert p.returncode==1 and len(calls)==1
    for changes in ({'MPPI_B':'0'},{'MPPI_B':'2.5'},{'MPPI_AVAILABLE_CORES':'32'},
                    {'MPPI_AVAILABLE_CORES':'9'},{'MPPI_AVAILABLE_CORES':'5'},
                    {'MPPI_AVAILABLE_CORES':'','SLURM_CPUS_PER_TASK':''}):
        p,calls=run('success',temp/'invalid',**changes);assert p.returncode==1 and not calls
    ready=temp/'ready';calls=temp/'hold.calls';env=environment('hold',calls);env['MPPI_TEST_READY']=str(ready)
    out=temp/'concurrent'
    first=subprocess.Popen(['bash',str(root/'scripts/run_all_server.sh'),str(root),str(out)],env=env,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
    try:
        deadline=time.monotonic()+5
        while not ready.exists() and time.monotonic()<deadline: time.sleep(.02)
        assert ready.exists()
        second,_=run('success',out);assert second.returncode==1
        first.communicate(timeout=10);assert first.returncode==0
    finally:
        if first.poll() is None: first.terminate();first.communicate(timeout=5)
print('LAUNCHER_TESTS_PASSED')
