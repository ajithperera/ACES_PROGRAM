#!/usr/bin/env python3
"""Run the pyaces API tests (api_cases.py), each in a fresh process and directory.

    python3 run_api_tests.py [test ...]

Work directories go under $PYACES_API_RUN_DIR, default pyaces/suite_run/api
(untracked). The extension is found as usual ($PYACES_EXT_DIR, else
pyaces/Main/build). Ends with the same "=== Summary ===" block as
pyaces_suite_runner.py, which run_tests.sh reads.
"""

import os
import shutil
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
PYACES = os.path.dirname(HERE)
GENBAS = os.path.join(os.path.dirname(PYACES), 'basis', 'GENBAS')
RUN_DIR = os.environ.get('PYACES_API_RUN_DIR', os.path.join(PYACES, 'suite_run', 'api'))
TIMEOUT = 300

sys.path.insert(0, HERE)
from api_cases import TESTS, ZMATS  # noqa: E402


def run_one(name):
    zmat, _ = TESTS[name]
    wd = os.path.join(RUN_DIR, name)
    shutil.rmtree(wd, ignore_errors=True)
    os.makedirs(wd)
    with open(os.path.join(wd, 'ZMAT'), 'w') as f:
        f.write(ZMATS[zmat])
    os.symlink(GENBAS, os.path.join(wd, 'GENBAS'))
    try:
        p = subprocess.run(['python3', '-u', os.path.join(HERE, 'api_cases.py'), name],
                           cwd=wd, capture_output=True, text=True, timeout=TIMEOUT)
    except subprocess.TimeoutExpired:
        return 'TIMEOUT', f'after {TIMEOUT} s', wd
    with open(os.path.join(wd, 'run.log'), 'w') as f:
        f.write('=== stdout ===\n' + p.stdout + '\n=== stderr ===\n' + p.stderr)
    if 'API_TEST_RESULT: PASS' in p.stdout and p.returncode == 0:
        return 'PASS', '', wd
    if 'API_TEST_RESULT: FAIL' in p.stdout:
        lines = [l for l in (p.stdout + p.stderr).splitlines() if l.startswith(('AssertionError', 'KeyError', 'ValueError', 'TypeError', 'AttributeError'))]
        return 'FAIL', lines[-1] if lines else '', wd
    how = f'signal {-p.returncode}' if p.returncode < 0 else f'exit code {p.returncode}'
    return 'CRASH', f'{how}, no result line', wd


def main():
    names = sys.argv[1:] or list(TESTS)
    unknown = [n for n in names if n not in TESTS]
    if unknown:
        sys.exit(f'unknown test(s): {unknown}; known: {list(TESTS)}')
    os.makedirs(RUN_DIR, exist_ok=True)
    counts = {}
    for i, name in enumerate(names, 1):
        t0 = time.time()
        status, detail, wd = run_one(name)
        counts[status] = counts.get(status, 0) + 1
        print(f'[{i}/{len(names)}] {name:20s} {status:8s} {detail} ({time.time() - t0:.1f}s)', flush=True)
        if status != 'PASS':
            print(f'        log: {os.path.join(wd, "run.log")}', flush=True)
    print('\n=== Summary ===')
    for status in sorted(counts):
        print(f'  {status}: {counts[status]}')
    print()
    return 0 if set(counts) == {'PASS'} else 1


if __name__ == '__main__':
    sys.exit(main())
