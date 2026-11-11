"""Run the numerically tested Ti–Ce research workflow from any working directory.

Outputs are overwritten reproducibly in figures/, outputs/, validation_figures/
validation_outputs/, extended_outputs/ and extended_figures/. The model is illustrative, not a validated digital twin.
"""
import argparse
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parent

def main() -> None:
    parser = argparse.ArgumentParser(description='Reproduce the Ti–Ce multiphysics numerical study')
    parser.add_argument('--full', action='store_true', help='Use full-resolution integrated solver (slower)')
    args = parser.parse_args()
    env = os.environ.copy()
    env.setdefault('OPENBLAS_NUM_THREADS','1')
    for args2 in (["run_integrated.py"] if args.full else ["run_integrated.py", "--fast"],
                  ["calibrate_transport.py"], ["design_tradeoffs.py"], ["run_extensions.py"], ["run_evidence_thermodynamics.py"]):
        print('Running:', ' '.join(args2), flush=True)
        subprocess.run([sys.executable, *args2], cwd=ROOT, env=env, check=True)
    print('All five numerical workflows completed. See figures/, outputs/, validation_figures/, validation_outputs/, extended_figures/, extended_outputs/, thermo_evidence/.')

if __name__ == '__main__':
    main()
