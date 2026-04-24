"""
run_pipeline.py
Run the full VirtualCompositeDesign pipeline in order:
  1. Baseline CLT analysis + angle sweep  (main.py)
  2. CLT vs FEA comparison               (compare_clt_fea.py)
  3. SA layup optimizer                  (layup_optimizer_sa.py)
  4. SA spot-check (top-3 candidates)    (sa_spotcheck.py)

Usage:
    python src/run_pipeline.py
"""
from __future__ import annotations
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PYTHON = sys.executable

STEPS = [
    ("Baseline CLT + angle sweep",   ROOT / "src" / "main.py"),
    ("CLT vs FEA comparison",        ROOT / "src" / "compare_clt_fea.py"),
    ("SA layup optimizer",           ROOT / "src" / "layup_optimizer_sa.py"),
    ("SA spot-check",                ROOT / "src" / "sa_spotcheck.py"),
]


def run_step(name: str, script: Path) -> bool:
    print(f"\n{'='*60}")
    print(f"  STEP: {name}")
    print(f"{'='*60}\n")
    result = subprocess.run([PYTHON, str(script)], cwd=str(ROOT))
    if result.returncode != 0:
        print(f"\n  STEP FAILED: {name} (exit code {result.returncode})")
        return False
    print(f"\n  STEP DONE: {name}")
    return True


def main() -> None:
    print("\n=== VirtualCompositeDesign — Full Pipeline ===")
    for name, script in STEPS:
        if not run_step(name, script):
            sys.exit(1)
    print("\n=== All steps complete ===\n")


if __name__ == "__main__":
    main()
