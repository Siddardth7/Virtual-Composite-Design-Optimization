"""
compare_clt_fea.py  — Checkpoint 1 (D06)
CLT-vs-FEA comparison pipeline for the [0/45/-45/90]s baseline laminate.

Steps
-----
1. Compute CLT Navier centre deflection and peak σₓₓ (bottom of 0° ply).
2. Load FEA results:
     a. If ccx binary is available → generate .inp, run CalculiX, parse .dat.
     b. Otherwise → load from fea/results/fea_summary.csv (pre-computed run).
3. Compute % errors for both metrics.
4. Print formatted comparison table.
5. Save table to data/fea_comparison.csv.
6. Generate side-by-side bar chart with % error annotations → figures/fea_comparison.png.
7. CHECKPOINT: report PASS / FAIL for both errors < 5%.

CLT reference
-------------
  Layup    : [0, 45, -45, 90]_s   (IM7/8552, t_total = 1 mm)
  Plate    : 300 × 300 mm SSSS
  Load     : q = 1 000 Pa uniform pressure
  Baseline : w_CLT = 0.064411 mm   (Navier series, validated in main.py)
"""

from __future__ import annotations
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ── repo layout ───────────────────────────────────────────────────────────────
ROOT     = Path(__file__).resolve().parents[1]
FEA_DIR  = ROOT / "fea"
DATA_DIR = ROOT / "data"
FIG_DIR  = ROOT / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

# ── add src/ to path for local imports ───────────────────────────────────────
sys.path.insert(0, str(Path(__file__).resolve().parent))
from utils import load_materials, deg2rad
from clt import (
    Ply,
    laminate_abd,
    navier_center_deflection,
    Q_matrix,
    Q_bar,
)

# ── plate / laminate constants (match main.py / generate_inp.py) ──────────────
LX        = 0.300          # m
LY        = 0.300          # m
Q_LOAD    = 1_000.0        # Pa
BASE_STACK = [0, 45, -45, 90, 90, -45, 45, 0]
NAVIER_ODD = 5             # Navier series terms (same as main.py)

# ── file paths ────────────────────────────────────────────────────────────────
INP_FILE    = FEA_DIR / "abaqus_inputs" / "plate_ss.inp"
DAT_FILE    = FEA_DIR / "abaqus_inputs" / "plate_ss.dat"
FEA_CSV     = FEA_DIR / "results" / "fea_summary.csv"
OUT_CSV     = DATA_DIR / "fea_comparison.csv"
OUT_PNG     = FIG_DIR  / "fea_comparison.png"

# ── error threshold for CHECKPOINT ───────────────────────────────────────────
ERROR_THRESHOLD_PCT = 5.0


# ════════════════════════════════════════════════════════════════════════════
# 1. CLT computation
# ════════════════════════════════════════════════════════════════════════════

def compute_clt_values() -> dict:
    """
    Return CLT centre deflection (m) and peak σₓₓ (Pa) for the baseline laminate.

    σₓₓ is taken as |stress| at the surface of the critical ply (0° bottom ply,
    z = −t/2).  The physical bending causes tension there under downward pressure;
    the CLT Navier formula uses a downward-positive convention so the raw value
    is negative — we return the magnitude.
    """
    mats = load_materials(DATA_DIR / "materials.csv")
    mat  = mats.iloc[0]
    E1   = float(mat["E1"])
    E2   = float(mat["E2"])
    G12  = float(mat["G12"])
    v12  = float(mat["v12"])
    tply = float(mat["t_ply"])

    plies  = [Ply(E1, E2, G12, v12, deg2rad(th), tply) for th in BASE_STACK]
    A, B, D, z = laminate_abd(plies)

    # ── centre deflection ─────────────────────────────────────────────────
    w_m = navier_center_deflection(D, LX, LY, Q_LOAD, max_odd=NAVIER_ODD)

    # ── curvatures at plate centre via Navier series ──────────────────────
    D11, D22, D12, D66 = D[0, 0], D[1, 1], D[0, 1], D[2, 2]
    kx = 0.0
    ky = 0.0
    for m in range(1, NAVIER_ODD + 1, 2):
        for n in range(1, NAVIER_ODD + 1, 2):
            mpa  = m * np.pi / LX
            npb  = n * np.pi / LY
            dnom = (D11 * mpa**4
                    + 2.0 * (D12 + 2.0 * D66) * mpa**2 * npb**2
                    + D22 * npb**4)
            Wmn  = (16.0 * Q_LOAD) / (np.pi**6 * m**2 * n**2 * dnom)
            sm   = np.sin(m * np.pi / 2)
            sn   = np.sin(n * np.pi / 2)
            kx  += Wmn * mpa**2 * sm * sn
            ky  += Wmn * npb**2 * sm * sn

    # ── σₓₓ at bottom of 0° ply (z = −t/2) ──────────────────────────────
    z_bot_ply0 = z[0]                             # most negative z interface
    Q0   = Q_matrix(E1, E2, G12, v12)
    Qb0  = Q_bar(Q0, 0.0)                         # 0° ply → no rotation
    kappa = np.array([kx, ky, 0.0])
    sig_bot = Qb0 @ (z_bot_ply0 * kappa)          # [σxx, σyy, τxy] in Pa
    sigma_xx_pa = abs(float(sig_bot[0]))           # magnitude (physical = tension)

    return {
        "deflection_m":    float(w_m),
        "sigma_xx_max_Pa": sigma_xx_pa,
    }


# ════════════════════════════════════════════════════════════════════════════
# 2. FEA results loader
# ════════════════════════════════════════════════════════════════════════════

def _find_ccx() -> str | None:
    """Return path to ccx binary if available on PATH or as fea/ccx symlink."""
    # check symlink in repo first
    repo_ccx = FEA_DIR / "ccx"
    if repo_ccx.exists() and not repo_ccx.is_symlink():
        return str(repo_ccx)
    if repo_ccx.is_symlink() and repo_ccx.resolve().exists():
        return str(repo_ccx)
    # fall back to PATH
    import shutil
    return shutil.which("ccx")


def _run_fea(ccx_bin: str) -> dict:
    """
    Generate .inp, run CalculiX, parse .dat, return FEA result dict.
    Must be called from fea/abaqus_inputs/ so ccx can find the .inp file.
    """
    from fea.generate_inp import generate_full_inp   # type: ignore[import]
    from fea.parse_ccx_results import parse_dat       # type: ignore[import]

    print("  Generating CalculiX input file …")
    generate_full_inp(out_path=INP_FILE)

    run_dir  = INP_FILE.parent
    job_name = INP_FILE.stem   # "plate_ss"
    cmd = [ccx_bin, job_name]
    print(f"  Running: {' '.join(cmd)}  (cwd={run_dir})")
    result = subprocess.run(cmd, cwd=str(run_dir), capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"CalculiX failed (rc={result.returncode}):\n{result.stderr[:800]}"
        )

    dat_path = run_dir / f"{job_name}.dat"
    print(f"  Parsing {dat_path.name} …")
    parsed = parse_dat(dat_path)

    return {
        "deflection_m":    abs(parsed["deflection_m"]),   # U3 is negative downward
        "sigma_xx_max_Pa": abs(parsed["sigma_xx_max_Pa"]),
        "notes":           parsed.get("notes", "ccx run"),
    }


def load_fea_values() -> dict:
    """
    Try to obtain FEA values by:
      1. Running CalculiX if ccx is on PATH.
      2. Loading from fea/results/fea_summary.csv (pre-computed).
    """
    ccx = _find_ccx()
    if ccx:
        print(f"  ccx found at: {ccx}")
        try:
            return _run_fea(ccx)
        except Exception as exc:
            print(f"  WARNING: ccx run failed — {exc}")
            print("  Falling back to stored FEA results.")

    if not FEA_CSV.exists():
        raise FileNotFoundError(
            f"No pre-computed FEA results found at {FEA_CSV}\n"
            "Run CalculiX manually: cd fea/abaqus_inputs && ccx plate_ss"
        )

    df  = pd.read_csv(FEA_CSV)
    row = df[df["case_id"] == "baseline"].iloc[0]
    print(f"  Loaded stored FEA results from {FEA_CSV.relative_to(ROOT)}")
    return {
        "deflection_m":    float(row["deflection_m"]),
        "sigma_xx_max_Pa": float(row["sigma_xx_max_Pa"]),
        "notes":           str(row.get("notes", "stored")),
    }


# ════════════════════════════════════════════════════════════════════════════
# 3. Error computation + table
# ════════════════════════════════════════════════════════════════════════════

def compute_errors(clt: dict, fea: dict) -> list[dict]:
    """Return per-metric comparison rows."""
    rows = []
    metrics = [
        ("Deflection",  "mm",  1e3,  "deflection_m"),
        ("σₓₓ (peak)",  "MPa", 1e-6, "sigma_xx_max_Pa"),
    ]
    for label, unit, scale, key in metrics:
        clt_val = clt[key] * scale
        fea_val = fea[key] * scale
        err_pct = abs(clt_val - fea_val) / fea_val * 100.0 if fea_val != 0.0 else float("nan")
        rows.append({
            "Metric":   label,
            "Unit":     unit,
            "CLT":      round(clt_val, 6),
            "FEA":      round(fea_val, 6),
            "Err_%":    round(err_pct, 3),
            "PASS":     err_pct < ERROR_THRESHOLD_PCT,
        })
    return rows


def print_table(rows: list[dict]) -> None:
    hdr = f"{'Metric':<18} {'Unit':<5} {'CLT':>12} {'FEA':>12} {'Err %':>8}  Status"
    sep = "─" * len(hdr)
    print(f"\n{sep}")
    print(hdr)
    print(sep)
    for r in rows:
        status = "✅ PASS" if r["PASS"] else "❌ FAIL"
        print(f"{r['Metric']:<18} {r['Unit']:<5} {r['CLT']:>12.6f} {r['FEA']:>12.6f} "
              f"{r['Err_%']:>7.3f}%  {status}")
    print(sep)


# ════════════════════════════════════════════════════════════════════════════
# 4. Save CSV
# ════════════════════════════════════════════════════════════════════════════

def save_csv(rows: list[dict], clt: dict, fea: dict) -> None:
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(rows)[["Metric", "Unit", "CLT", "FEA", "Err_%", "PASS"]]
    df.to_csv(OUT_CSV, index=False)
    print(f"\n  Saved: {OUT_CSV.relative_to(ROOT)}")


# ════════════════════════════════════════════════════════════════════════════
# 5. Bar chart
# ════════════════════════════════════════════════════════════════════════════

def save_chart(rows: list[dict]) -> None:
    labels   = [r["Metric"] for r in rows]
    clt_vals = [r["CLT"]    for r in rows]
    fea_vals = [r["FEA"]    for r in rows]
    err_pcts = [r["Err_%"]  for r in rows]
    units    = [r["Unit"]   for r in rows]

    x     = np.arange(len(labels))
    width = 0.35

    fig, ax = plt.subplots(figsize=(8, 5))
    bars_clt = ax.bar(x - width / 2, clt_vals, width, label="CLT (Navier)",
                      color="#2C7BB6", alpha=0.85, edgecolor="white")
    bars_fea = ax.bar(x + width / 2, fea_vals, width, label="FEA (CalculiX S4)",
                      color="#D7191C", alpha=0.85, edgecolor="white")

    # % error annotations centred above each pair
    for i, (bc, bf, err, unit) in enumerate(zip(bars_clt, bars_fea, err_pcts, units)):
        y_top = max(bc.get_height(), bf.get_height())
        ax.annotate(
            f"{err:.2f}%",
            xy=(x[i], y_top),
            xytext=(0, 8),
            textcoords="offset points",
            ha="center", va="bottom",
            fontsize=10, color="#555555",
        )

    ax.set_xticks(x)
    ax.set_xticklabels([f"{l}\n[{u}]" for l, u in zip(labels, units)], fontsize=11)
    ax.set_ylabel("Value", fontsize=11)
    ax.set_title(
        f"CLT vs FEA — [0/45/−45/90]s IM7/8552, 300×300 mm, q = 1 kPa",
        fontsize=12, pad=14,
    )
    ax.legend(fontsize=10)
    ax.yaxis.grid(True, linestyle="--", alpha=0.4)
    ax.set_axisbelow(True)

    # threshold line annotation
    ax.annotate(
        f"Error threshold: {ERROR_THRESHOLD_PCT:.0f}%",
        xy=(1, 1), xycoords="axes fraction",
        xytext=(-8, -10), textcoords="offset points",
        ha="right", va="top", fontsize=8, color="grey",
    )

    fig.tight_layout()
    fig.savefig(OUT_PNG, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {OUT_PNG.relative_to(ROOT)}")


# ════════════════════════════════════════════════════════════════════════════
# 6. CHECKPOINT review
# ════════════════════════════════════════════════════════════════════════════

def checkpoint_review(rows: list[dict]) -> bool:
    all_pass = all(r["PASS"] for r in rows)
    print(f"\n{'═'*50}")
    print(f"  CHECKPOINT 1 REVIEW — D06 — {__import__('datetime').date.today()}")
    print(f"{'═'*50}")
    for r in rows:
        status = "✅ PASS" if r["PASS"] else "❌ FAIL"
        print(f"  {r['Metric']:<18} error = {r['Err_%']:.3f}%  {status}")
    print(f"{'═'*50}")
    if all_pass:
        print("  ✅ CHECKPOINT 1 PASSED — both errors < 5%")
        print("  Phase 1+2 pipeline validated. Ready for Hashin implementation.")
    else:
        print("  ❌ CHECKPOINT 1 FAILED — investigate FEA mesh / BCs / material.")
    print(f"{'═'*50}\n")
    return all_pass


# ════════════════════════════════════════════════════════════════════════════
# main
# ════════════════════════════════════════════════════════════════════════════

def main() -> None:
    print("\n=== compare_clt_fea.py — Checkpoint 1 ===\n")

    # 1. CLT
    print("── Step 1: CLT Navier analysis ─────────────────────────")
    clt = compute_clt_values()
    print(f"  Deflection : {clt['deflection_m']*1e3:.6f} mm")
    print(f"  σₓₓ peak   : {clt['sigma_xx_max_Pa']/1e6:.4f} MPa")

    # 2. FEA
    print("\n── Step 2: FEA results ─────────────────────────────────")
    fea = load_fea_values()
    print(f"  Deflection : {fea['deflection_m']*1e3:.6f} mm")
    print(f"  σₓₓ peak   : {fea['sigma_xx_max_Pa']/1e6:.4f} MPa")
    print(f"  Notes      : {fea.get('notes','')}")

    # 3. Errors + table
    print("\n── Step 3: Comparison table ────────────────────────────")
    rows = compute_errors(clt, fea)
    print_table(rows)

    # 4. Save CSV
    print("\n── Step 4: Save outputs ────────────────────────────────")
    save_csv(rows, clt, fea)

    # 5. Save chart
    save_chart(rows)

    # 6. Checkpoint
    checkpoint_review(rows)


if __name__ == "__main__":
    main()
