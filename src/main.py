"""
main.py
Baseline driver for CLT:
- Build laminate from materials.csv
- Compute ABD and mid-plane response with full A–B–D coupling
- Estimate SSSS plate center deflection via Navier series
- Run an angle sweep and export CSV + plot

Assumes repo layout:
  data/materials.csv
  data/sweeps/
  figures/
"""

from __future__ import annotations
import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt

from utils import load_materials, deg2rad, DATA_DIR
from clt import (
    Ply,
    laminate_abd,
    midplane_response,
    ply_strains_stresses,
    navier_center_deflection,  # from Block 8
)

# --- paths --------------------------------------------------------------------
ROOT = Path(__file__).resolve().parents[1]
FIG_DIR = ROOT / "figures"
SWEEP_DIR = ROOT / "data" / "sweeps"
FIG_DIR.mkdir(parents=True, exist_ok=True)
SWEEP_DIR.mkdir(parents=True, exist_ok=True)

# --- configuration ------------------------------------------------------------
# Plate geometry and load for Navier check
A_LEN = 0.30  # m
B_LEN = 0.30  # m
Q_UNIFORM = 1_000.0  # N/m^2

# Baseline laminate
BASE_STACK = [0, 45, -45, 90, 90, -45, 45, 0]  # degrees
# Angle-sweep laminate: all plies at the same angle (keeps thickness constant)
SWEEP_STACK_COUNT = 8  # plies; each ply uses the same sweep angle θ

# ------------------------------------------------------------------------------
def build_plies(material_row: pd.Series, angles_deg: list[float]) -> list[Ply]:
    """Construct Ply objects from a materials.csv row and an angle list."""
    E1 = float(material_row["E1"])
    E2 = float(material_row["E2"])
    G12 = float(material_row["G12"])
    v12 = float(material_row["v12"])
    tply = float(material_row["t_ply"])
    return [Ply(E1, E2, G12, v12, deg2rad(th), tply) for th in angles_deg]

def effective_Ex_from_A(A: np.ndarray, total_t: float) -> float:
    """
    Effective Ex under pure membrane loading.
    Define Nx=1 N/m, Ny=Nxy=0, κ=0 ⇒ ε0 = A^{-1} N.
    Ex_eff = (σx/εx) with σx = Nx/t.
    """
    N = np.array([1.0, 0.0, 0.0], float)            # N/m
    eps0 = np.linalg.solve(A, N)                    # κ = 0 by definition
    ex = float(eps0[0])
    return (N[0] / total_t) / ex                    # Pa


def baseline():
    # --- load material ---------------------------------------------------------
    mats = load_materials(DATA_DIR / "materials.csv")
    mat = mats.iloc[0]  # use first material as baseline

    # --- build laminate --------------------------------------------------------
    plies = build_plies(mat, BASE_STACK)
    A, B, D, z = laminate_abd(plies)
    total_t = sum(p.t for p in plies)

    # --- mid-plane response under unit bending about x (diagnostic) ------------
    eps0_bx, kappa_bx = midplane_response(A, B, D, Nx=0.0, Ny=0.0, Nxy=0.0,
                                          Mx=1.0, My=0.0, Mxy=0.0, method="schur")

    # --- Navier SSSS center deflection under uniform pressure ------------------
    w0 = navier_center_deflection(D, A_LEN, B_LEN, Q_UNIFORM, max_odd=5)

    # --- print diagnostics -----------------------------------------------------
    print("=== Baseline Laminate =========================================")
    print(f"Stack: {BASE_STACK}")
    print(f"Total thickness t = {total_t*1e3:.3f} mm")
    print(f"||A|| = {np.linalg.norm(A):.6e}")
    print(f"||B|| = {np.linalg.norm(B):.6e}")
    print(f"||D|| = {np.linalg.norm(D):.6e}")
    print(f"kappa (unit Mx): {kappa_bx}")
    print(f"SSSS Navier center deflection at q=1 kPa: {w0*1e3:.6f} mm")

    # --- save a quick per-ply table at mid-surface for reference ----------------
    ply_out = ply_strains_stresses(plies, z, eps0_bx, kappa_bx)
    df_ply = pd.DataFrame([{
        "ply_id": po["k"] if "k" in po else i,
        "theta_deg": BASE_STACK[i],
        "z_bot_mm": 1e3*float(z[i]),
        "z_top_mm": 1e3*float(z[i+1]),
        "sig_top_xx_Pa": float(po["sig_top_xy"][0]),
        "sig_bot_xx_Pa": float(po["sig_bot_xy"][0]),
    } for i, po in enumerate(ply_out)])
    out_csv = SWEEP_DIR / "baseline_ply_summary.csv"
    df_ply.to_csv(out_csv, index=False)
    print(f"Saved: {out_csv}")

def angle_sweep():
    # --- load material ---------------------------------------------------------
    mats = load_materials(DATA_DIR / "materials.csv")
    mat = mats.iloc[0]

    # Sweep θ from 0→90 deg for an N-ply laminate with all plies at θ
    thetas = np.arange(0, 91, 5)
    rows = []
    for th in thetas:
        stack = [th]*SWEEP_STACK_COUNT
        plies = build_plies(mat, stack)
        A, B, D, z = laminate_abd(plies)
        t = sum(p.t for p in plies)

        Ex_eff = effective_Ex_from_A(A, t)
        rows.append({
            "angle_deg": th,
            "t_total_m": t,
            "normA": np.linalg.norm(A),
            "normB": np.linalg.norm(B),
            "normD": np.linalg.norm(D),
            "Ex_eff_Pa": Ex_eff,
        })

    df = pd.DataFrame(rows).sort_values("angle_deg")
    csv_path = SWEEP_DIR / "angle_sweep_ex.csv"
    df.to_csv(csv_path, index=False)

    # Plot Ex vs angle
    plt.figure()
    plt.plot(df["angle_deg"], df["Ex_eff_Pa"], marker="o")
    plt.xlabel("Ply angle θ [deg]")
    plt.ylabel("Effective Ex [Pa]")
    plt.title(f"Angle sweep: effective Ex with {SWEEP_STACK_COUNT}-ply [θ]_N laminate")
    figp = FIG_DIR / "angle_sweep_ex.png"
    plt.savefig(figp, dpi=200, bbox_inches="tight")

    print(f"Saved: {csv_path}")
    print(f"Saved: {figp}")

if __name__ == "__main__":
    baseline()
    angle_sweep()
