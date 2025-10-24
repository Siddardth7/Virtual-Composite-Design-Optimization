"""
main.py
One-command baseline: build ABD, compute response for a unit load,
run an angle sweep, and export CSV + plots.
"""
from __future__ import annotations
import numpy as np, pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt

from utils import load_materials, deg2rad, DATA_DIR
from clt import Ply, laminate_abd, midplane_response, ply_strains_stresses

FIG_DIR = Path(__file__).resolve().parents[1] / "figures"
SWEEP_DIR = Path(__file__).resolve().parents[1] / "data" / "sweeps"
FIG_DIR.mkdir(parents=True, exist_ok=True)
SWEEP_DIR.mkdir(parents=True, exist_ok=True)

def baseline():
    mats = load_materials(None)
    mat = mats.iloc[0]  # first material
    angles = [0, 15, 30, 45, 60, 75, 90]
    results = []
    # Geometry and loading: unit membrane load Nx = 1e3 N/m as example
    Nx = 1e3; Ny = 0.0; Nxy = 0.0
    # [theta/-theta]s symmetric 4-ply laminate with thickness from CSV
    tply = float(mat["t_ply"])
    for ang in angles:
        plies = [
            Ply(mat["E1"], mat["E2"], mat["G12"], mat["v12"], deg2rad(ang), tply),
            Ply(mat["E1"], mat["E2"], mat["G12"], mat["v12"], -deg2rad(ang), tply),
            Ply(mat["E1"], mat["E2"], mat["G12"], mat["v12"], -deg2rad(ang), tply),
            Ply(mat["E1"], mat["E2"], mat["G12"], mat["v12"], deg2rad(ang), tply),
        ]
        A,B,D,z = laminate_abd(plies)
        eps0,kappa = midplane_response(A,B,D,Nx=Nx,Ny=Ny,Nxy=Nxy)
        # Equivalent in-plane compliance Sx = eps0_x / Nx
        Sx = eps0[0]/Nx
        Ex_eff = 1.0/Sx if Sx != 0 else np.nan
        # Record
        results.append({"angle_deg": ang, "Ex_eff_Pa": Ex_eff})
    df = pd.DataFrame(results)
    csv_path = SWEEP_DIR / "angle_sweep_baseline.csv"
    df.to_csv(csv_path, index=False)
    # Plot
    plt.figure()
    plt.plot(df["angle_deg"], df["Ex_eff_Pa"], marker="o")
    plt.xlabel("Ply angle θ [deg]")
    plt.ylabel("Effective Ex [Pa]")
    plt.title("Angle sweep: effective Ex from CLT")
    figp = FIG_DIR / "angle_sweep_ex.png"
    plt.savefig(figp, dpi=200, bbox_inches="tight")
    print(f"Saved: {csv_path}")
    print(f"Saved: {figp}")

if __name__ == "__main__":
    baseline()
