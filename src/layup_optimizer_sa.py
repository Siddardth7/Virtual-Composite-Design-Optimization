"""
layup_optimizer_sa.py
Simulated Annealing optimizer for composite laminate stacking sequences.
Minimises: Navier centre deflection (via D11) + areal weight.
Material and ply geometry are loaded from data/materials.csv.
"""
from __future__ import annotations
import copy
import random
from pathlib import Path

import numpy as np

# ── project imports ───────────────────────────────────────────────────────────
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from clt import Q_matrix, Q_bar, ply_interfaces, abd_matrices, deg2rad
from utils import load_materials

# ── load material from CSV ────────────────────────────────────────────────────
_MAT_PATH = Path(__file__).resolve().parents[1] / "data" / "materials.csv"
_mat = load_materials(_MAT_PATH).iloc[0]

E1   = float(_mat["E1"])
E2   = float(_mat["E2"])
G12  = float(_mat["G12"])
nu12 = float(_mat["v12"])
t    = float(_mat["t_ply"])           # m per ply (0.000125 for IM7/8552)
RHO  = float(_mat["density"])         # kg/m³

# ── plate geometry (must match compare_clt_fea.py / sa_spotcheck.py) ─────────
LX = 0.300   # m
LY = 0.300   # m

# per-ply areal mass for this plate and material
_PLY_WEIGHT = RHO * t * LX * LY      # kg per ply

# ── allowed orientations and laminate size bounds ────────────────────────────
allowed_angles = [0, 45, -45, 90]
minPercent     = 0.1
min_full_plies = 4
max_full_plies = 40
min_half = min_full_plies // 2
max_half = max_full_plies // 2


# ── CLT helpers (thin wrappers over clt.py) ───────────────────────────────────

def compute_D_matrix(full_seq: list[int]) -> np.ndarray:
    """Return the 3×3 bending stiffness matrix D for *full_seq* using clt.py."""
    Q     = Q_matrix(E1, E2, G12, nu12)
    qbars = [Q_bar(Q, deg2rad(th)) for th in full_seq]
    z     = ply_interfaces(len(full_seq), len(full_seq) * t)
    _, _, D = abd_matrices(qbars, z)
    return D


def laminate_deflection_metric(full_seq: list[int]) -> float:
    """1/D11 — lower is stiffer (better). Penalty for non-physical D11."""
    D = compute_D_matrix(full_seq)
    D11 = D[0, 0]
    if D11 <= 0:
        return 1e12
    return 1.0 / D11


def laminate_weight(full_seq: list[int]) -> float:
    """Areal weight of laminate in kg (density × t_ply × plate area × n_plies)."""
    return len(full_seq) * _PLY_WEIGHT


def laminate_objective(full_seq: list[int]) -> float:
    return laminate_deflection_metric(full_seq) + laminate_weight(full_seq) + penalty(full_seq)


def penalty(full_seq: list[int]) -> float:
    pen = 0
    N = len(full_seq)
    half = N // 2
    if not np.allclose(full_seq[:half], full_seq[half:][::-1]):
        pen += 1e6
    count_p45 = np.sum(np.array(full_seq) == 45)
    count_m45 = np.sum(np.array(full_seq) == -45)
    pen += 1e6 * abs(count_p45 - count_m45)
    required = int(N * minPercent)
    count_0 = np.sum(np.array(full_seq) == 0)
    count_90 = np.sum(np.array(full_seq) == 90)
    count_45_total = count_p45 + count_m45
    if count_0 < required:
        pen += 1e6 * (required - count_0)
    if count_90 < required:
        pen += 1e6 * (required - count_90)
    if count_45_total < required:
        pen += 1e6 * (required - count_45_total)
    return pen


def objective(half_seq: list[int]) -> float:
    full_seq = half_seq + half_seq[::-1]
    return laminate_objective(full_seq)


def simulated_annealing(
    n_iterations: int = 10000,
    initial_temp: float = 1.0,
    cooling_rate: float = 0.999,
) -> tuple[list[int], float]:
    current_length = random.randint(min_half, max_half)
    current_half = [random.choice(allowed_angles) for _ in range(current_length)]
    current_obj = objective(current_half)
    best_half = copy.deepcopy(current_half)
    best_obj = current_obj
    T = initial_temp

    for it in range(n_iterations):
        moves = ['change']
        if len(current_half) < max_half:
            moves.append('insert')
        if len(current_half) > min_half:
            moves.append('delete')
        move = random.choice(moves)
        candidate_half = current_half.copy()
        if move == 'change':
            idx = random.randint(0, len(candidate_half) - 1)
            candidate_half[idx] = random.choice([a for a in allowed_angles if a != candidate_half[idx]])
        elif move == 'insert':
            idx = random.randint(0, len(candidate_half))
            candidate_half.insert(idx, random.choice(allowed_angles))
        elif move == 'delete':
            idx = random.randint(0, len(candidate_half) - 1)
            candidate_half.pop(idx)
        candidate_obj = objective(candidate_half)
        delta = candidate_obj - current_obj
        if delta < 0 or random.random() < np.exp(-delta / T):
            current_half = candidate_half
            current_obj = candidate_obj
            if current_obj < best_obj:
                best_half = copy.deepcopy(current_half)
                best_obj = current_obj
        T *= cooling_rate
        if (it + 1) % 1000 == 0:
            print(f"Iteration {it+1}: Best Obj = {best_obj:.3e}, Half sequence length = {len(best_half)}")

    best_full_seq = best_half + best_half[::-1]
    return best_full_seq, best_obj


if __name__ == "__main__":
    best_seq, best_value = simulated_annealing(n_iterations=10000)
    print("\nOptimized Layup Sequence (degrees):")
    print(best_seq)
    print("\nNumber of plies (full laminate):", len(best_seq))
    print("\nObjective Value (deflection metric + weight + penalties):", best_value)
