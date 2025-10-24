"""
clt.py
Classical Laminate Theory utilities.
References: Jones (1999), Daniel & Ishai (2006).
"""
from __future__ import annotations
import numpy as np
from dataclasses import dataclass
from typing import List, Tuple

@dataclass
class Ply:
    E1: float      # Pa
    E2: float      # Pa
    G12: float     # Pa
    v12: float     # -
    theta: float   # radians
    t: float       # m

def Q_matrix(E1: float, E2: float, G12: float, v12: float) -> np.ndarray:
    v21 = v12 * E2 / E1
    denom = 1.0 - v12 * v21
    Q11 = E1 / denom
    Q22 = E2 / denom
    Q12 = v12 * E2 / denom
    Q66 = G12
    return np.array([[Q11, Q12, 0.0],
                     [Q12, Q22, 0.0],
                     [0.0,  0.0,  Q66]], dtype=float)

def T_sigma(theta: float) -> np.ndarray:
    c = np.cos(theta); s = np.sin(theta)
    c2 = c*c; s2 = s*s; cs = c*s
    return np.array([[c2, s2, 2*cs],
                     [s2, c2, -2*cs],
                     [-cs, cs, c2 - s2]], dtype=float)

def T_epsilon(theta: float) -> np.ndarray:
    c = np.cos(theta); s = np.sin(theta)
    c2 = c*c; s2 = s*s; cs = c*s
    return np.array([[c2, s2, cs],
                     [s2, c2, -cs],
                     [-2*cs, 2*cs, c2 - s2]], dtype=float)

def Q_bar(Q: np.ndarray, theta: float) -> np.ndarray:
    T = T_sigma(theta)
    Tinv = np.linalg.inv(T)
    Teps = T_epsilon(theta)
    return Tinv @ Q @ Teps

def laminate_abd(plies: List[Ply]) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Returns A, B, D, and z-coordinates per ply interface.
    Stacking sequence defined from bottom to top in 'plies'.
    """
    total_t = sum(p.t for p in plies)
    z_bot = -total_t/2.0
    z = [z_bot]
    for p in plies:
        z.append(z[-1] + p.t)
    z = np.array(z)
    A = np.zeros((3,3)); B = np.zeros((3,3)); D = np.zeros((3,3))
    for k, p in enumerate(plies):
        Qb = Q_bar(Q_matrix(p.E1, p.E2, p.G12, p.v12), p.theta)
        zk, zk1 = z[k], z[k+1]
        A += Qb * (zk1 - zk)
        B += 0.5 * Qb * (zk1**2 - zk**2)
        D += (1/3) * Qb * (zk1**3 - zk**3)
    return A, B, D, z

def midplane_response(A: np.ndarray, B: np.ndarray, D: np.ndarray,
                      Nx: float=0, Ny: float=0, Nxy: float=0,
                      Mx: float=0, My: float=0, Mxy: float=0) -> Tuple[np.ndarray, np.ndarray]:
    """
    Solve for {eps0, kappa} from {N, M}:
    [A B; B D] {eps0; kappa} = {N; M}
    """
    K = np.block([[A, B],[B, D]])
    rhs = np.array([Nx,Ny,Nxy,Mx,My,Mxy], dtype=float)
    sol = np.linalg.solve(K, rhs)
    eps0 = sol[:3]; kappa = sol[3:]
    return eps0, kappa

def ply_strains_stresses(plies: List[Ply], z: np.ndarray, eps0: np.ndarray, kappa: np.ndarray) -> List[dict]:
    out = []
    for k, p in enumerate(plies):
        z_mid = 0.5*(z[k]+z[k+1])
        exy = eps0 + kappa * z_mid
        # Local stresses using Q_bar
        Qb = Q_bar(Q_matrix(p.E1, p.E2, p.G12, p.v12), p.theta)
        sxy = Qb @ exy
        out.append({"ply": k+1, "theta_rad": p.theta, "z_mid": z_mid,
                    "ex": exy[0], "ey": exy[1], "gxy": exy[2],
                    "sx": sxy[0], "sy": sxy[1], "txy": sxy[2]})
    return out
