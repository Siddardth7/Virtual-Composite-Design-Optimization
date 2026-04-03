"""
clt.py
Classical Laminate Theory utilities.
Inputs: material constants, ply angles, ply thicknesses, applied mid-plane loads (Nx, Ny, Nxy, Mx, My, Mxy)
Outputs: [A],[B],[D] matrices, mid-plane strains/curvatures, ply stresses, Tsai-Wu/Hashin indices.
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

# --- Block 1: angles and Q ---------------------------------
def deg2rad(theta_deg: float) -> float:
    return np.deg2rad(theta_deg)

def nu21(nu12: float, E1: float, E2: float) -> float:
    return nu12 * E2 / E1

def Q_matrix(E1: float, E2: float, G12: float, nu12: float) -> np.ndarray:
    """Plane-stress reduced stiffness in 1–2 axes. Units: Pa."""
    v21 = nu21(nu12, E1, E2)
    Delta = 1.0 - nu12 * v21
    Q11 = E1 / Delta
    Q22 = E2 / Delta
    Q12 = nu12 * E2 / Delta
    Q66 = G12
    return np.array([
        [Q11, Q12, 0.0],
        [Q12, Q22, 0.0],
        [0.0,  0.0,  Q66]
    ], dtype=float)


# --- Block 2: Qbar closed-form ------------------------------
def Q_bar(Q: np.ndarray, theta_rad: float) -> np.ndarray:
    """Transformed reduced stiffness in laminate x–y axes for angle theta."""
    m = np.cos(theta_rad)
    s = np.sin(theta_rad)
    m2, s2 = m*m, s*s
    m4, s4 = m2*m2, s2*s2
    ms2 = m2*s2

    Q11, Q22 = Q[0,0], Q[1,1]
    Q12, Q66 = Q[0,1], Q[2,2]

    Qb11 = Q11*m4 + 2*(Q12+2*Q66)*ms2 + Q22*s4
    Qb22 = Q11*s4 + 2*(Q12+2*Q66)*ms2 + Q22*m4
    Qb12 = (Q11+Q22-4*Q66)*ms2 + Q12*(m4 + s4)
    Qb16 = (Q11 - Q12 - 2*Q66)*m*m*m*s - (Q22 - Q12 - 2*Q66)*m*s*s*s
    Qb26 = (Q11 - Q12 - 2*Q66)*m*s*s*s - (Q22 - Q12 - 2*Q66)*m*m*m*s
    Qb66 = (Q11 + Q22 - 2*Q12 - 2*Q66)*ms2 + Q66*(m4 + s4)

    return np.array([
        [Qb11, Qb12, Qb16],
        [Qb12, Qb22, Qb26],
        [Qb16, Qb26, Qb66]
    ], dtype=float)

# --- Block 3: z interfaces and ABD --------------------------
def ply_interfaces(nplies: int, total_t: float) -> np.ndarray:
    """z from -t/2 to +t/2 with length nplies+1."""
    return np.linspace(-total_t/2.0, total_t/2.0, nplies+1, dtype=float)

def abd_matrices(qbar_list: list[np.ndarray], z_iface: np.ndarray) -> tuple[np.ndarray,np.ndarray,np.ndarray]:
    """Assemble A,B,D from per-ply Qbar and interface z."""
    A = np.zeros((3,3)); B = np.zeros((3,3)); D = np.zeros((3,3))
    for k in range(len(qbar_list)):
        zk, zk1 = z_iface[k], z_iface[k+1]
        dz = zk1 - zk
        A += qbar_list[k] * dz
        B += 0.5 * qbar_list[k] * (zk1**2 - zk**2)
        D += (1.0/3.0) * qbar_list[k] * (zk1**3 - zk**3)
    return A, B, D

# --- Block 4: coupled mid-plane solve -----------------------
def solve_midplane(A: np.ndarray, B: np.ndarray, D: np.ndarray,
                   N: np.ndarray, M: np.ndarray,
                   method: str = "schur") -> tuple[np.ndarray,np.ndarray]:
    """
    Solve [A B; B D]{eps0;kappa}={N;M}.
    N,M in [Nx, Ny, Nxy], [Mx, My, Mxy].
    """
    if method == "direct":
        K = np.block([[A,B],[B,D]])
        rhs = np.hstack([N, M])
        sol = np.linalg.solve(K, rhs)
        return sol[:3], sol[3:]
    # Schur complement
    AinvN = np.linalg.solve(A, N)
    AinvB = np.linalg.solve(A, B)
    S = D - B @ AinvB
    kappa = np.linalg.solve(S, M - B @ AinvN)
    eps0  = AinvN - AinvB @ kappa
    return eps0, kappa

# --- Block 5: ply-level recovery (global and local) ----------
def T_sigma(theta_rad: float) -> np.ndarray:
    """Stress transform (global->local) for plane stress."""
    m, s = np.cos(theta_rad), np.sin(theta_rad)
    return np.array([
        [ m*m,      s*s,      2*m*s],
        [ s*s,      m*m,     -2*m*s],
        [-m*s,       m*s,  m*m - s*s]
    ], dtype=float)

def T_epsilon(theta_rad: float) -> np.ndarray:
    """Strain transform (global->local) for engineering strains."""
    m, s = np.cos(theta_rad), np.sin(theta_rad)
    return np.array([
        [ m*m,      s*s,      m*s],
        [ s*s,      m*m,     -m*s],
        [-2*m*s,  2*m*s,  m*m - s*s]
    ], dtype=float)

def recover_ply_global(eps0: np.ndarray, kappa: np.ndarray,
                       z_top: float, z_bot: float,
                       qbar: np.ndarray) -> dict:
    """
    Return strains/stresses at top and bottom in GLOBAL axes (x,y,xy).
    """
    eps_top = eps0 + z_top * kappa
    eps_bot = eps0 + z_bot * kappa
    sig_top = qbar @ eps_top
    sig_bot = qbar @ eps_bot
    return {
        "eps_top_xy": eps_top, "eps_bot_xy": eps_bot,
        "sig_top_xy": sig_top, "sig_bot_xy": sig_bot
    }

def recover_ply_local(eps_xy: np.ndarray, sig_xy: np.ndarray, theta_rad: float) -> tuple[np.ndarray,np.ndarray]:
    """
    Transform a single point's (eps_xy, sig_xy) to LOCAL (1,2,12).
    """
    Te = T_epsilon(theta_rad)
    Ts = T_sigma(theta_rad)
    eps_12 = Te @ eps_xy
    sig_12 = Ts @ sig_xy
    return eps_12, sig_12

# --- Block 6: Hashin (1980) failure indices ------------------
def hashin(stress_local: np.ndarray,
           X_T: float, X_C: float,
           Y_T: float, Y_C: float,
           S12: float) -> dict:
    """
    Hashin 1980 failure indices in local (1-2) axes.
    stress_local = [sigma1, sigma2, tau12] (Pa). Strengths in Pa.
    Returns dict with FI_FT, FI_FC, FI_MT, FI_MC (each float; 0 when mode inactive).

    Fibre Tension  (FT, σ₁ ≥ 0): FI = (σ₁/X_T)² + (τ₁₂/S₁₂)²
    Fibre Compr.   (FC, σ₁ < 0): FI = (σ₁/X_C)²
    Matrix Tension (MT, σ₂ ≥ 0): FI = (σ₂/Y_T)² + (τ₁₂/S₁₂)²
    Matrix Compr.  (MC, σ₂ < 0): Hashin-Rotem simplified form
    """
    s1  = float(stress_local[0])
    s2  = float(stress_local[1])
    t12 = float(stress_local[2])

    FI_FT = (s1 / X_T)**2 + (t12 / S12)**2 if s1 >= 0.0 else 0.0
    FI_FC = (s1 / X_C)**2                   if s1 <  0.0 else 0.0
    FI_MT = (s2 / Y_T)**2 + (t12 / S12)**2  if s2 >= 0.0 else 0.0
    if s2 < 0.0:
        FI_MC = ((s2 / (2.0 * S12))**2
                 + ((Y_C / (2.0 * S12))**2 - 1.0) * (s2 / Y_C)
                 + (t12 / S12)**2)
    else:
        FI_MC = 0.0

    return {"FI_FT": FI_FT, "FI_FC": FI_FC, "FI_MT": FI_MT, "FI_MC": FI_MC}


# --- Block 7: Tsai–Wu failure index -------------------------
def tsai_wu(sig12: np.ndarray,
            Xt: float, Xc: float,
            Yt: float, Yc: float,
            S12: float) -> float:
    """
    Tsai-Wu failure index in local axes.
    sig12 = [sigma1, sigma2, tau12] (Pa)
    Strengths Xt,Xc,Yt,Yc,S12 in Pa. Return scalar FI.
    """
    s1, s2, t12 = sig12
    F1  = 1.0/Xt - 1.0/Xc
    F2  = 1.0/Yt - 1.0/Yc
    F11 = 1.0/(Xt*Xc)
    F22 = 1.0/(Yt*Yc)
    F66 = 1.0/(S12*S12)
    F12 = -0.5 * np.sqrt(F11 * F22)
    FI = F1*s1 + F2*s2 + F11*s1*s1 + F22*s2*s2 + F66*t12*t12 + 2.0*F12*s1*s2
    return float(FI)

# --- Block 7: end-to-end laminate evaluator ----------------
def evaluate_laminate(E1: float, E2: float, G12: float, nu12: float,
                      angles_deg: list[float], ply_t: float,
                      N: np.ndarray, M: np.ndarray,
                      strengths: dict | None = None):
    """
    Returns dict with A,B,D, eps0, kappa, and per-ply local/global results.

    Optional *strengths* dict enables Tsai-Wu and Hashin indices per ply:
        {"X_T", "X_C", "Y_T", "Y_C", "S12"}  (all in Pa)
    When provided, each ply entry gains:
        tsai_wu_top, tsai_wu_bot,
        hashin_FT_top, hashin_FC_top, hashin_MT_top, hashin_MC_top,
        hashin_FT_bot, hashin_FC_bot, hashin_MT_bot, hashin_MC_bot
    """
    Q = Q_matrix(E1,E2,G12,nu12)
    qbars = [Q_bar(Q, deg2rad(th)) for th in angles_deg]
    total_t = ply_t * len(angles_deg)
    z_iface = ply_interfaces(len(angles_deg), total_t)
    A,B,D = abd_matrices(qbars, z_iface)
    eps0, kappa = solve_midplane(A,B,D,N,M,method="schur")

    ply_results = []
    for k in range(len(angles_deg)):
        z_bot, z_top = z_iface[k], z_iface[k+1]
        theta = deg2rad(angles_deg[k])
        glob = recover_ply_global(eps0, kappa, z_top, z_bot, qbars[k])
        e12_top, s12_top = recover_ply_local(glob["eps_top_xy"], glob["sig_top_xy"], theta)
        e12_bot, s12_bot = recover_ply_local(glob["eps_bot_xy"], glob["sig_bot_xy"], theta)

        entry = {
            "k": k,
            "theta_deg": angles_deg[k],
            "z_bot": z_bot, "z_top": z_top,
            "sig_top_xy": glob["sig_top_xy"], "sig_bot_xy": glob["sig_bot_xy"],
            "sig_top_12": s12_top,            "sig_bot_12": s12_bot,
            "eps_top_xy": glob["eps_top_xy"], "eps_bot_xy": glob["eps_bot_xy"],
            "eps_top_12": e12_top,            "eps_bot_12": e12_bot
        }

        if strengths is not None:
            Xt  = strengths["X_T"]; Xc  = strengths["X_C"]
            Yt  = strengths["Y_T"]; Yc  = strengths["Y_C"]
            s12 = strengths["S12"]
            entry["tsai_wu_top"] = tsai_wu(s12_top, Xt, Xc, Yt, Yc, s12)
            entry["tsai_wu_bot"] = tsai_wu(s12_bot, Xt, Xc, Yt, Yc, s12)
            h_top = hashin(s12_top, Xt, Xc, Yt, Yc, s12)
            h_bot = hashin(s12_bot, Xt, Xc, Yt, Yc, s12)
            entry["hashin_FT_top"] = h_top["FI_FT"]; entry["hashin_FC_top"] = h_top["FI_FC"]
            entry["hashin_MT_top"] = h_top["FI_MT"]; entry["hashin_MC_top"] = h_top["FI_MC"]
            entry["hashin_FT_bot"] = h_bot["FI_FT"]; entry["hashin_FC_bot"] = h_bot["FI_FC"]
            entry["hashin_MT_bot"] = h_bot["FI_MT"]; entry["hashin_MC_bot"] = h_bot["FI_MC"]

        ply_results.append(entry)
    return {"A":A, "B":B, "D":D, "eps0":eps0, "kappa":kappa, "plies": ply_results}

# --- Block 8: Navier center deflection for SSSS orthotropic ----
def navier_center_deflection(D: np.ndarray, a: float, b: float, q: float,
                             max_odd: int = 5) -> float:
    """
    Center deflection w(a/2,b/2) for SSSS orthotropic plate using truncated
    Navier double-sine series. Uses D11,D22,D12,D66 from laminate D.
    a,b in m, q in N/m^2. Returns w in meters.
    Sum over odd m,n up to max_odd (e.g., 5).
    """
    D11, D22, D12, D66 = D[0,0], D[1,1], D[0,1], D[2,2]
    w = 0.0
    for m in range(1, max_odd+1, 2):
        for n in range(1, max_odd+1, 2):
            mpa = (m*np.pi)/a
            npb = (n*np.pi)/b
            denom = D11*mpa**4 + 2.0*(D12 + 2.0*D66)*mpa**2*npb**2 + D22*npb**4
            w += (16.0*q)/(np.pi**6 * m**2 * n**2) * (1.0/denom)
    return w

# === PUBLIC API 1: laminate_abd (matches main.py) ==========
def laminate_abd(plies: List[Ply]) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Assemble ABD for an arbitrary laminate and return (A,B,D,z_interfaces).
    Ply.theta is radians. Ply.t is thickness (m).
    """
    n = len(plies)
    total_t = sum(p.t for p in plies)
    # z from -t/2 to +t/2 with segment lengths equal to each ply thickness
    z = np.empty(n + 1, dtype=float)
    z[0] = -0.5 * total_t
    for k in range(n):
        z[k+1] = z[k] + plies[k].t

    A = np.zeros((3,3)); B = np.zeros((3,3)); D = np.zeros((3,3))
    for k, p in enumerate(plies):
        Q  = Q_matrix(p.E1, p.E2, p.G12, p.v12)
        Qb = Q_bar(Q, p.theta)               # theta already radians
        zk, zk1 = z[k], z[k+1]
        dz = zk1 - zk
        A += Qb * dz
        B += 0.5  * Qb * (zk1**2 - zk**2)
        D += (1.0/3.0) * Qb * (zk1**3 - zk**3)
    return A, B, D, z


# === PUBLIC API 2: midplane_response (matches main.py) =====
def midplane_response(A: np.ndarray, B: np.ndarray, D: np.ndarray,
                      Nx: float = 0.0, Ny: float = 0.0, Nxy: float = 0.0,
                      Mx: float = 0.0, My: float = 0.0, Mxy: float = 0.0,
                      method: str = "schur") -> Tuple[np.ndarray, np.ndarray]:
    """
    Return (eps0, kappa) given ABD and in-plane loads and moments.
    Defaults to membrane-only if M* = 0.
    """
    N = np.array([Nx, Ny, Nxy], dtype=float)
    M = np.array([Mx, My, Mxy], dtype=float)
    return solve_midplane(A, B, D, N, M, method=method)



# === PUBLIC API 3: ply_strains_stresses (matches main.py) ==
def ply_strains_stresses(plies: List[Ply],
                         z_interfaces: np.ndarray,
                         eps0: np.ndarray,
                         kappa: np.ndarray) -> List[dict]:
    """
    Compute per-ply strains/stresses at top and bottom.
    Returns a list of dicts with global xy and local 12 values.
    """
    out = []
    for k, p in enumerate(plies):
        z_bot, z_top = z_interfaces[k], z_interfaces[k+1]
        Q  = Q_matrix(p.E1, p.E2, p.G12, p.v12)
        Qb = Q_bar(Q, p.theta)

        # global at top/bottom
        g = recover_ply_global(eps0, kappa, z_top, z_bot, Qb)
        # local transforms
        e12_top, s12_top = recover_ply_local(g["eps_top_xy"], g["sig_top_xy"], p.theta)
        e12_bot, s12_bot = recover_ply_local(g["eps_bot_xy"], g["sig_bot_xy"], p.theta)

        out.append({
            "k": k,
            "theta_rad": p.theta,
            "z_bot": z_bot, "z_top": z_top,
            "eps_top_xy": g["eps_top_xy"], "eps_bot_xy": g["eps_bot_xy"],
            "sig_top_xy": g["sig_top_xy"], "sig_bot_xy": g["sig_bot_xy"],
            "eps_top_12": e12_top,         "eps_bot_12": e12_bot,
            "sig_top_12": s12_top,         "sig_bot_12": s12_bot
        })
    return out

