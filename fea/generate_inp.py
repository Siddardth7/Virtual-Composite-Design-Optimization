"""
fea/generate_inp.py
Part 1 — Mesh + Element generation for CalculiX (.inp) input file.
Part 2 — Material definition, shell composite section, BCs, load, analysis step.

Generates a complete CalculiX .inp file for a simply-supported composite
plate under uniform pressure, for validation against CLT Navier deflection.

Plate spec (from main.py / data/materials.csv):
  Lx = Ly = 0.300 m, q = 1 000 Pa
  Layup: [0, 45, -45, 90]_s  (IM7/8552, 1 mm total, 8 × 0.125 mm plies)
  CLT Navier centre deflection target: 0.064411 mm

Material — IM7/8552 (exact values from data/materials.csv):
  E1 = 161 GPa,  E2 = 11.4 GPa,  G12 = 5.17 GPa,  ν12 = 0.32,  ρ = 1600 kg/m³
  (E3, ν13 = E2, ν12 — transversely isotropic;  ν23 = 0.45 — IM7/8552 typical)

Element type: S4  (4-node bilinear shell, CalculiX)
Mesh density: nx × ny quads (default 20×20)

Node numbering convention (row-major, x varies fastest):
  node_id(i, j) = j * (nx+1) + i + 1   [1-indexed, i=0..nx, j=0..ny]

CCW element connectivity (viewed from +z):
  n1(i,j) → n2(i+1,j) → n3(i+1,j+1) → n4(i,j+1)

Outputs:
  fea/abaqus_inputs/mesh.inp      — Part 1 only (mesh + node sets)
  fea/abaqus_inputs/plate_ss.inp  — Complete model (mesh + material + BCs + step)
"""

from __future__ import annotations
import sys
from pathlib import Path
import numpy as np

# ── default plate / mesh parameters ─────────────────────────────────────────
DEFAULT_LX   = 0.300    # m  (plate x-length)
DEFAULT_LY   = 0.300    # m  (plate y-length)
DEFAULT_NX   = 20       # elements in x
DEFAULT_NY   = 20       # elements in y

# ── paths ────────────────────────────────────────────────────────────────────
FEA_DIR    = Path(__file__).resolve().parent
INP_DIR    = FEA_DIR / "abaqus_inputs"
INP_DIR.mkdir(parents=True, exist_ok=True)

# ── Part 2: material constants (IM7/8552 — source: data/materials.csv) ───────
MAT_NAME = "IM7_8552"

# Engineering constants — transversely isotropic assumption
E1   = 161.0e9    # Pa  longitudinal modulus
E2   = 11.4e9     # Pa  transverse modulus (in-plane)
E3   = 11.4e9     # Pa  transverse modulus (through-thickness) = E2
NU12 = 0.32       # —   major Poisson ratio
NU13 = 0.32       # —   = NU12 (transversely isotropic)
NU23 = 0.45       # —   transverse Poisson ratio (IM7/8552 typical)
G12  = 5.17e9     # Pa  in-plane shear modulus
G13  = 5.17e9     # Pa  = G12 (transversely isotropic)
G23  = E2 / (2.0 * (1.0 + NU23))   # Pa ≈ 3.931e9  (from NU23)
RHO  = 1600.0     # kg/m³  density

# ── Part 2: layup + load ─────────────────────────────────────────────────────
T_PLY      = 1.25e-4                        # m   (0.125 mm / ply — 8 plies = 1 mm)
LAYUP_DEG  = [0, 45, -45, 90, 90, -45, 45, 0]   # baseline [0/45/-45/90]_s
Q_PA       = 1_000.0                        # Pa  uniform pressure
PLATE_INP  = INP_DIR / "plate_ss.inp"       # complete model output


# ════════════════════════════════════════════════════════════════════════════
# Core mesh builder
# ════════════════════════════════════════════════════════════════════════════

def build_mesh(
    nx: int   = DEFAULT_NX,
    ny: int   = DEFAULT_NY,
    Lx: float = DEFAULT_LX,
    Ly: float = DEFAULT_LY,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Build node coordinates and element connectivity arrays.

    Returns
    -------
    nodes : (n_nodes, 4) float array
        Columns: [node_id, x, y, z]  (z = 0 for a flat mid-plane)
    elements : (n_elems, 5) int array
        Columns: [elem_id, n1, n2, n3, n4]  (1-indexed node IDs)
    """
    # ── node coordinates ─────────────────────────────────────────────────────
    xs = np.linspace(0.0, Lx, nx + 1)
    ys = np.linspace(0.0, Ly, ny + 1)

    n_nodes = (nx + 1) * (ny + 1)
    nodes = np.zeros((n_nodes, 4), dtype=float)

    idx = 0
    for j in range(ny + 1):          # y-direction outer loop
        for i in range(nx + 1):      # x-direction inner loop  (i varies fastest)
            nid = j * (nx + 1) + i + 1   # 1-indexed
            nodes[idx] = [nid, xs[i], ys[j], 0.0]
            idx += 1

    # ── element connectivity ─────────────────────────────────────────────────
    n_elems = nx * ny
    elements = np.zeros((n_elems, 5), dtype=int)

    eid = 0
    for je in range(ny):             # element row (y)
        for ie in range(nx):         # element col (x)
            n1 = je * (nx + 1) + ie + 1          # bottom-left
            n2 = je * (nx + 1) + (ie + 1) + 1    # bottom-right
            n3 = (je + 1) * (nx + 1) + (ie + 1) + 1  # top-right
            n4 = (je + 1) * (nx + 1) + ie + 1    # top-left
            elements[eid] = [eid + 1, n1, n2, n3, n4]
            eid += 1

    return nodes, elements


# ════════════════════════════════════════════════════════════════════════════
# Node-set builders
# ════════════════════════════════════════════════════════════════════════════

def build_node_sets(
    nodes: np.ndarray,
    nx: int,
    ny: int,
) -> dict[str, list[int]]:
    """
    Return named node sets needed for boundary conditions and output.

    Sets
    ----
    NALL        — every node
    EDGE_XMIN   — x = 0   (left edge)
    EDGE_XMAX   — x = Lx  (right edge)
    EDGE_YMIN   — y = 0   (bottom edge)
    EDGE_YMAX   — y = Ly  (top edge)
    CORNER_ALL  — four corner nodes
    CENTRE      — node nearest plate centre (for deflection output)
    """
    node_ids = nodes[:, 0].astype(int)    # all 1-indexed IDs

    # edges (by position in node numbering, not floats — more robust)
    edge_xmin = [j * (nx + 1) + 1         for j in range(ny + 1)]
    edge_xmax = [j * (nx + 1) + nx + 1    for j in range(ny + 1)]
    edge_ymin = [i + 1                     for i in range(nx + 1)]
    edge_ymax = [(ny) * (nx + 1) + i + 1  for i in range(nx + 1)]

    corners = [
        1,                          # (0, 0)
        nx + 1,                     # (Lx, 0)
        ny * (nx + 1) + 1,          # (0, Ly)
        ny * (nx + 1) + nx + 1,     # (Lx, Ly)
    ]

    # centre node — closest to (Lx/2, Ly/2)
    cx_idx = nx // 2
    cy_idx = ny // 2
    centre_nid = cy_idx * (nx + 1) + cx_idx + 1

    return {
        "NALL":      node_ids.tolist(),
        "EDGE_XMIN": edge_xmin,
        "EDGE_XMAX": edge_xmax,
        "EDGE_YMIN": edge_ymin,
        "EDGE_YMAX": edge_ymax,
        "CORNER_ALL": corners,
        "CENTRE":    [centre_nid],
    }


# ════════════════════════════════════════════════════════════════════════════
# .inp writers
# ════════════════════════════════════════════════════════════════════════════

def _write_heading(f, nx: int, ny: int, Lx: float, Ly: float) -> None:
    f.write("**\n")
    f.write("** CalculiX input — Part 1: Mesh + Elements\n")
    f.write(f"** Plate: {Lx*1e3:.1f} x {Ly*1e3:.1f} mm,  mesh: {nx}x{ny} S4 quads\n")
    f.write(f"** Nodes: {(nx+1)*(ny+1)},  Elements: {nx*ny}\n")
    f.write("** Generated by fea/generate_inp.py\n")
    f.write("**\n")


def _write_nodes(f, nodes: np.ndarray) -> None:
    f.write("*NODE, NSET=NALL\n")
    for row in nodes:
        nid = int(row[0])
        x, y, z = row[1], row[2], row[3]
        f.write(f"{nid:8d}, {x:18.10e}, {y:18.10e}, {z:18.10e}\n")
    f.write("**\n")


def _write_elements(f, elements: np.ndarray) -> None:
    f.write("*ELEMENT, TYPE=S4, ELSET=EALL\n")
    for row in elements:
        eid, n1, n2, n3, n4 = row
        f.write(f"{eid:8d}, {n1:8d}, {n2:8d}, {n3:8d}, {n4:8d}\n")
    f.write("**\n")


def _write_nsets(f, nsets: dict[str, list[int]]) -> None:
    for name, ids in nsets.items():
        f.write(f"*NSET, NSET={name}\n")
        # write 16 IDs per line (CalculiX allows up to 16)
        for chunk_start in range(0, len(ids), 16):
            chunk = ids[chunk_start:chunk_start + 16]
            f.write(", ".join(str(n) for n in chunk) + "\n")
    f.write("**\n")


# ════════════════════════════════════════════════════════════════════════════
# Public entry point
# ════════════════════════════════════════════════════════════════════════════

def generate_mesh_inp(
    out_path: Path | str | None = None,
    nx: int   = DEFAULT_NX,
    ny: int   = DEFAULT_NY,
    Lx: float = DEFAULT_LX,
    Ly: float = DEFAULT_LY,
) -> Path:
    """
    Generate CalculiX mesh .inp file (Part 1: nodes + elements + node sets).

    Parameters
    ----------
    out_path : path to write  (default: fea/abaqus_inputs/mesh.inp)
    nx, ny   : element counts in x and y
    Lx, Ly   : plate dimensions in metres

    Returns
    -------
    Path to the written file.
    """
    if out_path is None:
        out_path = INP_DIR / "mesh.inp"
    out_path = Path(out_path)

    nodes, elements = build_mesh(nx, ny, Lx, Ly)
    nsets = build_node_sets(nodes, nx, ny)

    with out_path.open("w") as f:
        _write_heading(f, nx, ny, Lx, Ly)
        _write_nodes(f, nodes)
        _write_elements(f, elements)
        _write_nsets(f, nsets)

    return out_path


# ════════════════════════════════════════════════════════════════════════════
# Part 2 — Material / Section / BCs / Step writers
# ════════════════════════════════════════════════════════════════════════════

def _write_material(f) -> None:
    """
    Write *MATERIAL block for IM7/8552.

    CalculiX ENGINEERING CONSTANTS format (9 values, split over 2 lines):
      Line 1 :  E1  E2  E3  NU12  NU13  NU23  G12  G13
      Line 2 :  G23
    """
    f.write("**\n")
    f.write("** ── Material: IM7/8552 (engineering constants) ─────────────────\n")
    f.write(f"*MATERIAL, NAME={MAT_NAME}\n")
    f.write("*ELASTIC, TYPE=ENGINEERING CONSTANTS\n")
    f.write(
        f" {E1:.6e}, {E2:.6e}, {E3:.6e},"
        f" {NU12:.4f}, {NU13:.4f}, {NU23:.4f},"
        f" {G12:.6e}, {G13:.6e}\n"
    )
    f.write(f" {G23:.6e}\n")
    f.write("*DENSITY\n")
    f.write(f" {RHO:.1f}\n")
    f.write("**\n")


def _write_shell_section(
    f,
    layup: list[int] | None = None,
    t_ply: float = T_PLY,
) -> None:
    """
    Write *SHELL SECTION, COMPOSITE for the baseline layup.

    Each ply line:  thickness, n_integration_pts, material_name, angle_deg
    n_int = 3  (Simpson rule — CalculiX default for composites)
    Ply ordering: first entry = bottom ply (most negative z).
    """
    if layup is None:
        layup = LAYUP_DEG
    f.write("**\n")
    f.write("** ── Shell composite section (bottom → top) ─────────────────────\n")
    f.write("*SHELL SECTION, ELSET=EALL, COMPOSITE, OFFSET=MID\n")
    n_int = 3
    for angle in layup:
        f.write(f" {t_ply:.6e}, {n_int}, {MAT_NAME}, {float(angle):.1f}\n")
    f.write("**\n")


def _write_boundary_conditions(f) -> None:
    """
    Write simply-supported BCs consistent with the Navier SSSS solution.

    Navier BCs for a plate in the x-y plane, pressure in -z direction:
      x-edges (XMIN, XMAX) : u = 0  (in-plane normal),  w = 0  (transverse)
      y-edges (YMIN, YMAX) : v = 0  (in-plane normal),  w = 0  (transverse)
    Rotational DOFs 4 (UR1) and 5 (UR2) are left FREE → zero edge moments.

    CalculiX DOF convention for shells:
      1 = U1 (x),  2 = U2 (y),  3 = U3 (z = w)
      4 = UR1,     5 = UR2,     6 = UR3
    """
    f.write("**\n")
    f.write("** ── Boundary conditions — simply-supported SSSS (Navier) ────────\n")
    f.write("*BOUNDARY\n")
    for nset in ("EDGE_XMIN", "EDGE_XMAX"):
        f.write(f"  {nset}, 1, 1, 0.0\n")   # u = 0
        f.write(f"  {nset}, 3, 3, 0.0\n")   # w = 0
    for nset in ("EDGE_YMIN", "EDGE_YMAX"):
        f.write(f"  {nset}, 2, 2, 0.0\n")   # v = 0
        f.write(f"  {nset}, 3, 3, 0.0\n")   # w = 0
    f.write("**\n")


def _write_step(f, q: float = Q_PA) -> None:
    """
    Write *STEP block: linear static analysis + pressure load + output requests.

    Output files generated by CalculiX:
      plate_ss.frd  — binary results (U, S, E) for CGX post-processing
      plate_ss.dat  — text dump of CENTRE node displacement (for Python parsing)
    """
    f.write("**\n")
    f.write("** ── Analysis step: linear static ───────────────────────────────\n")
    f.write("*STEP, NLGEOM=NO\n")
    f.write("*STATIC\n")
    f.write("**\n")
    f.write(f"** Uniform pressure  q = {q:.1f} Pa  (positive = into plate, -z direction)\n")
    f.write("*DLOAD\n")
    f.write(f"  EALL, P, {q:.4e}\n")
    f.write("**\n")
    f.write("** .frd file output (CGX / post-processing)\n")
    f.write("*NODE FILE\n")
    f.write("  U\n")
    f.write("*EL FILE\n")
    f.write("  S, E\n")
    f.write("**\n")
    f.write("** .dat file output (text — for Python post-processing)\n")
    f.write("*NODE PRINT, NSET=CENTRE\n")
    f.write("  U\n")
    f.write("*EL PRINT, ELSET=EALL\n")
    f.write("  S\n")
    f.write("**\n")
    f.write("*END STEP\n")


# ════════════════════════════════════════════════════════════════════════════
# Part 2 public entry point
# ════════════════════════════════════════════════════════════════════════════

def generate_full_inp(
    out_path: Path | str | None = None,
    nx: int   = DEFAULT_NX,
    ny: int   = DEFAULT_NY,
    Lx: float = DEFAULT_LX,
    Ly: float = DEFAULT_LY,
    layup: list[int] | None = None,
    t_ply: float = T_PLY,
    q: float = Q_PA,
) -> Path:
    """
    Generate the complete CalculiX .inp file for the composite plate.

    Combines Part 1 (mesh + node sets) with Part 2 (material definition,
    composite shell section, simply-supported BCs, and linear-static step).

    Parameters
    ----------
    out_path : output file path  (default: fea/abaqus_inputs/plate_ss.inp)
    nx, ny   : element counts in x and y
    Lx, Ly   : plate dimensions in metres
    layup    : ply angles in degrees, bottom-to-top  (default: LAYUP_DEG)
    t_ply    : ply thickness in metres  (default: T_PLY = 0.125 mm)
    q        : uniform pressure in Pa  (default: Q_PA = 1000 Pa)

    Returns
    -------
    Path to the written file.
    """
    if out_path is None:
        out_path = PLATE_INP
    if layup is None:
        layup = LAYUP_DEG
    out_path = Path(out_path)

    nodes, elements = build_mesh(nx, ny, Lx, Ly)
    nsets = build_node_sets(nodes, nx, ny)

    with out_path.open("w") as f:
        # ── Part 1: mesh ────────────────────────────────────────────────────
        _write_heading(f, nx, ny, Lx, Ly)
        _write_nodes(f, nodes)
        _write_elements(f, elements)
        _write_nsets(f, nsets)
        # ── Part 2: material + section + BCs + step ─────────────────────────
        _write_material(f)
        _write_shell_section(f, layup, t_ply)
        _write_boundary_conditions(f)
        _write_step(f, q)

    return out_path


# ════════════════════════════════════════════════════════════════════════════
# CLI
# ════════════════════════════════════════════════════════════════════════════

def _report(nodes: np.ndarray, elements: np.ndarray, nsets: dict, out: Path) -> None:
    print(f"  Written : {out}")
    print(f"  Nodes   : {len(nodes)}")
    print(f"  Elements: {len(elements)}")
    print()
    # spot-check corners
    corners = [0, -1]
    print("  Corner nodes (first / last):")
    for c in corners:
        nid = int(nodes[c, 0])
        x, y = nodes[c, 1], nodes[c, 2]
        print(f"    Node {nid:5d}  x={x*1e3:8.3f} mm  y={y*1e3:8.3f} mm")
    # centre node
    cid = nsets["CENTRE"][0]
    cn  = nodes[cid - 1]
    print(f"\n  Centre node {cid}: x={cn[1]*1e3:.3f} mm  y={cn[2]*1e3:.3f} mm")
    print()
    print("  Node set sizes:")
    for k, v in nsets.items():
        print(f"    {k:<12s}: {len(v)}")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Generate CalculiX .inp files for composite plate (Parts 1 + 2)"
    )
    parser.add_argument("--nx",     type=int,   default=DEFAULT_NX,  help="elements in x")
    parser.add_argument("--ny",     type=int,   default=DEFAULT_NY,  help="elements in y")
    parser.add_argument("--Lx",     type=float, default=DEFAULT_LX,  help="plate x-length [m]")
    parser.add_argument("--Ly",     type=float, default=DEFAULT_LY,  help="plate y-length [m]")
    parser.add_argument("--out",    type=str,   default=None,        help="mesh.inp path override")
    parser.add_argument("--full",   type=str,   default=None,        help="plate_ss.inp path override")
    parser.add_argument("--mesh-only", action="store_true",
                        help="write mesh.inp only (skip Part 2)")
    args = parser.parse_args()

    print(f"\n=== generate_inp.py  Parts 1 + 2 ===")
    print(f"  Plate    : {args.Lx*1e3:.0f} x {args.Ly*1e3:.0f} mm")
    print(f"  Mesh     : {args.nx} x {args.ny}  ({args.nx*args.ny} S4 elements)")
    print(f"  Layup    : {LAYUP_DEG}")
    print(f"  t_ply    : {T_PLY*1e3:.3f} mm   total: {len(LAYUP_DEG)*T_PLY*1e3:.3f} mm")
    print(f"  Pressure : {Q_PA:.0f} Pa\n")

    nodes, elements = build_mesh(args.nx, args.ny, args.Lx, args.Ly)
    nsets = build_node_sets(nodes, args.nx, args.ny)

    # ── Part 1: mesh only ──────────────────────────────────────────────────
    mesh_out = generate_mesh_inp(args.out, args.nx, args.ny, args.Lx, args.Ly)
    print("── Part 1 (mesh.inp) ──────────────────────────────")
    _report(nodes, elements, nsets, mesh_out)

    if args.mesh_only:
        sys.exit(0)

    # ── Parts 1+2: complete model ──────────────────────────────────────────
    full_out = generate_full_inp(args.full, args.nx, args.ny, args.Lx, args.Ly)
    print("\n── Part 2 (plate_ss.inp) ──────────────────────────")
    print(f"  Written  : {full_out}")
    print(f"  Material : {MAT_NAME}")
    print(f"    E1={E1/1e9:.1f} GPa  E2={E2/1e9:.1f} GPa  G12={G12/1e9:.2f} GPa  ν12={NU12}")
    print(f"    G23={G23/1e9:.3f} GPa  ν23={NU23}  ρ={RHO} kg/m³")
    print(f"  Section  : COMPOSITE, OFFSET=MID, {len(LAYUP_DEG)} plies × {T_PLY*1e3:.3f} mm")
    print(f"  BCs      : SSSS (Navier — u/w=0 on x-edges, v/w=0 on y-edges)")
    print(f"  Load     : DLOAD P = {Q_PA:.0f} Pa on EALL")
    print(f"\n  Run with:  ccx plate_ss")
    print(f"  (from:     fea/abaqus_inputs/)")
    print()
