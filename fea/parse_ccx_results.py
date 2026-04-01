"""
fea/parse_ccx_results.py
Parse CalculiX .dat text output for composite plate validation.

Extracts:
  - Centre-node U3 deflection (w at plate centre, in metres)
  - S11 (σₓₓ) at element integration points (max absolute value across all
    plies / integration points → represents peak in-plane normal stress)

CalculiX .dat format reference (relevant sections):
  *  NODE PRINT block — one line per node per variable:
       <node_id>   <U1>   <U2>   <U3>
  *  EL PRINT block — one line per integration point per variable,
     preceded by a header that identifies the element set and variable name.
     S11 appears in the column labelled SXX in the stress output.

Usage (CLI):
  python fea/parse_ccx_results.py fea/results/plate_ss.dat

Returns (when imported):
  parse_dat(path) -> dict with keys:
      deflection_m    : float  — magnitude |U3| at CENTRE node (m)
      sigma_xx_max_Pa : float  — max |S11| across all output points (Pa)
      notes           : str
"""

from __future__ import annotations
import re
import sys
from pathlib import Path


# ─────────────────────────────────────────────────────────────────────────────
# Internal helpers
# ─────────────────────────────────────────────────────────────────────────────

def _is_float(s: str) -> bool:
    try:
        float(s)
        return True
    except ValueError:
        return False


def _parse_node_print(lines: list[str], start: int) -> dict[int, list[float]]:
    """
    Scan forward from *start* collecting node-data rows until the next keyword
    or blank-header line.  Returns {node_id: [U1, U2, U3, ...]}
    """
    data: dict[int, list[float]] = {}
    i = start
    while i < len(lines):
        line = lines[i].strip()
        if line.startswith("*") or line == "":
            i += 1
            if line.startswith("*") and not line.startswith("** "):
                break
            continue
        parts = line.split()
        if parts and parts[0].isdigit():
            nid = int(parts[0])
            vals = [float(v) for v in parts[1:] if _is_float(v)]
            if vals:
                data[nid] = vals
        i += 1
    return data


def _parse_el_print(lines: list[str], start: int) -> list[float]:
    """
    Scan forward collecting numeric data rows (element / integration-point output).
    Returns a flat list of all floats found (column ordering preserved per row).
    The S11 column position is determined by the preceding column header.
    """
    s11_col: int | None = None
    s11_values: list[float] = []
    i = start

    while i < len(lines):
        line = lines[i]
        stripped = line.strip()

        # detect header row for stress output (contains SXX or S11-like labels)
        if re.search(r"\bS11\b|\bSXX\b", stripped, re.IGNORECASE):
            # find which column (0-indexed in whitespace-split tokens after first)
            header_parts = stripped.split()
            for col_idx, token in enumerate(header_parts):
                if re.match(r"S11|SXX", token, re.IGNORECASE):
                    s11_col = col_idx  # 0-indexed in the header tokens
                    break
            i += 1
            continue

        # stop at next keyword section
        if stripped.startswith("*") and not stripped.startswith("** "):
            break

        parts = stripped.split()
        if len(parts) >= 2 and parts[0].isdigit() and s11_col is not None:
            # row: elem_id [ip_id] val1 val2 ...
            # data tokens start after element (and optional integration-point) ids
            numeric_start = 0
            for idx, p in enumerate(parts):
                if not p.isdigit():
                    numeric_start = idx
                    break
            vals = parts[numeric_start:]
            if s11_col < len(vals):
                try:
                    s11_values.append(float(vals[s11_col]))
                except ValueError:
                    pass
        i += 1

    return s11_values


# ─────────────────────────────────────────────────────────────────────────────
# Public API
# ─────────────────────────────────────────────────────────────────────────────

def parse_dat(dat_path: str | Path) -> dict:
    """
    Parse a CalculiX .dat file and return a results dict.

    Parameters
    ----------
    dat_path : path to the .dat file produced by CalculiX

    Returns
    -------
    dict with keys:
        deflection_m      : float   |U3| at CENTRE node (metres)
        sigma_xx_max_Pa   : float   max |S11| over all output points (Pa)
        centre_node_id    : int     node ID of the CENTRE output node
        notes             : str
    """
    path = Path(dat_path)
    if not path.exists():
        raise FileNotFoundError(f"CalculiX .dat file not found: {path}")

    with path.open() as f:
        lines = f.readlines()

    deflection_m: float | None = None
    sigma_xx_max_Pa: float | None = None
    centre_node_id: int | None = None
    all_s11: list[float] = []

    i = 0
    while i < len(lines):
        line = lines[i].strip()

        # ── NODE PRINT section ─────────────────────────────────────────────
        if re.match(r"\*NODE PRINT", line, re.IGNORECASE):
            node_data = _parse_node_print(lines, i + 1)
            for nid, vals in node_data.items():
                # U3 is the 3rd displacement component (index 2)
                if len(vals) >= 3:
                    u3 = vals[2]
                    if deflection_m is None or abs(u3) > abs(deflection_m):
                        deflection_m = u3
                        centre_node_id = nid

        # ── EL PRINT section ──────────────────────────────────────────────
        elif re.match(r"\*EL PRINT", line, re.IGNORECASE):
            s11_vals = _parse_el_print(lines, i + 1)
            all_s11.extend(s11_vals)

        i += 1

    if all_s11:
        sigma_xx_max_Pa = max(all_s11, key=abs)

    notes = "parsed from .dat"
    if deflection_m is None:
        notes += "; WARN: no U3 found"
    if sigma_xx_max_Pa is None:
        notes += "; WARN: no S11 found"

    return {
        "deflection_m":    float(deflection_m)    if deflection_m    is not None else None,
        "sigma_xx_max_Pa": float(sigma_xx_max_Pa) if sigma_xx_max_Pa is not None else None,
        "centre_node_id":  centre_node_id,
        "notes":           notes,
    }


# ─────────────────────────────────────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────────────────────────────────────

def _cli_main(argv: list[str]) -> None:
    if len(argv) < 2:
        print("Usage: python fea/parse_ccx_results.py <path/to/plate_ss.dat>")
        sys.exit(1)

    dat_path = Path(argv[1])
    print(f"\n=== parse_ccx_results.py ===")
    print(f"  File: {dat_path}\n")

    try:
        r = parse_dat(dat_path)
    except FileNotFoundError as exc:
        print(f"  ERROR: {exc}")
        sys.exit(1)

    defl = r["deflection_m"]
    s11  = r["sigma_xx_max_Pa"]

    print(f"  Centre-node deflection  U3 = {defl*1e3:+.6f} mm   (node {r['centre_node_id']})")
    print(f"  Peak σₓₓ              S11 = {s11/1e6:+.3f} MPa")
    print(f"  Notes: {r['notes']}")

    # Sanity checks
    if defl is not None and defl >= 0:
        print("\n  ⚠  WARNING: U3 is non-negative (expected downward = negative).")
    if s11 is not None and s11 < 0:
        print("\n  ⚠  WARNING: peak S11 is negative — check ply surface (expect tension on bottom).")
    print()


if __name__ == "__main__":
    _cli_main(sys.argv)
