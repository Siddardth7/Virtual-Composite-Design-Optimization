"""
utils.py
General utilities for materials I/O and angles.
Normalizes materials.csv columns to:
name, E1, E2, G12, v12, density, t_ply  (all SI)
"""
from __future__ import annotations
import pandas as pd
from pathlib import Path
import math

DATA_DIR = Path(__file__).resolve().parents[1] / "data"

# allow common header variants → canonical names
_COL_MAPS = [
    # canonical : possible variants in CSV
    ("name",    ["name", "material", "material_name"]),
    ("E1",      ["E1", "E1_Pa"]),
    ("E2",      ["E2", "E2_Pa"]),
    ("G12",     ["G12", "G12_Pa"]),
    ("v12",     ["v12", "nu12", "nu_12", "nu12[-]"]),
    ("density", ["density", "rho", "rho_kgm3", "rho_[kg/m3]"]),
    ("t_ply",   ["t_ply", "ply_t", "ply_t_m", "tply_m"]),
]

_REQUIRED = {k for k, _ in _COL_MAPS}

def _normalize_columns(df: pd.DataFrame) -> pd.DataFrame:
    cols = {c: c.strip() for c in df.columns}
    df = df.rename(columns=cols)
    # build rename dict by first matching variant present
    rename = {}
    for canon, variants in _COL_MAPS:
        for v in variants:
            if v in df.columns:
                rename[v] = canon
                break
    df = df.rename(columns=rename)
    # check required present
    missing = _REQUIRED - set(df.columns)
    if missing:
        raise ValueError(f"materials.csv missing columns: {missing}")
    return df

def _to_float(df: pd.DataFrame, col: str) -> pd.Series:
    # handles numeric strings like "161e9"
    s = pd.to_numeric(df[col], errors="coerce")
    if s.isna().any():
        bad = df.loc[s.isna(), col].tolist()[:5]
        raise ValueError(f"Non-numeric values in column '{col}': {bad}")
    return s.astype(float)

def load_materials(csv_path: Path | None = None) -> pd.DataFrame:
    """
    Load and normalize a materials CSV into canonical columns:
      name, E1, E2, G12, v12, density, t_ply   (SI units)
    Accepts variants like E1_Pa, rho_kgm3, ply_t_m, etc., and renames them.
    """
    path = Path(csv_path) if csv_path else (DATA_DIR / "materials.csv")
    df = pd.read_csv(path)
    df = _normalize_columns(df)

    # convert to float
    for c in ["E1", "E2", "G12", "v12", "density", "t_ply"]:
        df[c] = _to_float(df, c)

    return df[["name","E1","E2","G12","v12","density","t_ply"]]

def deg2rad(angle_deg: float) -> float:
    return math.radians(angle_deg)

def rad2deg(angle_rad: float) -> float:
    return math.degrees(angle_rad)
