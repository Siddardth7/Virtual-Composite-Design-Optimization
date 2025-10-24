"""
utils.py
General utilities for materials I/O and plotting.
"""
from __future__ import annotations
import pandas as pd
from pathlib import Path
from typing import Dict

DATA_DIR = Path(__file__).resolve().parents[1] / "data"

def load_materials(csv_path: Path | None = None) -> pd.DataFrame:
    """
    Load a materials CSV with columns:
    name, E1[Pa], E2[Pa], G12[Pa], v12[-], density[kg/m3], t_ply[m]
    """
    path = Path(csv_path) if csv_path else (DATA_DIR / "materials.csv")
    df = pd.read_csv(path)
    required = {"name","E1","E2","G12","v12","density","t_ply"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"materials.csv missing columns: {missing}")
    return df

def deg2rad(angle_deg: float) -> float:
    import math
    return angle_deg * math.pi / 180.0
