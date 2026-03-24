"""
compare_clt_fea.py
Placeholder hooks to import an FEA baseline for comparison.
Populate `read_fea_results` to parse solver outputs.
"""
from __future__ import annotations
import pandas as pd
from pathlib import Path

def read_fea_results(results_dir: str | Path) -> pd.DataFrame:
    """
    Expected columns: case_id, deflection_m, sigma_xx_max_Pa, notes
    """
    p = Path(results_dir)
    csv = p / "fea_summary.csv"
    if not csv.exists():
        # Create placeholder
        df = pd.DataFrame([{"case_id":"baseline","deflection_m":None,"sigma_xx_max_Pa":None,"notes":"TODO"}])
        df.to_csv(csv, index=False)
    return pd.read_csv(csv)
