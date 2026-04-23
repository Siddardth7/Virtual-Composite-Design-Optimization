import sys
from pathlib import Path
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from utils import load_materials

ROOT = Path(__file__).resolve().parents[1]


class TestLoadMaterials:
    def test_loads_without_error(self):
        df = load_materials(ROOT / "data" / "materials.csv")
        assert len(df) >= 1

    def test_canonical_columns_present(self):
        df = load_materials(ROOT / "data" / "materials.csv")
        for col in ["name", "E1", "E2", "G12", "v12", "density", "t_ply"]:
            assert col in df.columns, f"Missing column: {col}"

    def test_strength_columns_present(self):
        df = load_materials(ROOT / "data" / "materials.csv")
        for col in ["X_T", "X_C", "Y_T", "Y_C", "S12"]:
            assert col in df.columns, f"Missing strength column: {col}"

    def test_im7_8552_e1_value(self):
        df = load_materials(ROOT / "data" / "materials.csv")
        im7 = df[df["name"] == "IM7_8552"].iloc[0]
        assert abs(im7["E1"] - 161e9) / 161e9 < 1e-3

    def test_all_values_positive(self):
        df = load_materials(ROOT / "data" / "materials.csv")
        for col in ["E1", "E2", "G12", "v12", "density", "t_ply"]:
            assert (df[col] > 0).all(), f"Non-positive value in {col}"
