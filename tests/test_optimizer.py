import sys
from pathlib import Path
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))


class TestOptimizerMaterial:
    def test_optimizer_uses_correct_material(self):
        """SA optimizer material must match materials.csv IM7/8552 values."""
        import layup_optimizer_sa as sa
        from utils import load_materials
        ROOT = Path(__file__).resolve().parents[1]
        mat = load_materials(ROOT / "data" / "materials.csv").iloc[0]
        assert abs(sa.E1 - float(mat["E1"])) / float(mat["E1"]) < 1e-3, (
            f"SA E1={sa.E1} != materials.csv E1={mat['E1']}"
        )
        assert abs(sa.E2 - float(mat["E2"])) / float(mat["E2"]) < 1e-3
        assert abs(sa.t - float(mat["t_ply"])) / float(mat["t_ply"]) < 1e-3

    def test_optimizer_uses_clt_engine(self):
        """SA optimizer must import Q_bar from clt.py, not define its own."""
        src_file = Path(__file__).resolve().parents[1] / "src" / "layup_optimizer_sa.py"
        content = src_file.read_text()
        assert "from clt import" in content, "SA optimizer must use clt.py for CLT computations"

    def test_symmetric_laminate_from_sa(self):
        """Full sequence returned by SA must be symmetric."""
        import random
        import layup_optimizer_sa as sa
        random.seed(0)
        seq, _ = sa.simulated_annealing(n_iterations=500)
        n = len(seq)
        half = n // 2
        assert seq[:half] == seq[half:][::-1], "SA result is not symmetric"
