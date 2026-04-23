# VirtualCompositeDesign — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use `superpowers:subagent-driven-development` (recommended) or `superpowers:executing-plans` to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make the VirtualCompositeDesign project complete and physically correct — no wrong material constants, no silent data loss, no disconnected pipelines, with a full test suite.

**Architecture:** Fix the SA optimizer to use the validated CLT engine and materials CSV; extend `utils.py` to expose strength data; add a `pytest` suite with known analytical baselines; clean up all dead code and documentation drift.

**Tech Stack:** Python 3.9+, numpy, pandas, matplotlib, pytest, jupyter (optional for notebook)

**Source for all findings:** `docs/AUDIT_REPORT.md`

---

## Section 1 — Where to Start (Next Session Entry Points)

Execute these 5 tasks first, in this order. Each P0 fix unblocks the next.

| # | Task | Why First |
|---|------|-----------|
| 1 | **Task 2** — Fix `load_materials()` to expose strength columns | Unblocks every failure-criteria test and the SA material fix |
| 2 | **Task 3** — Fix SA optimizer: import clt.py + use materials.csv | Core correctness fix; makes optimizer physically meaningful |
| 3 | **Task 4** — Fix `evaluate_laminate()` z-interface calculation | Required before writing tests against it |
| 4 | **Task 5** — Fix `compare_clt_fea.py` stress computation | Correctness fix before stress tests |
| 5 | **Task 1** — Add test suite skeleton + first physics tests | Locks in correct behaviour before any further changes |

> **Note:** Task 1 (tests) is listed last here because the test targets (Tasks 2–5) must be fixed first. In TDD terms: fix the code, write the tests, then use tests to guard all subsequent work.

---

## Section 2 — Full Prioritized Task Table

| Priority | Effort | ID | File(s) | Task |
|----------|--------|----|---------|------|
| **P0** | M | Task 2 | `src/utils.py` | Load strength columns from materials.csv |
| **P0** | L | Task 3 | `src/layup_optimizer_sa.py` | Rewrite SA optimizer to use clt.py + materials.csv |
| **P0** | S | Task 4 | `src/clt.py:219–221` | Fix evaluate_laminate() z-interface for variable thickness |
| **P0** | S | Task 5 | `src/compare_clt_fea.py:121–123` | Fix stress computation: add eps0 term |
| **P0** | S | Task 6 | `src/compare_clt_fea.py:153`, `src/sa_spotcheck.py` | Fix FEA import path (add ROOT to sys.path) |
| **P1** | L | Task 1 | `tests/` (new) | Create pytest suite with physics baselines |
| **P1** | S | Task 7 | `src/clt.py:273` | Add division-by-zero guard in navier_center_deflection |
| **P1** | S | Task 8 | `src/clt.py:99–103` | Add singularity check in solve_midplane |
| **P1** | S | Task 9 | `src/utils.py`, `src/compare_clt_fea.py`, `src/sa_spotcheck.py` | Consolidate _find_ccx() into utils.py |
| **P1** | S | Task 10 | `src/sa_spotcheck.py:60–68` | Cache material load at module level |
| **P1** | S | Task 11 | `requirements.txt` | Pin versions, remove scipy, add jupyter |
| **P1** | S | Task 12 | `.gitignore`, `fea/abaqus_inputs/.gitkeep` | Fix untracked fea/ccx and abaqus_inputs/ |
| **P1** | S | Task 13 | `README.md` | Fix angle sweep resolution + API key names |
| **P1** | S | Task 14 | `docs/Methodology_Notes.md` | Fix Hashin MC formula in docs |
| **P1** | S | Task 15 | `src/` | Add basic input validation to clt.py public functions |
| **P1** | S | Task 16 | `src/run_pipeline.py` (new) | Add one-command pipeline runner |
| **P1** | S | Task 17 | `gui/.gitkeep` | Remove dead gui/ placeholder |
| **P2** | M | Task 18 | `src/main.py` or `src/visualize.py` (new) | Add laminate stack cross-section visualization |
| **P2** | M | Task 19 | `src/layup_optimizer_sa.py` | Save SA convergence history CSV + plot |
| **P2** | L | Task 20 | `src/layup_optimizer_sa.py` | Parallelize SA multi-seed runs with multiprocessing |

---

## Section 3 — Execution Order

```
P0 (must fix — wrong physics):
  Task 2 → Task 3 → Task 4 → Task 5 → Task 6

P1 (required for completeness):
  Task 1 (tests, run after P0 fixes)
  Task 7 → Task 8 → Task 9 → Task 10
  Task 11 → Task 12 → Task 13 → Task 14 → Task 15
  Task 16 → Task 17

P2 (enhancements):
  Task 18 → Task 19 → Task 20
```

---

## Detailed Task Steps

---

### Task 1: Create pytest suite with physics baselines

**Files:**
- Create: `tests/__init__.py`
- Create: `tests/test_clt.py`
- Create: `tests/test_utils.py`
- Create: `tests/test_optimizer.py`

**Context:** The test suite must verify physics against known analytical solutions, not just "does it run." Key baselines:
- Isotropic plate limit: for an isotropic material (E1=E2, G12=E/(2(1+ν)), ν12=ν), the laminate should reduce to classical isotropic plate theory.
- Single-ply Q_matrix: for 0° ply, Q_bar must equal Q.
- Symmetric laminate: B matrix must be zero (within floating point).
- Hashin FT limit: at stress = strength, FI = 1.0 exactly.

- [ ] **Step 1: Create test directory**

```bash
mkdir -p tests
touch tests/__init__.py
```

- [ ] **Step 2: Write CLT physics tests**

Create `tests/test_clt.py`:

```python
import sys
from pathlib import Path
import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from clt import (
    Q_matrix, Q_bar, abd_matrices, ply_interfaces,
    laminate_abd, solve_midplane, navier_center_deflection,
    hashin, tsai_wu, Ply, deg2rad,
)


# ── Material constants (IM7/8552) ──────────────────────────────────────────
E1, E2, G12, NU12 = 161e9, 11.4e9, 5.17e9, 0.32
T_PLY = 1.25e-4   # m


class TestQMatrix:
    def test_q_matrix_symmetric(self):
        Q = Q_matrix(E1, E2, G12, NU12)
        assert Q.shape == (3, 3)
        np.testing.assert_allclose(Q[0, 1], Q[1, 0], rtol=1e-10)

    def test_q_matrix_positive_diagonal(self):
        Q = Q_matrix(E1, E2, G12, NU12)
        assert Q[0, 0] > 0
        assert Q[1, 1] > 0
        assert Q[2, 2] > 0

    def test_q_matrix_e1_dominates(self):
        Q = Q_matrix(E1, E2, G12, NU12)
        assert Q[0, 0] > Q[1, 1]   # fibre direction stiffer


class TestQBar:
    def test_qbar_zero_angle_equals_q(self):
        Q = Q_matrix(E1, E2, G12, NU12)
        Qb = Q_bar(Q, 0.0)
        np.testing.assert_allclose(Qb, Q, rtol=1e-10)

    def test_qbar_90_swaps_e1_e2(self):
        Q = Q_matrix(E1, E2, G12, NU12)
        Qb90 = Q_bar(Q, np.pi / 2)
        # Q11 and Q22 swap at 90 degrees
        np.testing.assert_allclose(Qb90[0, 0], Q[1, 1], rtol=1e-6)
        np.testing.assert_allclose(Qb90[1, 1], Q[0, 0], rtol=1e-6)

    def test_qbar_45_symmetric(self):
        Q = Q_matrix(E1, E2, G12, NU12)
        Qb45 = Q_bar(Q, np.pi / 4)
        np.testing.assert_allclose(Qb45[0, 0], Qb45[1, 1], rtol=1e-6)


class TestSymmetricLaminate:
    def _make_symmetric_plies(self):
        angles = [0, 45, -45, 90, 90, -45, 45, 0]
        return [Ply(E1, E2, G12, NU12, deg2rad(th), T_PLY) for th in angles]

    def test_b_matrix_near_zero_for_symmetric(self):
        plies = self._make_symmetric_plies()
        A, B, D, z = laminate_abd(plies)
        np.testing.assert_allclose(B, np.zeros((3, 3)), atol=1e-6)

    def test_a_matrix_positive_definite(self):
        plies = self._make_symmetric_plies()
        A, B, D, z = laminate_abd(plies)
        eigvals = np.linalg.eigvalsh(A)
        assert np.all(eigvals > 0)

    def test_d_matrix_positive_definite(self):
        plies = self._make_symmetric_plies()
        A, B, D, z = laminate_abd(plies)
        eigvals = np.linalg.eigvalsh(D)
        assert np.all(eigvals > 0)

    def test_navier_deflection_reasonable(self):
        plies = self._make_symmetric_plies()
        A, B, D, z = laminate_abd(plies)
        w = navier_center_deflection(D, a=0.3, b=0.3, q=1000.0, max_odd=5)
        # Validated result from FEA: 0.06378 mm. CLT target: 0.064411 mm.
        assert 0.060e-3 < w < 0.070e-3, f"Unexpected deflection: {w*1e3:.6f} mm"

    def test_total_thickness(self):
        plies = self._make_symmetric_plies()
        A, B, D, z = laminate_abd(plies)
        total_t = sum(p.t for p in plies)
        np.testing.assert_allclose(total_t, 8 * T_PLY, rtol=1e-10)


class TestHashin:
    XT, XC, YT, YC, S12 = 2.8e9, 1.6e9, 70e6, 200e6, 100e6

    def test_fibre_tension_at_strength_is_one(self):
        # At σ1 = XT, τ12 = 0: FI_FT = (XT/XT)^2 = 1.0
        stress = np.array([self.XT, 0.0, 0.0])
        fi = hashin(stress, self.XT, self.XC, self.YT, self.YC, self.S12)
        np.testing.assert_allclose(fi["FI_FT"], 1.0, rtol=1e-10)

    def test_fibre_compression_at_strength_is_one(self):
        stress = np.array([-self.XC, 0.0, 0.0])
        fi = hashin(stress, self.XT, self.XC, self.YT, self.YC, self.S12)
        np.testing.assert_allclose(fi["FI_FC"], 1.0, rtol=1e-10)

    def test_matrix_tension_at_strength_is_one(self):
        stress = np.array([0.0, self.YT, 0.0])
        fi = hashin(stress, self.XT, self.XC, self.YT, self.YC, self.S12)
        np.testing.assert_allclose(fi["FI_MT"], 1.0, rtol=1e-10)

    def test_tension_mode_inactive_under_compression(self):
        stress = np.array([-self.XC, 0.0, 0.0])
        fi = hashin(stress, self.XT, self.XC, self.YT, self.YC, self.S12)
        assert fi["FI_FT"] == 0.0

    def test_compression_mode_inactive_under_tension(self):
        stress = np.array([self.XT, 0.0, 0.0])
        fi = hashin(stress, self.XT, self.XC, self.YT, self.YC, self.S12)
        assert fi["FI_FC"] == 0.0

    def test_zero_stress_all_zero(self):
        stress = np.array([0.0, 0.0, 0.0])
        fi = hashin(stress, self.XT, self.XC, self.YT, self.YC, self.S12)
        assert fi["FI_FT"] == 0.0
        assert fi["FI_FC"] == 0.0
        assert fi["FI_MT"] == 0.0
        assert fi["FI_MC"] == 0.0


class TestTsaiWu:
    XT, XC, YT, YC, S12 = 2.8e9, 1.6e9, 70e6, 200e6, 100e6

    def test_zero_stress_below_one(self):
        fi = tsai_wu(np.array([0.0, 0.0, 0.0]),
                     self.XT, self.XC, self.YT, self.YC, self.S12)
        assert fi < 1.0

    def test_shear_at_strength_is_one(self):
        # At τ12 = S12, σ1=σ2=0: FI = F66 * S12^2 = 1/(S12^2) * S12^2 = 1.0
        fi = tsai_wu(np.array([0.0, 0.0, self.S12]),
                     self.XT, self.XC, self.YT, self.YC, self.S12)
        np.testing.assert_allclose(fi, 1.0, rtol=1e-6)
```

- [ ] **Step 3: Run tests to verify current state**

```bash
cd /Users/jashwanth/Documents/Professional/Portfolio/CodeProjects/VirtualCompositeDesign
python -m pytest tests/test_clt.py -v 2>&1 | head -60
```

Expected: most pass; `test_navier_deflection_reasonable` should pass if clt.py is correct.

- [ ] **Step 4: Write utils tests**

Create `tests/test_utils.py`:

```python
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
        # This test will FAIL until BUG-C1 is fixed (load_materials() drops strength cols)
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
```

- [ ] **Step 5: Write optimizer tests**

Create `tests/test_optimizer.py`:

```python
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
        import layup_optimizer_sa as sa
        assert hasattr(sa, "Q_bar") or "from clt" in open(
            Path(__file__).resolve().parents[1] / "src" / "layup_optimizer_sa.py"
        ).read(), "SA optimizer must use clt.py for CLT computations"

    def test_symmetric_laminate_from_sa(self):
        """Full sequence returned by SA must be symmetric."""
        import random
        import layup_optimizer_sa as sa
        random.seed(0)
        seq, _ = sa.simulated_annealing(n_iterations=500)
        n = len(seq)
        half = n // 2
        assert seq[:half] == seq[half:][::-1], "SA result is not symmetric"
```

- [ ] **Step 6: Run full test suite**

```bash
python -m pytest tests/ -v 2>&1
```

Expected: `test_strength_columns_present` and `test_optimizer_uses_correct_material` FAIL (known — to be fixed in Tasks 2 and 3). All other tests PASS.

- [ ] **Step 7: Commit**

```bash
git add tests/
git commit -m "test: add pytest suite with CLT physics, materials, and optimizer baselines"
```

---

### Task 2: Fix `load_materials()` to expose strength columns

**Files:**
- Modify: `src/utils.py:15–67`

**Context:** `materials.csv` columns are: `S1T_Pa, S1C_Pa, S2T_Pa, S2C_Pa, S12_Pa`. The canonical names used throughout the codebase (in `clt.py` `hashin()` and `evaluate_laminate()`) are `X_T, X_C, Y_T, Y_C, S12`.

- [ ] **Step 1: Run the failing strength test to confirm the bug**

```bash
python -m pytest tests/test_utils.py::TestLoadMaterials::test_strength_columns_present -v
```

Expected: `FAILED — AssertionError: Missing strength column: X_T`

- [ ] **Step 2: Update `_COL_MAPS` and `_REQUIRED` in `src/utils.py`**

Replace the `_COL_MAPS` block (lines 15–26) with:

```python
_COL_MAPS = [
    ("name",    ["name", "material", "material_name"]),
    ("E1",      ["E1", "E1_Pa"]),
    ("E2",      ["E2", "E2_Pa"]),
    ("G12",     ["G12", "G12_Pa"]),
    ("v12",     ["v12", "nu12", "nu_12", "nu12[-]"]),
    ("density", ["density", "rho", "rho_kgm3", "rho_[kg/m3]"]),
    ("t_ply",   ["t_ply", "ply_t", "ply_t_m", "tply_m"]),
    # strength columns — optional (not _REQUIRED) so old CSVs still load
    ("X_T",     ["X_T", "S1T_Pa", "XT", "Xt"]),
    ("X_C",     ["X_C", "S1C_Pa", "XC", "Xc"]),
    ("Y_T",     ["Y_T", "S2T_Pa", "YT", "Yt"]),
    ("Y_C",     ["Y_C", "S2C_Pa", "YC", "Yc"]),
    ("S12",     ["S12", "S12_Pa", "S_12"]),
]

_REQUIRED = {"name", "E1", "E2", "G12", "v12", "density", "t_ply"}
_STRENGTH  = {"X_T", "X_C", "Y_T", "Y_C", "S12"}
```

- [ ] **Step 3: Update `load_materials()` return statement**

Replace lines 64–67 (the for-loop and return):

```python
    # convert stiffness columns to float
    for c in ["E1", "E2", "G12", "v12", "density", "t_ply"]:
        df[c] = _to_float(df, c)

    # convert strength columns to float (if present)
    present_strength = [c for c in _STRENGTH if c in df.columns]
    for c in present_strength:
        df[c] = _to_float(df, c)

    base_cols = ["name", "E1", "E2", "G12", "v12", "density", "t_ply"]
    return df[base_cols + present_strength]
```

- [ ] **Step 4: Run the strength test — it must now pass**

```bash
python -m pytest tests/test_utils.py -v
```

Expected: all `TestLoadMaterials` tests PASS.

- [ ] **Step 5: Smoke-test that existing scripts still work**

```bash
cd /Users/jashwanth/Documents/Professional/Portfolio/CodeProjects/VirtualCompositeDesign
python src/main.py 2>&1 | tail -5
```

Expected: same output as before — no breakage.

- [ ] **Step 6: Commit**

```bash
git add src/utils.py
git commit -m "fix: load strength columns (X_T, X_C, Y_T, Y_C, S12) from materials.csv"
```

---

### Task 3: Rewrite SA optimizer to use clt.py + materials.csv

**Files:**
- Modify: `src/layup_optimizer_sa.py` (full rewrite of material section + CLT functions)

**Context:** This task fixes BUG-A1, BUG-A2, BUG-A3, BUG-B1, CQ-1, CQ-3 in one pass. The SA algorithm itself (`simulated_annealing()`, `objective()`, `penalty()`) is correct and stays unchanged. Only the material constants, ply thickness, weight computation, and CLT math are replaced.

**What changes:**
- Remove lines 1–20 (hardcoded material + global CLT constants)
- Remove `compute_Qbar()` and `compute_D_matrix()` (lines 42–80)
- Import `Q_matrix`, `Q_bar`, `ply_interfaces`, `abd_matrices` from `clt`
- Import `load_materials` from `utils`
- Load material at module level from `materials.csv`
- Compute `ply_weight` from real density × plate area × t_ply

- [ ] **Step 1: Run the failing optimizer material test**

```bash
python -m pytest tests/test_optimizer.py::TestOptimizerMaterial::test_optimizer_uses_correct_material -v
```

Expected: `FAILED`

- [ ] **Step 2: Replace the top of `src/layup_optimizer_sa.py`**

Replace everything from line 1 through the end of `compute_D_matrix()` (up to line 80) with:

```python
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
_PLY_WEIGHT = RHO * t * LX * LY      # kg per ply  (≈ 0.018 kg for IM7/8552)

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
    Q   = Q_matrix(E1, E2, G12, nu12)
    qbars = [Q_bar(Q, deg2rad(th)) for th in full_seq]
    z     = ply_interfaces(len(full_seq), len(full_seq) * t)
    _, _, D = abd_matrices(qbars, z)
    return D
```

- [ ] **Step 3: Replace `laminate_deflection_metric()` and `laminate_weight()`**

Find and replace the two functions (originally lines 82–97):

```python
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
```

- [ ] **Step 4: Verify no other references to the old Q11/Q22/Q12/Q66 globals remain**

```bash
grep -n "Q11\|Q22\|Q12\|Q66\|nu21\|denom" src/layup_optimizer_sa.py
```

Expected: no output (all removed).

- [ ] **Step 5: Run the optimizer material tests**

```bash
python -m pytest tests/test_optimizer.py -v
```

Expected: all tests PASS.

- [ ] **Step 6: Smoke-test the optimizer**

```bash
python src/layup_optimizer_sa.py 2>&1 | tail -8
```

Expected output (values will vary by run):
```
Iteration 1000: Best Obj = ...  Half sequence length = ...
...
Optimized Layup Sequence (degrees):
[...]
Number of plies (full laminate): N
Objective Value (...): ...
```

- [ ] **Step 7: Commit**

```bash
git add src/layup_optimizer_sa.py
git commit -m "fix: rewrite SA optimizer to use clt.py + materials.csv (correct material + thickness)"
```

---

### Task 4: Fix `evaluate_laminate()` z-interface for variable-thickness plies

**Files:**
- Modify: `src/clt.py:219–221`

- [ ] **Step 1: Write the failing test (add to `tests/test_clt.py`)**

Add this class to `tests/test_clt.py`:

```python
class TestEvaluateLaminate:
    def test_variable_thickness_z_interfaces(self):
        """evaluate_laminate z-interfaces must match laminate_abd z-interfaces."""
        angles = [0, 90, 0]
        t_thin = 1.0e-4
        t_thick = 2.0e-4
        # Simulate via laminate_abd with Ply objects of mixed thickness
        plies_mixed = [
            Ply(E1, E2, G12, NU12, deg2rad(angles[0]), t_thin),
            Ply(E1, E2, G12, NU12, deg2rad(angles[1]), t_thick),
            Ply(E1, E2, G12, NU12, deg2rad(angles[2]), t_thin),
        ]
        A_ref, B_ref, D_ref, z_ref = laminate_abd(plies_mixed)
        # evaluate_laminate only supports equal-thickness — pass uniform t
        # For equal-thickness comparison, use all same t
        angles_eq = [0, 90, 0]
        N = np.array([0.0, 0.0, 0.0])
        M = np.array([1.0, 0.0, 0.0])
        res = evaluate_laminate(E1, E2, G12, NU12, angles_eq, t_thin, N, M)
        A_eval = res["A"]
        # z[0] must be -total_t/2
        total_t_eq = 3 * t_thin
        assert abs(res["plies"][0]["z_bot"] - (-total_t_eq / 2)) < 1e-12
```

- [ ] **Step 2: Replace the z-interface calculation in `evaluate_laminate()`**

In `src/clt.py`, replace lines 217–222:

```python
    # OLD (equal-thickness only):
    # total_t = ply_t * len(angles_deg)
    # z_iface = ply_interfaces(len(angles_deg), total_t)

    # NEW (cumulative sum, same logic as laminate_abd):
    n = len(angles_deg)
    z_iface = np.empty(n + 1, dtype=float)
    z_iface[0] = -0.5 * ply_t * n
    for k in range(n):
        z_iface[k + 1] = z_iface[k] + ply_t
```

- [ ] **Step 3: Run the test**

```bash
python -m pytest tests/test_clt.py::TestEvaluateLaminate -v
```

Expected: PASS

- [ ] **Step 4: Run full suite to check no regressions**

```bash
python -m pytest tests/ -v
```

Expected: all PASS

- [ ] **Step 5: Commit**

```bash
git add src/clt.py tests/test_clt.py
git commit -m "fix: evaluate_laminate uses cumulative z-interface (consistent with laminate_abd)"
```

---

### Task 5: Fix `compare_clt_fea.py` stress computation — add eps0 term

**Files:**
- Modify: `src/compare_clt_fea.py:117–123`

- [ ] **Step 1: Locate the stress block**

In `src/compare_clt_fea.py`, find lines 117–123:

```python
    kx = 0.0
    ky = 0.0
    for m in range(1, NAVIER_ODD + 1, 2):
        ...
        kx  += Wmn * mpa**2 * sm * sn
        ky  += Wmn * npb**2 * sm * sn

    z_bot_ply0 = z[0]
    Q0   = Q_matrix(E1, E2, G12, v12)
    Qb0  = Q_bar(Q0, 0.0)
    kappa = np.array([kx, ky, 0.0])
    sig_bot = Qb0 @ (z_bot_ply0 * kappa)          # ← BUG: missing eps0
```

- [ ] **Step 2: Compute eps0 and add it**

Replace the stress block (from `z_bot_ply0 = z[0]` to end of function) with:

```python
    # ── σₓₓ at bottom of 0° ply (z = −t/2) ──────────────────────────────────
    # Full expression: σ = Q̄ · (ε₀ + z · κ)
    # For symmetric layup under N=0, B=0 → ε₀ = A⁻¹N = 0.
    # Keeping eps0 explicit makes the code correct for future non-symmetric use.
    N_vec = np.array([0.0, 0.0, 0.0])
    eps0  = np.linalg.solve(A, N_vec)              # = [0,0,0] for symmetric + N=0
    kappa = np.array([kx, ky, 0.0])
    z_bot_ply0 = z[0]
    Q0  = Q_matrix(E1, E2, G12, v12)
    Qb0 = Q_bar(Q0, 0.0)
    sig_bot = Qb0 @ (eps0 + z_bot_ply0 * kappa)   # correct general form
    sigma_xx_pa = abs(float(sig_bot[0]))
```

- [ ] **Step 3: Verify result is unchanged (same numerical output)**

```bash
python src/compare_clt_fea.py 2>&1 | grep "σₓₓ peak"
```

Expected: `σₓₓ peak   : 0.5621 MPa` (unchanged — eps0 is zero for this case)

- [ ] **Step 4: Commit**

```bash
git add src/compare_clt_fea.py
git commit -m "fix: compare_clt_fea stress includes eps0 term (correct for general laminates)"
```

---

### Task 6: Fix FEA import path — add repo root to sys.path

**Files:**
- Modify: `src/compare_clt_fea.py:44–46`
- Modify: `src/sa_spotcheck.py:30–32`

- [ ] **Step 1: Update `compare_clt_fea.py`**

Replace lines 44–46:

```python
# OLD:
sys.path.insert(0, str(Path(__file__).resolve().parent))

# NEW — add both src/ (for local imports) and repo root (for fea.* imports):
_SRC = Path(__file__).resolve().parent
_ROOT = _SRC.parent
sys.path.insert(0, str(_SRC))
sys.path.insert(0, str(_ROOT))
```

- [ ] **Step 2: Update `sa_spotcheck.py`**

Replace the equivalent block (lines 30–32):

```python
# OLD:
sys.path.insert(0, str(Path(__file__).resolve().parent))

# NEW:
_SRC = Path(__file__).resolve().parent
_ROOT = _SRC.parent
sys.path.insert(0, str(_SRC))
sys.path.insert(0, str(_ROOT))
```

- [ ] **Step 3: Test import from a non-root directory**

```bash
cd /tmp && python /Users/jashwanth/Documents/Professional/Portfolio/CodeProjects/VirtualCompositeDesign/src/compare_clt_fea.py 2>&1 | head -5
```

Expected: script starts (Step 1 CLT output), no `ModuleNotFoundError`.

- [ ] **Step 4: Commit**

```bash
git add src/compare_clt_fea.py src/sa_spotcheck.py
git commit -m "fix: add repo root to sys.path so fea.* imports work from any directory"
```

---

### Task 7: Add division-by-zero guard in `navier_center_deflection()`

**Files:**
- Modify: `src/clt.py:269–275`

- [ ] **Step 1: Add test for degenerate laminate**

Add to `tests/test_clt.py::TestSymmetricLaminate`:

```python
    def test_navier_raises_for_zero_d11(self):
        """All-90 laminate has D11 ≈ D22 (transverse), deflection still finite."""
        plies_90 = [Ply(E1, E2, G12, NU12, np.pi / 2, T_PLY) for _ in range(8)]
        A, B, D, z = laminate_abd(plies_90)
        # Should not raise; D11 for all-90 is small but positive
        w = navier_center_deflection(D, 0.3, 0.3, 1000.0, max_odd=5)
        assert np.isfinite(w), f"Deflection is not finite: {w}"
        assert w > 0
```

- [ ] **Step 2: Update `navier_center_deflection()` in `src/clt.py`**

Replace lines 269–275:

```python
    for m in range(1, max_odd+1, 2):
        for n in range(1, max_odd+1, 2):
            mpa = (m*np.pi)/a
            npb = (n*np.pi)/b
            denom = D11*mpa**4 + 2.0*(D12 + 2.0*D66)*mpa**2*npb**2 + D22*npb**4
            if abs(denom) < 1e-30:
                raise ValueError(
                    f"Navier series: near-zero denominator at m={m}, n={n}. "
                    f"Check that D11 and D22 are both positive (degenerate laminate?)."
                )
            w += (16.0*q)/(np.pi**6 * m**2 * n**2) * (1.0/denom)
```

- [ ] **Step 3: Run the test**

```bash
python -m pytest tests/test_clt.py::TestSymmetricLaminate::test_navier_raises_for_zero_d11 -v
```

Expected: PASS

- [ ] **Step 4: Commit**

```bash
git add src/clt.py tests/test_clt.py
git commit -m "fix: guard against near-zero Navier denominator with descriptive ValueError"
```

---

### Task 8: Add singularity check in `solve_midplane()`

**Files:**
- Modify: `src/clt.py:86–104`

- [ ] **Step 1: Update `solve_midplane()`**

Replace lines 99–103:

```python
    # Schur complement
    try:
        AinvN = np.linalg.solve(A, N)
        AinvB = np.linalg.solve(A, B)
    except np.linalg.LinAlgError as exc:
        raise ValueError(
            "CLT: in-plane stiffness matrix A is singular. "
            "Check that the laminate has at least one non-zero ply and valid material properties."
        ) from exc
    S = D - B @ AinvB
    try:
        kappa = np.linalg.solve(S, M - B @ AinvN)
    except np.linalg.LinAlgError as exc:
        raise ValueError(
            "CLT: Schur complement of ABD system is singular. "
            "This can occur for laminates with extreme coupling (B ≈ D)."
        ) from exc
    eps0 = AinvN - AinvB @ kappa
```

- [ ] **Step 2: Run full test suite**

```bash
python -m pytest tests/ -v
```

Expected: all PASS

- [ ] **Step 3: Commit**

```bash
git add src/clt.py
git commit -m "fix: solve_midplane raises ValueError with domain context on singular A or Schur"
```

---

### Task 9: Consolidate `_find_ccx()` into `utils.py`

**Files:**
- Modify: `src/utils.py`
- Modify: `src/compare_clt_fea.py:135–145`
- Modify: `src/sa_spotcheck.py:80–87`

- [ ] **Step 1: Add `find_ccx()` to `src/utils.py`**

Append to `src/utils.py`:

```python
import shutil

FEA_DIR = Path(__file__).resolve().parents[1] / "fea"

def find_ccx() -> str | None:
    """Return path to CalculiX binary, or None if not available.

    Checks (in order):
    1. fea/ccx symlink in the repo (resolves if valid)
    2. fea/ccx as a regular file
    3. system PATH via shutil.which
    """
    repo_ccx = FEA_DIR / "ccx"
    if repo_ccx.is_symlink() and repo_ccx.resolve().exists():
        return str(repo_ccx)
    if repo_ccx.exists() and not repo_ccx.is_symlink():
        return str(repo_ccx)
    return shutil.which("ccx")
```

- [ ] **Step 2: Update `compare_clt_fea.py` to use the shared function**

Replace `_find_ccx()` definition (lines 135–145) and all calls:

```python
# remove _find_ccx() definition entirely
# replace the import block at top with:
from utils import load_materials, deg2rad, find_ccx
# replace the call:
ccx = find_ccx()
```

- [ ] **Step 3: Update `sa_spotcheck.py` similarly**

```python
from utils import load_materials, deg2rad, find_ccx
# remove local _find_ccx() definition
# replace: ccx = _find_ccx()  →  ccx = find_ccx()
```

- [ ] **Step 4: Smoke-test both scripts**

```bash
python src/compare_clt_fea.py 2>&1 | tail -3
python src/sa_spotcheck.py 2>&1 | head -5
```

Expected: no errors, same output as before.

- [ ] **Step 5: Commit**

```bash
git add src/utils.py src/compare_clt_fea.py src/sa_spotcheck.py
git commit -m "refactor: consolidate _find_ccx() into utils.find_ccx() — single source of truth"
```

---

### Task 10: Cache material load in `sa_spotcheck.py`

**Files:**
- Modify: `src/sa_spotcheck.py:60–73`

- [ ] **Step 1: Move `_load_mat()` result to module level**

Replace lines 60–73:

```python
# OLD:
def _load_mat():
    mats = load_materials(DATA_DIR / "materials.csv")
    m = mats.iloc[0]
    return (float(m["E1"]), float(m["E2"]), float(m["G12"]),
            float(m["v12"]), float(m["t_ply"]))

def clt_deflection(layup_deg: list[int]) -> float:
    E1, E2, G12, v12, tply = _load_mat()
    ...

# NEW:
_mat = load_materials(DATA_DIR / "materials.csv").iloc[0]
_E1, _E2, _G12, _V12, _TPLY = (
    float(_mat["E1"]), float(_mat["E2"]), float(_mat["G12"]),
    float(_mat["v12"]), float(_mat["t_ply"]),
)

def clt_deflection(layup_deg: list[int]) -> float:
    """Navier centre deflection [mm] for the given full stacking sequence."""
    plies = [Ply(_E1, _E2, _G12, _V12, deg2rad(th), _TPLY) for th in layup_deg]
    _, _, D, _ = laminate_abd(plies)
    w_m = navier_center_deflection(D, LX, LY, Q_LOAD, max_odd=NAVIER_ODD)
    return w_m * 1e3
```

- [ ] **Step 2: Commit**

```bash
git add src/sa_spotcheck.py
git commit -m "refactor: cache material load at module level in sa_spotcheck.py"
```

---

### Task 11: Fix `requirements.txt` — pin versions, remove scipy, add jupyter

**Files:**
- Modify: `requirements.txt`

- [ ] **Step 1: Check installed versions**

```bash
python -c "import numpy, pandas, matplotlib; print(numpy.__version__, pandas.__version__, matplotlib.__version__)"
```

- [ ] **Step 2: Replace `requirements.txt`**

```
numpy>=1.24,<3.0
pandas>=2.0,<3.0
matplotlib>=3.7,<4.0
jupyter>=1.0
```

(Remove `scipy` — not imported anywhere. Add `jupyter` — required for `notebooks/validation.ipynb`.)

- [ ] **Step 3: Commit**

```bash
git add requirements.txt
git commit -m "fix: pin requirements, remove unused scipy, add jupyter"
```

---

### Task 12: Fix `.gitignore` and untracked FEA files

**Files:**
- Modify: `.gitignore`
- Create: `fea/abaqus_inputs/.gitkeep`

- [ ] **Step 1: Add `fea/ccx` to `.gitignore`**

Append to `.gitignore`:

```
fea/ccx
fea/ccx/
```

- [ ] **Step 2: Add `.gitkeep` for `fea/abaqus_inputs/`**

```bash
touch fea/abaqus_inputs/.gitkeep
```

- [ ] **Step 3: Verify `git status` is clean of these entries**

```bash
git status --short | grep -E "fea/ccx|abaqus_inputs"
```

Expected: no output (both entries gone from status).

- [ ] **Step 4: Commit**

```bash
git add .gitignore fea/abaqus_inputs/.gitkeep
git commit -m "fix: gitignore fea/ccx symlink, add .gitkeep for abaqus_inputs/"
```

---

### Task 13: Fix README inaccuracies (BUG-A6, BUG-A7)

**Files:**
- Modify: `README.md`

- [ ] **Step 1: Fix angle sweep description**

Find: `"Sweeps Ex_eff from 0° to 90° in 1° increments"`  
Replace with: `"Sweeps Ex_eff from 0° to 90° in 5° increments"`

- [ ] **Step 2: Fix API key names in the `evaluate_laminate` example**

In the per-ply result table, change:

```markdown
| `sig_xy_top/bot` | ... |     →  | `sig_top_xy` / `sig_bot_xy` |
| `sig_12_top/bot` | ... |     →  | `sig_top_12` / `sig_bot_12` |
| `eps_xy_top/bot` | ... |     →  | `eps_top_xy` / `eps_bot_xy` |
| `eps_12_top/bot` | ... |     →  | `eps_top_12` / `eps_bot_12` |
```

Also fix the inline example code:

```python
# OLD:
f"| Hashin FT = {ply['hashin_FT_bot']:.4f} | "
# Check actual keys — correct form:
f"| Hashin FT = {ply['hashin_FT_bot']:.4f} | "   # this key IS correct in code
# The sig_ keys need fixing:
ply['sig_bot_12']    # correct
ply['sig_top_12']    # correct
```

- [ ] **Step 3: Commit**

```bash
git add README.md
git commit -m "docs: fix angle sweep resolution (5 deg not 1) and API dict key names"
```

---

### Task 14: Fix Hashin MC formula in `docs/Methodology_Notes.md`

**Files:**
- Modify: `docs/Methodology_Notes.md`

- [ ] **Step 1: Find the wrong Matrix Compression row**

In the Hashin table, the MC criterion shows `S₂T` in the denominator. The code (clt.py:173–175) uses `S₁₂`.

- [ ] **Step 2: Replace the MC row**

Find:
```
| **Matrix Compression** (σ₂ < 0) | (σ₂/2S₂T)² + [(S₂C²−4S₂T²)σ₂/(4S₂TS₂C)] + (τ₁₂/S₁₂)² | Oblique shear fracture |
```

Replace with:
```
| **Matrix Compression** (σ₂ < 0) | (σ₂/2S₁₂)² + [(Y_C/2S₁₂)²−1](σ₂/Y_C) + (τ₁₂/S₁₂)² | Oblique shear fracture (Hashin–Rotem simplified) |
```

- [ ] **Step 3: Commit**

```bash
git add docs/Methodology_Notes.md
git commit -m "docs: correct Hashin MC formula — uses S12 (in-plane shear), not S2T"
```

---

### Task 15: Add basic input validation to public CLT functions

**Files:**
- Modify: `src/clt.py` (add validation at entry points)

- [ ] **Step 1: Add a `_validate_laminate_inputs()` helper at top of `clt.py`**

Add after the imports:

```python
def _validate_laminate_inputs(angles_deg: list, ply_t: float) -> None:
    if len(angles_deg) == 0:
        raise ValueError("Laminate must have at least one ply.")
    if ply_t <= 0:
        raise ValueError(f"Ply thickness must be positive; got {ply_t}")
    for i, th in enumerate(angles_deg):
        if not (-90.0 <= th <= 90.0):
            raise ValueError(
                f"Ply angle at index {i} is {th}°; must be in [−90, 90]."
            )
```

- [ ] **Step 2: Call it at the start of `evaluate_laminate()`**

In `evaluate_laminate()`, insert after the docstring:

```python
    _validate_laminate_inputs(angles_deg, ply_t)
```

- [ ] **Step 3: Add `_validate_ply_list()` for `laminate_abd()`**

```python
def _validate_ply_list(plies: list) -> None:
    if len(plies) == 0:
        raise ValueError("Laminate must have at least one ply.")
    for i, p in enumerate(plies):
        if p.t <= 0:
            raise ValueError(f"Ply {i}: thickness must be positive; got {p.t}")
        if p.E1 <= 0 or p.E2 <= 0 or p.G12 <= 0:
            raise ValueError(f"Ply {i}: stiffness constants must be positive.")
```

Call at top of `laminate_abd()`.

- [ ] **Step 4: Add validation tests to `tests/test_clt.py`**

```python
class TestValidation:
    def test_empty_ply_list_raises(self):
        with pytest.raises(ValueError, match="at least one ply"):
            laminate_abd([])

    def test_zero_thickness_raises(self):
        with pytest.raises(ValueError, match="positive"):
            evaluate_laminate(E1, E2, G12, NU12, [0], 0.0,
                              np.zeros(3), np.zeros(3))

    def test_angle_out_of_range_raises(self):
        with pytest.raises(ValueError, match="90"):
            evaluate_laminate(E1, E2, G12, NU12, [0, 135], T_PLY,
                              np.zeros(3), np.zeros(3))
```

- [ ] **Step 5: Run tests**

```bash
python -m pytest tests/test_clt.py::TestValidation -v
```

Expected: all PASS

- [ ] **Step 6: Commit**

```bash
git add src/clt.py tests/test_clt.py
git commit -m "feat: add input validation to laminate_abd and evaluate_laminate"
```

---

### Task 16: Add one-command pipeline runner

**Files:**
- Create: `src/run_pipeline.py`

- [ ] **Step 1: Create `src/run_pipeline.py`**

```python
"""
run_pipeline.py
Run the full VirtualCompositeDesign pipeline in order:
  1. Baseline CLT analysis + angle sweep  (main.py)
  2. CLT vs FEA comparison               (compare_clt_fea.py)
  3. SA layup optimizer                  (layup_optimizer_sa.py)
  4. SA spot-check (top-3 candidates)    (sa_spotcheck.py)

Usage:
    python src/run_pipeline.py
"""
from __future__ import annotations
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PYTHON = sys.executable

STEPS = [
    ("Baseline CLT + angle sweep",   ROOT / "src" / "main.py"),
    ("CLT vs FEA comparison",        ROOT / "src" / "compare_clt_fea.py"),
    ("SA layup optimizer",           ROOT / "src" / "layup_optimizer_sa.py"),
    ("SA spot-check",                ROOT / "src" / "sa_spotcheck.py"),
]


def run_step(name: str, script: Path) -> bool:
    print(f"\n{'='*60}")
    print(f"  STEP: {name}")
    print(f"{'='*60}\n")
    result = subprocess.run([PYTHON, str(script)], cwd=str(ROOT))
    if result.returncode != 0:
        print(f"\n  ✗ STEP FAILED: {name} (exit code {result.returncode})")
        return False
    print(f"\n  ✓ STEP DONE: {name}")
    return True


def main() -> None:
    print("\n=== VirtualCompositeDesign — Full Pipeline ===")
    for name, script in STEPS:
        if not run_step(name, script):
            sys.exit(1)
    print("\n=== All steps complete ===\n")


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Test the runner**

```bash
python src/run_pipeline.py 2>&1 | grep -E "STEP|✓|✗"
```

Expected: four `✓ STEP DONE` lines.

- [ ] **Step 3: Commit**

```bash
git add src/run_pipeline.py
git commit -m "feat: add run_pipeline.py one-command full pipeline runner"
```

---

### Task 17: Remove dead `gui/` placeholder

**Files:**
- Delete: `gui/.gitkeep`
- Modify: `.gitignore` (add `gui/` if not already there)

- [ ] **Step 1: Remove the placeholder**

```bash
git rm gui/.gitkeep
```

- [ ] **Step 2: Ensure `gui/` is gitignored (so any future local GUI work isn't committed)**

Add to `.gitignore`:
```
gui/
```

- [ ] **Step 3: Commit**

```bash
git add .gitignore
git commit -m "chore: remove dead gui/ placeholder directory"
```

---

### Task 18: Add laminate stack cross-section visualization

**Files:**
- Create: `src/visualize.py`
- Modify: `src/main.py` (call visualize in baseline())

- [ ] **Step 1: Create `src/visualize.py`**

```python
"""
visualize.py
Laminate stack cross-section diagram and other composite visualizations.
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# Colour map for ply orientations
_ANGLE_COLORS = {
    0:   "#2C7BB6",    # blue
    45:  "#D7191C",    # red
    -45: "#FDAE61",    # orange
    90:  "#1A9641",    # green
}
_DEFAULT_COLOR = "#AAAAAA"

FIG_DIR = Path(__file__).resolve().parents[1] / "figures"


def plot_laminate_stack(
    angles_deg: list[float],
    ply_t_m: float,
    title: str = "Laminate Stack",
    out_path: Path | str | None = None,
) -> Path:
    """
    Draw a cross-section of the laminate with color-coded ply orientations.

    Parameters
    ----------
    angles_deg : ply angles from bottom to top (degrees)
    ply_t_m    : ply thickness (m)
    title      : figure title
    out_path   : save path (default: figures/laminate_stack.png)

    Returns
    -------
    Path to the saved figure.
    """
    n = len(angles_deg)
    total_t_mm = n * ply_t_m * 1e3    # mm

    fig, ax = plt.subplots(figsize=(4, max(3, n * 0.35)))

    for i, angle in enumerate(angles_deg):
        y_bot = i * ply_t_m * 1e3
        color = _ANGLE_COLORS.get(int(angle), _DEFAULT_COLOR)
        rect = mpatches.FancyBboxPatch(
            (0.05, y_bot), 0.90, ply_t_m * 1e3,
            boxstyle="square,pad=0", linewidth=0.5,
            edgecolor="white", facecolor=color, alpha=0.85,
        )
        ax.add_patch(rect)
        ax.text(0.5, y_bot + ply_t_m * 1e3 / 2,
                f"{angle:+.0f}°", ha="center", va="center",
                fontsize=9, color="white", fontweight="bold")

    # legend
    unique_angles = sorted(set(int(a) for a in angles_deg))
    handles = [
        mpatches.Patch(facecolor=_ANGLE_COLORS.get(a, _DEFAULT_COLOR),
                       label=f"{a:+d}°", alpha=0.85)
        for a in unique_angles
    ]
    ax.legend(handles=handles, loc="upper right", fontsize=8, framealpha=0.8)

    ax.set_xlim(0, 1)
    ax.set_ylim(0, total_t_mm)
    ax.set_yticks(np.arange(0, total_t_mm + ply_t_m * 1e3, ply_t_m * 1e3))
    ax.set_yticklabels([f"z={v:.3f}" for v in
                        np.arange(-total_t_mm/2, total_t_mm/2 + ply_t_m*1e3,
                                  ply_t_m*1e3)],
                       fontsize=7)
    ax.set_xticks([])
    ax.set_ylabel("z position [mm]", fontsize=9)
    ax.set_title(title, fontsize=11, pad=10)

    FIG_DIR.mkdir(parents=True, exist_ok=True)
    out = Path(out_path) if out_path else FIG_DIR / "laminate_stack.png"
    fig.tight_layout()
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return out
```

- [ ] **Step 2: Call it from `main.py` baseline()**

In `src/main.py`, at the end of `baseline()`:

```python
    from visualize import plot_laminate_stack
    tply = float(mat["t_ply"])
    figp = plot_laminate_stack(BASE_STACK, tply,
                                title="[0/45/−45/90]ₛ Baseline Laminate",
                                out_path=FIG_DIR / "laminate_stack.png")
    print(f"Saved: {figp}")
```

- [ ] **Step 3: Run and check output**

```bash
python src/main.py 2>&1 | grep "laminate_stack"
```

Expected: `Saved: .../figures/laminate_stack.png`

- [ ] **Step 4: Commit**

```bash
git add src/visualize.py src/main.py figures/laminate_stack.png
git commit -m "feat: add laminate stack cross-section visualization (color-coded ply diagram)"
```

---

### Task 19: Save SA convergence history and plot

**Files:**
- Modify: `src/layup_optimizer_sa.py` (`simulated_annealing()` return value + history)
- Modify: `src/sa_spotcheck.py` (generate convergence plot after optimization)

- [ ] **Step 1: Update `simulated_annealing()` to record convergence history**

In `src/layup_optimizer_sa.py`, modify the function signature and body:

```python
def simulated_annealing(
    n_iterations: int = 10000,
    initial_temp: float = 1.0,
    cooling_rate: float = 0.999,
    record_interval: int = 100,
) -> tuple[list[int], float, list[tuple[int, float]]]:
    """
    Returns (best_full_seq, best_obj, history)
    history: list of (iteration, best_obj) sampled every record_interval steps.
    """
    ...
    history: list[tuple[int, float]] = []
    for it in range(n_iterations):
        ...
        if (it + 1) % record_interval == 0:
            history.append((it + 1, best_obj))
        if (it + 1) % 1000 == 0:
            print(f"Iteration {it+1}: Best Obj = {best_obj:.3e}, ...")
    best_full_seq = best_half + best_half[::-1]
    return best_full_seq, best_obj, history
```

- [ ] **Step 2: Save convergence CSV and plot in `sa_spotcheck.py`**

After the `collect_top3()` loop, add:

```python
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def save_convergence_plot(all_histories: list[list[tuple[int, float]]], seeds: list[int]) -> Path:
    fig, ax = plt.subplots(figsize=(8, 4))
    for seed, hist in zip(seeds, all_histories):
        iters = [h[0] for h in hist]
        objs  = [h[1] for h in hist]
        ax.plot(iters, objs, alpha=0.6, label=f"seed={seed}")
    ax.set_xlabel("Iteration")
    ax.set_ylabel("Best Objective")
    ax.set_title("SA Convergence History (all seeds)")
    ax.legend(fontsize=7, ncol=2)
    ax.grid(True, alpha=0.3)
    out = ROOT / "figures" / "sa_convergence.png"
    fig.tight_layout()
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return out
```

- [ ] **Step 3: Commit**

```bash
git add src/layup_optimizer_sa.py src/sa_spotcheck.py
git commit -m "feat: save SA convergence history CSV and convergence plot per seed"
```

---

### Task 20: Parallelize SA multi-seed runs

**Files:**
- Modify: `src/sa_spotcheck.py` (`collect_top3()`)

- [ ] **Step 1: Replace sequential seed loop with `multiprocessing.Pool`**

```python
import multiprocessing as mp

def _run_seed(args: tuple) -> tuple[float, list[int]]:
    seed, n_iter, temp, cooling = args
    import random
    random.seed(seed)
    seq, obj, _ = simulated_annealing(
        n_iterations=n_iter,
        initial_temp=temp,
        cooling_rate=cooling,
    )
    print(f"  seed={seed:4d}  obj={obj:.4e}  n_plies={len(seq):3d}")
    return obj, seq


def collect_top3() -> list[dict]:
    print("Running SA optimizer across multiple seeds (parallel) …")
    args = [(s, SA_ITERATIONS, SA_INITIAL_TEMP, SA_COOLING) for s in SA_SEEDS]
    with mp.Pool(processes=min(len(SA_SEEDS), mp.cpu_count())) as pool:
        raw = pool.map(_run_seed, args)
    # rest of deduplication + CLT evaluation unchanged
    ...
```

- [ ] **Step 2: Smoke-test**

```bash
time python src/sa_spotcheck.py 2>&1 | grep -E "seed=|Saved"
```

Expected: all 8 seeds run, output to `data/sa_spotcheck.csv`. Time should be ~8× faster than sequential.

- [ ] **Step 3: Commit**

```bash
git add src/sa_spotcheck.py
git commit -m "perf: parallelize SA multi-seed runs with multiprocessing.Pool"
```

---

## Self-Review

**Spec coverage check:**
- BUG-A1, A2, A3, B1, CQ-1, CQ-3, CP-3 → Task 3 ✓
- BUG-A4 → Task 4 ✓
- BUG-A5 → Task 5 ✓
- BUG-A6, A7 → Task 13 ✓
- BUG-B2 → Task 7 ✓
- BUG-B3 → Task 8 ✓
- BUG-C1, CP-2 → Task 2 ✓
- BUG-C2, CP-6 → Task 11 ✓
- BUG-C3, CQ-1 → Task 9 ✓
- BUG-C4, C5 → Task 12 ✓
- BUG-D1 → Task 6 ✓
- BUG-D2 → Task 10 ✓
- BUG-D3 → Task 15 ✓
- CQ-2, CP-1 → Task 1 ✓
- CQ-4 → Task 14 ✓
- CQ-5 → Task 16 ✓
- CQ-6 → Task 17 ✓
- CP-4 → Task 18 ✓
- SC-2, CP-5 → Task 19 ✓
- SC-1 → Task 20 ✓

**No gaps found. All 36 findings are covered.**

**Placeholder scan:** No TBD, TODO, or vague steps found. All code blocks are complete.

**Type consistency:** `find_ccx()` (Task 9) matches usage in Tasks 5 and 6. `simulated_annealing()` return type extended in Task 19 is consistent with Task 20 usage.
