# VirtualCompositeDesign — Full Project Audit Report

**Date:** 2026-04-23  
**Auditor:** Claude Code  
**Branch:** main @ 4f8ed8f  
**Goal:** Complete, correct, production-quality composite laminate analysis toolkit.

---

## Executive Summary

- **The CLT core engine (`clt.py`) and FEA pipeline (`compare_clt_fea.py`, `fea/`) are physically correct** — the Navier deflection, ABD assembly, Hashin 1980, and Tsai–Wu implementations have been validated against CalculiX FEA with < 1% deflection error and < 3% stress error.
- **The SA optimizer (`layup_optimizer_sa.py`) is completely disconnected from the validated CLT engine.** It uses hardcoded material constants from a different fiber system (T300/Epoxy instead of IM7/8552), a ply thickness 8× too large (1 mm vs 0.125 mm), and duplicates CLT math independently. Its optimization results are physically meaningless relative to the CLT validation pipeline.
- **The material loader (`utils.py`) silently discards all five strength columns** (S1T, S1C, S2T, S2C, S12) present in `materials.csv`, forcing every downstream script to hardcode strength values manually — making the CSV data partially useless.
- **There is no automated test suite.** The validation notebook is a manual check only. No script catches regressions.
- **The project is structurally sound** — repo layout, docstrings, FEA integration, and README quality are all strong. The gap between "mostly working" and "complete" is about 15 targeted fixes.

---

## Section 1 — Correctness Bugs (Wrong Physics / Wrong Math)

### BUG-A1 — Wrong material constants in SA optimizer ⚠️ CRITICAL

| | |
|---|---|
| **File** | `src/layup_optimizer_sa.py`, lines 9–15 |
| **Impact** | Critical — optimizer produces meaningless results |

**What it does:**
The optimizer hardcodes `E1=89e9, E2=8e9, G12=4.5e9, nu12=0.3` — properties consistent with T300/Epoxy or similar generic carbon/epoxy, not the project's IM7/8552.

`materials.csv` specifies IM7/8552: `E1=161e9, E2=11.4e9, G12=5.17e9, nu12=0.32`.

**Why it matters:**
The SA optimizer computes its `D11` stiffness metric using the wrong fiber system. A layup optimized for T300/Epoxy is not optimal for IM7/8552 — the anisotropy ratios (E1/E2 ≈ 11 for T300 vs ≈ 14 for IM7) produce different optimal ply-angle distributions. The CLT validation in `sa_spotcheck.py` then evaluates those candidate layups using the correct IM7/8552 properties — making the SA objective and the CLT ranking fundamentally incoherent.

---

### BUG-A2 — Wrong ply thickness in SA optimizer ⚠️ CRITICAL

| | |
|---|---|
| **File** | `src/layup_optimizer_sa.py`, line 24 |
| **Impact** | Critical — bending stiffness off by up to 512× |

**What it does:**
`t = 0.001` (1 mm per ply) is hardcoded. `materials.csv` defines `t_ply = 0.000125` m (0.125 mm per ply).

**Why it matters:**
Bending stiffness `D ∝ t³`. A ply that is 8× too thick produces a D matrix that is up to `8³ = 512×` too stiff. The optimizer's objective `1/D11` is therefore 512× smaller than it should be for any given layup, making the SA objective scores completely non-representative of the actual structural performance.

---

### BUG-A3 — Magic ply weight constant

| | |
|---|---|
| **File** | `src/layup_optimizer_sa.py`, lines 95–97 |
| **Impact** | High — weight objective is physically groundless |

**What it does:**
`ply_weight = 0.005` kg is hardcoded with no derivation.

**Correct value** for IM7/8552, 300×300 mm plate:
`mass_per_ply = density × t_ply × Lx × Ly = 1600 × 0.000125 × 0.3 × 0.3 = 0.018 kg`

The hardcoded value is 3.6× too small and ignores plate geometry entirely.

---

### BUG-A4 — `evaluate_laminate()` assumes equal-thickness plies

| | |
|---|---|
| **File** | `src/clt.py`, lines 219–221 |
| **Impact** | Medium — silent failure for mixed-thickness laminates |

**What it does:**
```python
total_t = ply_t * len(angles_deg)          # assumes all plies are ply_t thick
z_iface = ply_interfaces(len(angles_deg), total_t)  # divides total_t equally
```
`ply_interfaces()` uses `np.linspace` which only works for equal-thickness plies.

**Contrast with `laminate_abd()`** (lines 285–289):
```python
z[0] = -0.5 * total_t
for k in range(n):
    z[k+1] = z[k] + plies[k].t   # cumulative, handles variable thickness
```

If a user calls `evaluate_laminate()` with a mixed-thickness stack, the z-interface positions are wrong, which produces wrong ply strains and stresses without any error or warning.

---

### BUG-A5 — `compare_clt_fea.py` stress computation omits membrane term

| | |
|---|---|
| **File** | `src/compare_clt_fea.py`, lines 121–123 |
| **Impact** | Medium — silently wrong for non-symmetric or loaded laminates |

**What it does:**
```python
kappa = np.array([kx, ky, 0.0])
sig_bot = Qb0 @ (z_bot_ply0 * kappa)   # ← missing eps0
```

The correct expression is `sig_bot = Qb0 @ (eps0 + z_bot * kappa)`.

**Why it is currently not wrong:** For this specific case (symmetric layup, `B=0`, no in-plane loads `N=0`), `eps0 = A⁻¹ N = 0`. The omitted term happens to be zero. But the code contains an implicit assumption that is not documented and would silently fail for:
- Asymmetric laminates (B ≠ 0 → coupling creates nonzero eps0 under pure bending)
- Any applied in-plane load Nx, Ny, or Nxy

---

### BUG-A6 — README documents wrong angle sweep resolution

| | |
|---|---|
| **File** | `README.md`, Usage section step 1 |
| **Impact** | Low — documentation misleads users |

README states: *"Sweeps Ex_eff from 0° to 90° in **1° increments**"*

Code (`src/main.py`, line 117): `thetas = np.arange(0, 91, 5)` → **5° steps**, producing 19 data points, not 91.

---

### BUG-A7 — README API reference uses wrong dict key names

| | |
|---|---|
| **File** | `README.md`, API Reference section |
| **Impact** | High — users get `KeyError` when copying the example |

README example shows:
```python
ply['sig_12_bot']    # ← does not exist
ply['eps_12_top']    # ← does not exist
```

Actual keys in `evaluate_laminate()` output:
```python
ply['sig_bot_12']    # ← correct
ply['sig_top_12']    # ← correct
ply['eps_bot_12']    # ← correct
ply['eps_top_12']    # ← correct
```

---

## Section 2 — Numerical / Solver Bugs

### BUG-B1 — SA optimizer duplicates CLT math with wrong constants

| | |
|---|---|
| **File** | `src/layup_optimizer_sa.py`, lines 42–80 |
| **Impact** | Critical — optimizer stiffness metric is not from the validated engine |

`compute_Qbar()` (lines 42–64) and `compute_D_matrix()` (lines 66–80) re-implement the same math already present and validated in `clt.py`'s `Q_bar()` and `abd_matrices()`. The duplicate code:
- Uses the wrong material (BUG-A1)
- Uses the wrong ply thickness (BUG-A2)
- Diverges from the validated implementation with no test coverage

The optimizer's D matrix has no relationship to the CLT engine that validates the results.

---

### BUG-B2 — No division-by-zero guard in `navier_center_deflection()`

| | |
|---|---|
| **File** | `src/clt.py`, line 273 |
| **Impact** | Medium — crash or silent `inf` for degenerate laminates |

```python
denom = D11*mpa**4 + 2.0*(D12 + 2.0*D66)*mpa**2*npb**2 + D22*npb**4
w += (16.0*q)/(np.pi**6 * m**2 * n**2) * (1.0/denom)   # ← no guard
```

For an all-90° laminate, `D11 ≈ 0` and `denom → 0`. This produces `ZeroDivisionError` or `inf` with no actionable message.

---

### BUG-B3 — `solve_midplane()` has no singularity check

| | |
|---|---|
| **File** | `src/clt.py`, lines 99–103 |
| **Impact** | Medium — cryptic `LinAlgError` for degenerate inputs |

```python
AinvN = np.linalg.solve(A, N)    # raises LinAlgError if A is singular
AinvB = np.linalg.solve(A, B)
```

A singular `A` matrix (e.g., zero-stiffness layup, or a layup accidentally passed all-zero) raises `numpy.linalg.LinAlgError: Singular matrix` with no domain context for the user.

---

## Section 3 — I/O and Data Bugs

### BUG-C1 — `load_materials()` silently discards strength columns ⚠️ HIGH

| | |
|---|---|
| **File** | `src/utils.py`, lines 26, 64–67 |
| **Impact** | High — strength data in CSV is never used automatically |

`materials.csv` contains five strength columns: `S1T_Pa, S1C_Pa, S2T_Pa, S2C_Pa, S12_Pa`.

`load_materials()` defines `_REQUIRED` and the return statement to include only:
`["name", "E1", "E2", "G12", "v12", "density", "t_ply"]`

All strength data is silently dropped. Every script that needs failure criteria must hardcode strength values manually (e.g., `sa_spotcheck.py`, `notebooks/validation.ipynb`). Adding a new material row to the CSV has no effect on failure calculations.

---

### BUG-C2 — `requirements.txt` has no version pinning and a dead dependency

| | |
|---|---|
| **File** | `requirements.txt` |
| **Impact** | Medium — reproducibility risk; dead dep wastes install time |

```
numpy       # no version
scipy       # never imported in any source file
pandas      # no version
matplotlib  # no version
            # jupyter missing — needed for validation notebook
```

---

### BUG-C3 — `_find_ccx()` is duplicated with inconsistent logic

| | |
|---|---|
| **Files** | `src/compare_clt_fea.py` lines 138–141; `src/sa_spotcheck.py` lines 81–86 |
| **Impact** | Low — maintenance risk |

Two copies of the same function with reversed check order. Should be a single shared utility in `utils.py`.

---

### BUG-C4 — `fea/ccx` broken symlink is untracked and not gitignored

| | |
|---|---|
| **File** | `.gitignore` |
| **Impact** | Low — noisy `git status` on all clones |

`fea/ccx` appears as `?? fea/ccx` in `git status`. It is a broken symlink (created during development) that should be gitignored.

---

### BUG-C5 — `fea/abaqus_inputs/` directory appears as untracked

| | |
|---|---|
| **File** | `.gitignore`, `fea/abaqus_inputs/` |
| **Impact** | Low — noisy `git status` |

The directory appears as `?? fea/abaqus_inputs/`. The `.gitignore` covers `fea/abaqus_inputs/*.inp` but the directory itself needs a `.gitkeep` committed so the directory exists in the repo (consistent with `fea/results/.gitkeep`).

---

## Section 4 — Interface / Integration Bugs

### BUG-D1 — `compare_clt_fea.py` FEA import fails outside repo root

| | |
|---|---|
| **File** | `src/compare_clt_fea.py`, line 153 |
| **Impact** | Medium — breaks live FEA path for any non-root invocation |

```python
sys.path.insert(0, str(Path(__file__).resolve().parent))  # adds src/
# later:
from fea.generate_inp import generate_full_inp   # needs repo root in sys.path
```

Running `python src/compare_clt_fea.py` from `src/` or any subdirectory raises `ModuleNotFoundError: No module named 'fea'`. Fix: also insert the repo root.

---

### BUG-D2 — `sa_spotcheck.py` re-reads CSV on every deflection call

| | |
|---|---|
| **File** | `src/sa_spotcheck.py`, lines 60–64 |
| **Impact** | Trivial — minor inefficiency |

`_load_mat()` opens and parses `materials.csv` once per `clt_deflection()` call. Should be cached at module level.

---

### BUG-D3 — No input validation anywhere

| | |
|---|---|
| **Files** | All `src/` scripts |
| **Impact** | Medium — cryptic errors for bad inputs |

No script validates:
- Negative or zero ply thickness
- Empty ply list
- Ply angles outside [−90, 90]°
- Missing required CSV columns

Errors surface as bare `numpy` stack traces with no domain context.

---

## Section 5 — Code Quality Issues

| ID | File | Issue |
|----|------|-------|
| CQ-1 | `src/layup_optimizer_sa.py` | Does not import `clt.py` or `utils.py`; duplicates CLT math with wrong constants |
| CQ-2 | (absent) | No automated test suite — zero `pytest` files |
| CQ-3 | `src/layup_optimizer_sa.py` lines 16–20 | Global-scope material computation at import time; untestable in isolation |
| CQ-4 | `docs/Methodology_Notes.md` | Hashin MC formula shows wrong denominator (`S₂T` instead of `S₁₂`); code is correct, doc is wrong |
| CQ-5 | (absent) | No one-command pipeline runner — must run 4 scripts manually in order |
| CQ-6 | `gui/.gitkeep` | Dead placeholder directory; no GUI code exists or is planned |

---

## Section 6 — Scaling Opportunities

| ID | Area | Description |
|----|------|-------------|
| SC-1 | Performance | SA 8-seed runs are sequential; `multiprocessing.Pool` gives ~8× speedup |
| SC-2 | Observability | No SA convergence history saved; adding a convergence CSV/plot validates optimizer behaviour |
| SC-3 | Multi-material | Once BUG-C1 fixed, `sa_spotcheck.py` can sweep across multiple materials in the CSV |
| SC-4 | Architecture | FEM backend (CalculiX) is not abstracted; a thin adapter would make the solver swappable |
| SC-5 | Optimization | Single scalar objective (1/D₁₁ + weight); Pareto front (weight vs. deflection) is more physically meaningful |
| SC-6 | Robustness | No uncertainty quantification; Monte Carlo on material scatter (E₁ ± 5%) would show sensitivity |

---

## Section 7 — Completeness Gaps

| ID | Gap | Why It Blocks "Complete" |
|----|-----|--------------------------|
| CP-1 | No automated test suite | Cannot call an engineering library complete without tests against known analytical solutions |
| CP-2 | Strength data unused from CSV | `materials.csv` has strength columns that are silently ignored; failure calculations hardcode values |
| CP-3 | SA optimizer disconnected from CLT engine | Core pipeline has a broken link between optimization and validation |
| CP-4 | No laminate stack visualization | Standard output for any composites design tool; missing from all scripts |
| CP-5 | No SA convergence plot | No way to verify optimizer convergence without visual history |
| CP-6 | `requirements.txt` incomplete | `jupyter` missing; `scipy` listed but unused; no version pins |

---

## Finding Count Summary

| Category | Count | Critical | High | Medium | Low |
|----------|-------|----------|------|--------|-----|
| Correctness Bugs | 7 | 2 | 2 | 2 | 1 |
| Numerical/Solver | 3 | 1 | 0 | 2 | 0 |
| I/O and Data | 5 | 0 | 1 | 1 | 3 |
| Interface/Integration | 3 | 0 | 0 | 2 | 1 |
| Code Quality | 6 | — | — | — | — |
| Scaling | 6 | — | — | — | — |
| Completeness | 6 | — | — | — | — |
| **Total** | **36** | **3** | **3** | **7** | **5** |
