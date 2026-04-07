# Composite Laminate Analysis & Layup Optimization Toolkit

> **Python implementation of Classical Laminate Theory (CLT) with CalculiX FEA validation and Simulated Annealing layup optimization.**

A complete engineering toolkit for analysing, validating, and optimizing composite laminate structures.
Built from first principles to support the 2024–25 SAMPE Fuselage Competition at USC — then extended
into a fully validated, open-source composites analysis library.

**What you can do with this toolkit:**
- Compute full **[A], [B], [D] stiffness matrices** for any symmetric or unsymmetric laminate
- Recover **ply-level stresses and strains** in global and local material axes
- Evaluate **Tsai–Wu** and **Hashin 1980** failure indices at every ply interface
- Predict **SSSS plate centre deflection** via Navier double-sine series
- **Validate CLT against CalculiX FEA** with automated PASS/FAIL checkpoints
- **Optimise stacking sequences** for minimum weight + maximum bending stiffness using Simulated Annealing

---

## Table of Contents

1. [Theory Background](#theory-background)
2. [Repository Structure](#repository-structure)
3. [Installation](#installation)
4. [Usage](#usage)
   - [Baseline Analysis & Angle Sweep](#1-baseline-analysis--angle-sweep)
   - [CLT vs FEA Comparison](#2-clt-vs-fea-comparison)
   - [Simulated Annealing Optimizer](#3-simulated-annealing-optimizer)
   - [SA Spot-Check](#4-sa-spot-check)
   - [Validation Notebook](#5-validation-notebook)
5. [API Reference](#api-reference)
6. [Validation Results](#validation-results)
7. [SA Optimisation Results](#sa-optimisation-results)
8. [Material Properties](#material-properties-im78552)
9. [Problem Setup](#problem-setup)
10. [Dependencies](#dependencies)
11. [References](#references)

---

## Theory Background

### Classical Laminate Theory (CLT)

CLT models a composite laminate as a stack of thin, linearly elastic, orthotropic plies
bonded under **Kirchhoff–Love** plate kinematics (small deflections, perfect bonding, plane stress
per ply). The fundamental constitutive relation is the **ABD system**:

```
⌈ N ⌉   ⌈ A  B ⌉ ⌈ ε₀ ⌉
⌊ M ⌋ = ⌊ B  D ⌋ ⌊ κ  ⌋
```

| Matrix | Name | Role |
|--------|------|------|
| **A** (3×3) | In-plane stiffness | Relates in-plane force resultants N to mid-plane strains ε₀ |
| **B** (3×3) | Coupling stiffness | Couples bending to in-plane response (zero for symmetric laminates) |
| **D** (3×3) | Bending stiffness | Relates moment resultants M to curvatures κ |

Each ply contributes to A, B, D through its **transformed reduced stiffness Q̄(θ)** — a function
of the ply orientation angle θ, assembled by integrating over ply z-interfaces.

### Navier Series Centre Deflection

For a simply-supported (SSSS) rectangular plate under uniform pressure *q*, the centre deflection is:

```
w(Lx/2, Ly/2) = Σₘ Σₙ  Wmn · sin(mπ/2) · sin(nπ/2)

where:  Wmn = 16q / [π⁶ m² n² · (D₁₁α⁴ + 2(D₁₂+2D₆₆)α²β² + D₂₂β⁴)]
              α = mπ/Lx,  β = nπ/Ly
              sum over odd m, n = 1, 3, 5, …
```

Five terms per direction (m, n = 1, 3, 5, 7, 9) achieve < 0.1% convergence error for typical
composite plates, confirmed against a 9-term reference.

### Failure Criteria

**Tsai–Wu** — scalar quadratic polynomial in stress space.  
F₁₂ interaction term uses Tsai–Hahn: `F₁₂ = −0.5 √(F₁₁ · F₂₂)`.

**Hashin 1980** — four physically distinct failure modes in local (1-2) ply axes:

| Mode | Active when | Failure Index |
|------|-------------|---------------|
| Fibre Tension (FT) | σ₁ ≥ 0 | (σ₁/X_T)² + (τ₁₂/S₁₂)² |
| Fibre Compression (FC) | σ₁ < 0 | (σ₁/X_C)² |
| Matrix Tension (MT) | σ₂ ≥ 0 | (σ₂/Y_T)² + (τ₁₂/S₁₂)² |
| Matrix Compression (MC) | σ₂ < 0 | (σ₂/2S₁₂)² + [(Y_C/2S₁₂)²−1](σ₂/Y_C) + (τ₁₂/S₁₂)² |

FI ≥ 1 indicates initiation of that failure mode. The MC expression is the Hashin–Rotem simplified form.

---

## Repository Structure

```
composite-laminate-clt/
│
├── src/                            # All Python source modules
│   ├── clt.py                      # Core CLT engine (Q, Q̄, ABD, solve, recovery, Hashin, Navier)
│   ├── main.py                     # Entry-point: baseline laminate + parametric angle sweep
│   ├── utils.py                    # Material CSV loader + unit conversion helpers
│   ├── compare_clt_fea.py          # CLT vs CalculiX FEA comparison pipeline
│   ├── layup_optimizer_sa.py       # Simulated Annealing stacking sequence optimizer
│   └── sa_spotcheck.py             # SA top-3 spot-check with CLT/FEA ranking
│
├── fea/                            # Finite Element Analysis pipeline
│   ├── generate_inp.py             # Generates CalculiX .inp (mesh + material + BCs + step)
│   ├── parse_ccx_results.py        # Parses CalculiX .dat output (U3 deflection + S11 stress)
│   ├── abaqus_inputs/              # Runtime .inp / .dat / .frd files (gitignored)
│   └── results/
│       └── fea_summary.csv         # Pre-computed FEA reference results (baseline laminate)
│
├── data/
│   ├── materials.csv               # Material property database (IM7/8552)
│   ├── fea_comparison.csv          # CLT vs FEA comparison output table
│   ├── sa_spotcheck.csv            # SA top-3 candidate sequences with CLT deflections
│   └── sweeps/
│       ├── baseline_ply_summary.csv    # Per-ply stresses, strains, failure indices
│       └── angle_sweep_ex.csv          # Ex_eff vs ply orientation (0°→90°)
│
├── notebooks/
│   └── validation.ipynb            # 5-cell Jupyter validation notebook (Checkpoint 2)
│
├── docs/
│   ├── Clt_theory.pdf              # Full CLT derivation (theory reference)
│   └── Methodology_Notes.md        # Assumptions, CLT/Hashin theory, CalculiX setup, validation log
│
├── figures/                        # Plots generated by scripts and notebook
│   ├── angle_sweep_ex.png          # Ex_eff vs angle sweep
│   ├── fea_comparison.png          # CLT vs FEA bar chart (Checkpoint 1)
│   └── failure_indices.png         # Per-ply Hashin index bar chart (Checkpoint 2)
│
├── requirements.txt
└── .gitignore
```

---

## Installation

### Prerequisites

- **Python 3.9 or higher**
- *(Optional)* **CalculiX** — only needed for live FEA runs. Scripts fall back to pre-computed
  stored results automatically if `ccx` is not found.

### Step 1 — Clone the repository

```bash
git clone https://github.com/Siddardth7/Virtual-Composite-Design-Optimization.git
cd Virtual-Composite-Design-Optimization
```

### Step 2 — Create a virtual environment

```bash
python -m venv venv

# macOS / Linux
source venv/bin/activate

# Windows
venv\Scripts\activate
```

### Step 3 — Install Python dependencies

```bash
pip install -r requirements.txt
```

### Step 4 *(Optional)* — Install CalculiX for live FEA

```bash
# macOS via conda-forge (recommended)
conda install -c conda-forge calculix

# Ubuntu / Debian
sudo apt-get install calculix-ccx

# Verify installation
ccx          # should print usage / version info
```

> **Without CalculiX**, `compare_clt_fea.py` loads pre-computed results from
> `fea/results/fea_summary.csv`, and `sa_spotcheck.py` runs in CLT-only mode.
> All CLT analysis, failure criteria, SA optimizer, and angle sweep work fully offline.

---

## Usage

### 1. Baseline Analysis & Angle Sweep

**Script:** `src/main.py`

Runs the full CLT pipeline on the baseline `[0/45/−45/90]ₛ` laminate and performs
a parametric sweep of ply orientation from 0° to 90°.

```bash
python src/main.py
```

**What it does:**

| Step | Output |
|------|--------|
| Assembles [A], [B], [D] matrices for the baseline layup | Printed to terminal |
| Solves mid-plane strains and curvatures under unit moment Mx | Printed to terminal |
| Computes Navier SSSS centre deflection at q = 1 kPa | Printed to terminal |
| Exports per-ply stress/strain/failure index table | `data/sweeps/baseline_ply_summary.csv` |
| Sweeps Ex_eff from 0° to 90° in 1° increments | `data/sweeps/angle_sweep_ex.csv` |
| Saves angle sweep plot | `figures/angle_sweep_ex.png` |

**Expected terminal output:**

```
=== Baseline Laminate =========================================
Stack: [0, 45, -45, 90, 90, -45, 45, 0]
Total thickness t = 1.000 mm
||A|| = 1.045116e+08
||B|| = 3.087388e-13
||D|| = 1.029148e+01
kappa (unit Mx): [ 0.1170325  -0.0688027  -0.01706776]
SSSS Navier center deflection at q=1 kPa: 0.064411 mm
Saved: data/sweeps/baseline_ply_summary.csv
Saved: data/sweeps/angle_sweep_ex.csv
Saved: figures/angle_sweep_ex.png
```

> **Note:** `||B|| ≈ 3e-13` (essentially zero) confirms symmetry — the B matrix
> vanishes for symmetric laminates, as expected from theory.

---

### 2. CLT vs FEA Comparison

**Script:** `src/compare_clt_fea.py`

Validates the CLT Navier solution against CalculiX finite element analysis.
This is **Checkpoint 1** of the project validation plan.

```bash
python src/compare_clt_fea.py
```

**What it does:**

1. Computes CLT Navier centre deflection and peak σₓₓ for the baseline laminate
2. Loads FEA results (runs CalculiX live if `ccx` is available, otherwise loads `fea/results/fea_summary.csv`)
3. Computes % error for both metrics
4. Prints a formatted comparison table with PASS/FAIL status
5. Saves table to `data/fea_comparison.csv`
6. Saves side-by-side bar chart to `figures/fea_comparison.png`

**Expected terminal output:**

```
=== compare_clt_fea.py — Checkpoint 1 ===

── Step 1: CLT Navier analysis ─────────────────────────
  Deflection : 0.064411 mm
  σₓₓ peak   : 0.5621 MPa

── Step 2: FEA results ─────────────────────────────────
  Deflection : 0.063780 mm
  σₓₓ peak   : 0.5482 MPa
  Notes      : CalculiX 2.23 S4 20x20 mesh plate_ss.inp

── Step 3: Comparison table ────────────────────────────

───────────────────────────────────────────────────────────────────
Metric             Unit           CLT          FEA    Err %  Status
───────────────────────────────────────────────────────────────────
Deflection         mm        0.064411     0.063780   0.990%  ✅ PASS
σₓₓ (peak)         MPa       0.562070     0.548200   2.530%  ✅ PASS
───────────────────────────────────────────────────────────────────

══════════════════════════════════════════════════
  CHECKPOINT 1 REVIEW — D06
══════════════════════════════════════════════════
  Deflection         error = 0.990%  ✅ PASS
  σₓₓ (peak)         error = 2.530%  ✅ PASS
══════════════════════════════════════════════════
  ✅ CHECKPOINT 1 PASSED — both errors < 5%
══════════════════════════════════════════════════
```

**PASS threshold:** both metrics must be < 5% error. Both are well within bounds.

---

### 3. Simulated Annealing Optimizer

**Script:** `src/layup_optimizer_sa.py`

Optimises a symmetric laminate stacking sequence for **minimum weight + maximum bending stiffness**
using Simulated Annealing over the discrete orientation space {0°, ±45°, 90°}.

```bash
python src/layup_optimizer_sa.py
```

**Design variables:**
- Half-laminate ply sequence (mirrored to form the full symmetric laminate)
- Variable total ply count (4 to 40 plies)

**Constraints enforced via penalty:**
- Symmetry (full laminate = half + reversed half)
- Balance (equal ±45° ply count)
- Minimum 10% of each orientation (0°, 90°, ±45°) for laminates with ≥ 10 plies

**Objective:** `1/D₁₁ + weight_penalty`  — lower is better (stiffer, lighter)

**SA moves at each iteration:**
- `change` — swap one ply to a different orientation
- `insert` — add a ply at a random position (if below maximum)
- `delete` — remove a ply at a random position (if above minimum)

**Expected output (truncated):**
```
Iteration 1000: Best Obj = 2.622e-02, Half sequence length = 2
Iteration 2000: Best Obj = 2.622e-02, Half sequence length = 2
...

Optimized Layup Sequence (degrees):
[-45, 45, 45, -45]

Number of plies (full laminate): 4

Objective Value (deflection metric + weight + penalties): 0.026217
```

---

### 4. SA Spot-Check

**Script:** `src/sa_spotcheck.py`

Runs the SA optimizer across **8 different random seeds**, collects the top 3 unique
candidate layups, and ranks them by **CLT Navier centre deflection** (ascending = stiffer = better).
If CalculiX is available, also generates and runs FEA for each candidate to confirm CLT ranking.

```bash
python src/sa_spotcheck.py
```

**What it does:**

1. Runs SA (8 000 iterations) for each of 8 random seeds
2. Deduplicates candidates by exact stacking sequence
3. Computes CLT Navier deflection for each unique candidate using IM7/8552 material
4. Runs CalculiX FEA for each candidate if `ccx` is available
5. Prints ranked comparison table
6. Saves results to `data/sa_spotcheck.csv`

**Expected output:**

```
══════════════════════════════════════════════════════════════════════════
  SA Spot-Check — Top 3 Candidate Layups
══════════════════════════════════════════════════════════════════════════
  Rank  n    w_CLT [mm]    w_FEA [mm]   CLT→FEA err  Stacking Sequence
────────────────────────────────────────────────────────────────────────
  #1    24    0.002031          N/A    N/A (no ccx)   [0,-45,45,...,0]
  #2    20    0.003495          N/A    N/A (no ccx)   [0,-45,-45,...,0]
  #3    4     0.387517          N/A    N/A (no ccx)   [-45,45,45,-45]
══════════════════════════════════════════════════════════════════════════
  Saved: data/sa_spotcheck.csv
```

> **Without CalculiX**, CLT ranking is used directly. This is validated as reliable —
> Checkpoint 1 showed < 1% CLT/FEA deflection error for the same plate and material.

---

### 5. Validation Notebook

**File:** `notebooks/validation.ipynb`

A 5-cell Jupyter notebook providing an **interactive end-to-end validation** of the
full CLT-FEA pipeline, including ABD matrices, per-ply failure tables, and publication-quality plots.

```bash
jupyter notebook notebooks/validation.ipynb
```

Then select **Kernel → Restart & Run All** to execute all cells.

| Cell | Content |
|------|---------|
| **Cell 1** | ABD matrix assembly + per-ply Tsai–Wu and Hashin index table (all 4 modes) |
| **Cell 2** | CalculiX subprocess call (auto-falls back to stored results) |
| **Cell 3** | Inline CLT vs FEA comparison table |
| **Cell 4** | Annotated side-by-side bar chart → `figures/fea_comparison.png` |
| **Cell 5** | Per-ply Hashin failure index bar chart → `figures/failure_indices.png` |

---

## API Reference

### `src/clt.py` — Core CLT Engine

The core module. Import directly for custom analysis workflows.

```python
import sys
sys.path.insert(0, "src")
from clt import Ply, laminate_abd, navier_center_deflection, evaluate_laminate
from utils import load_materials, deg2rad
```

#### `Ply` dataclass

```python
@dataclass
class Ply:
    E1:    float   # Longitudinal modulus (Pa)
    E2:    float   # Transverse modulus (Pa)
    G12:   float   # In-plane shear modulus (Pa)
    v12:   float   # Major Poisson's ratio
    theta: float   # Ply orientation angle (radians)
    t:     float   # Ply thickness (m)
```

#### `laminate_abd(plies) → (A, B, D, z)`

Assembles the full ABD stiffness matrices and returns the ply z-interface array.

```python
plies = [Ply(161e9, 11.4e9, 5.17e9, 0.32, deg2rad(th), 1.25e-4)
         for th in [0, 45, -45, 90, 90, -45, 45, 0]]

A, B, D, z = laminate_abd(plies)
print(f"D11 = {D[0,0]:.4e} N·m")   # Primary bending stiffness
print(f"||B|| = {np.linalg.norm(B):.3e}")  # Should be ~0 for symmetric
```

#### `navier_center_deflection(D, Lx, Ly, q, max_odd=9) → float`

Computes SSSS plate centre deflection via Navier series.

```python
w_m = navier_center_deflection(D, Lx=0.3, Ly=0.3, q=1000.0, max_odd=5)
print(f"Centre deflection: {w_m * 1e3:.6f} mm")
# → Centre deflection: 0.064411 mm
```

#### `evaluate_laminate(E1, E2, G12, nu12, angles, t_ply, N, M, strengths=None) → dict`

High-level function: builds the laminate, solves, recovers per-ply stresses, and
(if `strengths` is provided) evaluates all failure indices.

```python
strengths = {
    "X_T": 2.8e9,   # Fibre tensile strength (Pa)
    "X_C": 1.6e9,   # Fibre compressive strength (Pa)
    "Y_T": 70e6,    # Matrix tensile strength (Pa)
    "Y_C": 200e6,   # Matrix compressive strength (Pa)
    "S12": 100e6,   # In-plane shear strength (Pa)
}

result = evaluate_laminate(
    E1=161e9, E2=11.4e9, G12=5.17e9, nu12=0.32,
    angles=[0, 45, -45, 90, 90, -45, 45, 0],
    t_ply=1.25e-4,
    N=np.array([0, 0, 0]),
    M=np.array([1, 0, 0]),   # unit Mx
    strengths=strengths,
)

# Access per-ply results
for i, ply in enumerate(result["plies"]):
    print(f"Ply {i+1:2d} ({ply['angle_deg']:4.0f}°): "
          f"σ₁ = {ply['sig_12_bot'][0]/1e6:8.3f} MPa | "
          f"Hashin FT = {ply['hashin_FT_bot']:.4f} | "
          f"Hashin MT = {ply['hashin_MT_bot']:.4f}")
```

**`result["plies"]` — keys available per ply:**

| Key | Description |
|-----|-------------|
| `angle_deg` | Ply orientation in degrees |
| `sig_xy_top/bot` | Global stress [σₓ, σᵧ, τₓᵧ] at top/bottom interface (Pa) |
| `sig_12_top/bot` | Local stress [σ₁, σ₂, τ₁₂] at top/bottom interface (Pa) |
| `eps_xy_top/bot` | Global strain at top/bottom interface |
| `eps_12_top/bot` | Local strain at top/bottom interface |
| `tsai_wu_top/bot` | Tsai–Wu failure index (FI ≥ 1 → failure) |
| `hashin_FT_top/bot` | Hashin fibre tension index |
| `hashin_FC_top/bot` | Hashin fibre compression index |
| `hashin_MT_top/bot` | Hashin matrix tension index |
| `hashin_MC_top/bot` | Hashin matrix compression index |

#### `hashin(stress_local, X_T, X_C, Y_T, Y_C, S12) → dict`

Low-level Hashin 1980 evaluation at a single point.

```python
from clt import hashin

fi = hashin(
    stress_local=np.array([150e6, 5e6, 10e6]),  # [σ₁, σ₂, τ₁₂] in Pa
    X_T=2.8e9, X_C=1.6e9,
    Y_T=70e6,  Y_C=200e6, S12=100e6,
)
# fi = {"FI_FT": 0.00289, "FI_FC": 0.0, "FI_MT": 0.00510, "FI_MC": 0.0}
```

### `src/utils.py` — Material Loader

```python
from utils import load_materials

mats = load_materials("data/materials.csv")
mat  = mats.iloc[0]   # IM7/8552

print(mat["E1"])      # 161000000000.0  (Pa)
print(mat["t_ply"])   # 0.000125        (m)
```

`load_materials` normalises common CSV header variants (`E1_Pa`, `rho_kgm3`, `ply_t_m`, etc.)
into canonical names (`E1`, `density`, `t_ply`) for consistent downstream access.

---

## Validation Results

### Checkpoint 1 — CLT vs FEA Baseline

**Laminate:** [0/45/−45/90]ₛ | **Material:** IM7/8552 | **Plate:** 300×300 mm SSSS | **Load:** 1 kPa

| Metric | CLT (Navier) | FEA (CalculiX S4, 20×20) | Error | Status |
|--------|-------------|--------------------------|-------|--------|
| Centre deflection | 0.064411 mm | 0.063780 mm | **0.99%** | ✅ PASS |
| Peak σₓₓ (0° ply bottom) | 0.5621 MPa | 0.5482 MPa | **2.53%** | ✅ PASS |

Both metrics below the 5% threshold. **Checkpoint 1 PASSED.**

> The FEA model uses **S4** (4-node bilinear shell) elements with 3-point Simpson
> through-thickness integration per ply, simply-supported boundary conditions on
> all four edges, and the same IM7/8552 engineering constants as the CLT model.

### Checkpoint 2 — Validation Notebook

All 5 cells execute cleanly with `Kernel → Restart & Run All`.

| Cell | Test | Result |
|------|------|--------|
| 1 | ABD matrix assembly + per-ply failure index table | ✅ |
| 2 | CalculiX subprocess (fallback to stored results) | ✅ |
| 3 | CLT vs FEA comparison table (both errors < 5%) | ✅ |
| 4 | Annotated comparison bar chart rendered | ✅ |
| 5 | Hashin failure index bar chart rendered | ✅ |

**Checkpoint 2 PASSED.**

---

## SA Optimisation Results

SA optimizer run across **8 random seeds** (8 000 iterations each), top 3 unique candidates
ranked by CLT Navier centre deflection (ascending = stiffer = better):

| CLT Rank | Plies | w_CLT [mm] | SA Objective | Stacking Sequence |
|----------|-------|-----------|-------------|-------------------|
| **#1** | 24 | **0.002031** | 0.1200 | `[0,−45,45,−45,45,−45,90,−45,45,−45,45,45,45,45,−45,45,−45,90,−45,45,−45,45,−45,0]` |
| **#2** | 20 | **0.003495** | 0.1000 | `[0,−45,−45,45,45,45,−45,45,−45,90,90,−45,45,−45,45,45,45,−45,−45,0]` |
| **#3** | 4  | **0.387517** | 0.0262 | `[−45,45,45,−45]` |

**Key insight:** The 4-ply `[±45]ₛ` candidate (#3) has the lowest SA objective but the
*highest* CLT deflection — it is the least stiff in bending. This occurs because the penalty
for minimum ply-fraction is zero at 4 plies (`required = int(4 × 0.1) = 0`), so the SA optimizer
converges to a pure ±45 solution that maximises in-plane shear stiffness rather than bending
stiffness D₁₁. The full CLT Navier evaluation (using IM7/8552 properties) correctly identifies
the 24-ply quasi-isotropic-rich sequence as the stiffest candidate.

CLT ranking validated as a reliable proxy for FEA ranking — Checkpoint 1 demonstrated
< 1% CLT/FEA deflection error for this plate geometry and material.

---

## Material Properties (IM7/8552)

Hexcel IM7 carbon fibre / 8552 epoxy — aerospace-grade unidirectional prepreg.
Source: `data/materials.csv`

| Property | Symbol | Value | Unit |
|----------|--------|-------|------|
| Longitudinal modulus | E₁ | 161 | GPa |
| Transverse modulus | E₂ | 11.4 | GPa |
| In-plane shear modulus | G₁₂ | 5.17 | GPa |
| Major Poisson's ratio | ν₁₂ | 0.32 | — |
| Density | ρ | 1600 | kg/m³ |
| Ply thickness | t_ply | 0.125 | mm |
| Fibre tensile strength | X_T | 2800 | MPa |
| Fibre compressive strength | X_C | 1600 | MPa |
| Matrix tensile strength | Y_T | 70 | MPa |
| Matrix compressive strength | Y_C | 200 | MPa |
| In-plane shear strength | S₁₂ | 100 | MPa |

To add a new material, append a row to `data/materials.csv` with the same column headers.
`utils.load_materials()` accepts common header name variants automatically.

---

## Problem Setup

All scripts share the same baseline problem definition:

| Parameter | Value | Notes |
|-----------|-------|-------|
| Plate dimensions | 300 × 300 mm | Square plate |
| Baseline layup | [0/+45/−45/90]ₛ | 8 plies, 1 mm total |
| Ply thickness | 0.125 mm | 8 × 0.125 = 1.000 mm |
| Boundary conditions | SSSS | Simply supported on all four edges |
| Applied load | 1 kPa uniform | Transverse pressure over full plate |
| FEA mesh | 20 × 20 S4 | 400 elements, 441 nodes |

The simply-supported BCs match the Navier series assumptions exactly:
- x-edges: `u = 0`, `w = 0` (in-plane + transverse displacement)
- y-edges: `v = 0`, `w = 0`
- Edge rotations left **free** → zero bending moments at boundaries

---

## Dependencies

```
numpy >= 1.24
scipy >= 1.10
matplotlib >= 3.7
pandas >= 2.0
```

Install all dependencies:
```bash
pip install -r requirements.txt
```

**Optional:**
- `jupyter` — for `notebooks/validation.ipynb`
- `CalculiX (ccx)` — for live FEA; install via `conda install -c conda-forge calculix`

---

## References

- Jones, R.M. (1999). *Mechanics of Composite Materials*, 2nd ed. Taylor & Francis.
- Daniel, I.M. & Ishai, O. (2006). *Engineering Mechanics of Composite Materials*. Oxford University Press.
- Tsai, S.W. & Wu, E.M. (1971). "A General Theory of Strength for Anisotropic Materials." *Journal of Composite Materials*, 5(1), 58–80.
- Hashin, Z. (1980). "Failure Criteria for Unidirectional Fiber Composites." *Journal of Applied Mechanics*, 47(2), 329–334.
- Reddy, J.N. (2004). *Mechanics of Laminated Composite Plates and Shells*, 2nd ed. CRC Press.
- MIL-HDBK-17-1F. *Polymer Matrix Composites — Guidelines for Characterization of Structural Materials.*
- ASTM D3039, D7264, D6641 — Standard test methods for composite mechanical properties.
