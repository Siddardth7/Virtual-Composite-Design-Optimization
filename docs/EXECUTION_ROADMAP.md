# Virtual Composite Design Optimization — Full Execution Roadmap
### Priority #2 Project — Completion Sprint Blueprint
**Project:** `VirtualCompositeDesign` (composite-laminate-clt)
**Resume Sprint Start:** March 27, 2026
**Target Launch:** April 9, 2026
**Owner:** Siddardth | M.S. Aerospace Engineering, UIUC
**Constraint:** 1 hour/day (5–6 PM), remote | Project is ~50% complete

---

## DAY 01 AUDIT — March 27, 2026 ✅ COMPLETE

### Code Review
| File | Status | Notes |
|------|--------|-------|
| `src/clt.py` | ✔ Verified | Q-matrix, ABD, Tsai-Wu, Navier all present. **Hashin not yet implemented** (block 6 ends at tsai_wu). evaluate_laminate() returns per-ply stresses/strains but no failure indices. Public API (laminate_abd, midplane_response, ply_strains_stresses) confirmed. |
| `src/main.py` | ✔ Verified | Builds [0/45/-45/90]s laminate, runs Navier deflection, angle sweep, saves CSV + PNG. Clean run confirmed. |
| `src/layup_optimizer_sa.py` | ✔ Verified | Standalone SA optimizer. Uses hardcoded material constants (not materials.csv). Optimizes symmetric laminates with balance + min-ply-% constraints. Note: uses its own Qbar/D matrix, NOT clt.py — separate implementation. |
| `src/utils.py` | ✔ Verified | Robust material loader with column aliasing. |

### Pipeline Run — `python src/main.py`
```
Stack: [0, 45, -45, 90, 90, -45, 45, 0]
Total thickness t = 1.000 mm
||A|| = 1.045116e+08  ||B|| = 3.087388e-13  ||D|| = 1.029148e+01
SSSS Navier center deflection at q=1 kPa: 0.064411 mm   ← BASELINE REFERENCE
```
Outputs saved: `data/sweeps/baseline_ply_summary.csv`, `data/sweeps/angle_sweep_ex.csv`, `figures/angle_sweep_ex.png`

### CalculiX Install
- **Platform:** Linux aarch64 (ARM) — pre-built x86 binary from dhondt.de is incompatible
- **Solution:** Installed via conda-forge (`miniconda` at `/sessions/vibrant-exciting-gates/miniconda/`)
- **Version:** CalculiX 2.23
- **ccx binary:** `/sessions/vibrant-exciting-gates/miniconda/bin/ccx`
- **Project symlink:** `fea/ccx` → points to the binary
- **Verify:** `./fea/ccx -v` → `This is Version 2.23` ✅

### Key Observations for Day 02
1. `src/clt.py` `evaluate_laminate()` does NOT yet call `tsai_wu()` inline — Hashin integration will be straightforward to add alongside it (Day 7–8).
2. `layup_optimizer_sa.py` is self-contained; does not use `clt.py`. For SA spot-check FEA (Day 11), the top-3 candidate layups will need to be fed into `evaluate_laminate()` from `clt.py`.
3. CLT baseline deflection = **0.064411 mm** at q=1 kPa for [0/45/-45/90]s, 300×300 mm plate. This is the reference FEA must match to within 5%.
4. The `fea/` directory is now created and ready for Day 02 work (`generate_inp.py`).

---

## TABLE OF CONTENTS

1. [Project Definition](#1-project-definition)
2. [Current State Audit](#2-current-state-audit)
3. [System Architecture](#3-system-architecture)
4. [Work Breakdown Structure (WBS)](#4-work-breakdown-structure-wbs)
5. [Master Timeline](#5-master-timeline)
6. [Weekly Execution Plan](#6-weekly-execution-plan)
7. [Daily Execution Plan](#7-daily-execution-plan-1-hourday)
8. [Task Tracking System](#8-task-tracking-system)
9. [Tools & Stack](#9-tools--stack)
10. [Deliverable Structure](#10-deliverable-structure-github)
11. [Quality Control System](#11-quality-control-system)
12. [Risks & Failure Points](#12-risks--failure-points)
13. [Final Launch Plan](#13-final-launch-plan)
14. [Execution Rules](#14-execution-rules)

---

## 1. PROJECT DEFINITION

### Technical Objective
Complete a Python-based computational framework for composite laminate analysis and optimal layup design. The core CLT engine and SA optimiser are fully built. The goal of this sprint is to close the only remaining gap: a working CalculiX FEA validation pipeline, Hashin failure criteria implementation, and a clean end-to-end validation notebook — making the project publishable.

### Scope

**REMAINING (this sprint):**
- CalculiX installation + configuration
- `fea/generate_inp.py` — CalculiX `.inp` file generator for composite shell plate
- `fea/parse_ccx_results.py` — parser for CalculiX `.dat` output (deflection + stress)
- `compare_clt_fea.py` — real CLT-vs-FEA comparison pipeline (not a stub)
- Hashin failure criteria in `clt.py` (4 failure modes)
- `notebooks/validation.ipynb` — complete 5-cell end-to-end validation notebook
- FEA spot-check of top 3 SA optimiser candidates
- README update with results + screenshots

**ALREADY COMPLETE (~50%):**
- `src/clt.py` — full CLT engine (Q-matrix, ABD, mid-plane solver, Tsai-Wu, Navier deflection)
- `src/main.py` — baseline analysis driver + parametric angle sweep
- `src/layup_optimizer_sa.py` — Simulated Annealing optimiser with symmetry/balance constraints
- `src/utils.py` + `data/materials.csv` — material loading utilities
- `docs/` — CLT theory PDF, LaTeX source, Methodology Notes
- `requirements.txt` — pinned dependencies

**EXCLUDED (out of scope):**
- GUI / web application
- Multi-material systems (IM7/8552 only)
- Dynamic / buckling analysis
- Any paid tool or paid cloud deployment

### Final Deliverables (What "Done" Looks Like)
1. `fea/generate_inp.py` — CalculiX input file generator
2. `fea/parse_ccx_results.py` — result parser (deflection + σₓₓ extraction)
3. `fea/baseline.inp` — generated CalculiX input for reference layup
4. `compare_clt_fea.py` — real comparison pipeline with bar chart output
5. `src/clt.py` — updated with Hashin failure criteria (4 modes)
6. `notebooks/validation.ipynb` — complete 5-cell validation notebook
7. `README.md` — updated with validation results, accuracy metrics, and usage
8. GitHub repository: `VirtualCompositeDesign` (public)

### Success Metrics
| Metric | Target |
|--------|--------|
| CLT vs FEA deflection error | < 5% on centre deflection |
| CLT vs FEA stress error | < 5% on σₓₓ at critical ply |
| Hashin failure modes | All 4 modes implemented (FT, FC, MT, MC) |
| SA candidate FEA validation | Top 3 candidates rank-ordered correctly vs CLT predictions |
| Notebook completeness | All 5 cells run clean, outputs visible |
| GitHub status | Public repo, README has accuracy table + usage instructions |

---

## 2. CURRENT STATE AUDIT

> **Day 01 is reserved for a complete review of existing code before any new work begins.**

### What Is Complete ✔

| File | Status | Notes |
|------|--------|-------|
| `src/clt.py` | ✔ Done | Q-matrix, ABD, Tsai-Wu, Navier — all functions implemented |
| `src/main.py` | ✔ Done | Baseline + angle sweep + CSV/PNG outputs |
| `src/layup_optimizer_sa.py` | ✔ Done | SA search with symmetry/balance constraints |
| `src/utils.py` | ✔ Done | Material loading from CSV |
| `data/materials.csv` | ✔ Done | Full IM7/8552 property database |
| `requirements.txt` | ✔ Done | Pinned NumPy, SciPy, Matplotlib, Pandas |
| `docs/` | ✔ Done | CLT theory PDF, LaTeX source, Methodology Notes |

### What Is Incomplete / Blocked ⚠️

| File | Status | Blocker |
|------|--------|---------|
| `fea/generate_inp.py` | ❌ Missing | Core FEA file to build — highest priority |
| `fea/parse_ccx_results.py` | ❌ Missing | Depends on generate_inp.py |
| `fea/baseline.inp` | ❌ Missing | Generated by generate_inp.py |
| `compare_clt_fea.py` | ⚠️ Stub | TODO placeholder only, no real pipeline |
| `src/clt.py` → Hashin | ⚠️ Missing | Hashin not implemented (Tsai-Wu only) |
| `notebooks/validation.ipynb` | ⚠️ Incomplete | CLT cells work; FEA cells are placeholders |

### The Abaqus → CalculiX Switch
- Abaqus is unavailable on macOS without a licence server
- **Solution: CalculiX** — free, open-source, Homebrew install, uses `.inp` syntax nearly identical to Abaqus, supports `*SHELL SECTION, COMPOSITE`
- Install: `brew install calculix` → verify with `ccx -v`

---

## 3. SYSTEM ARCHITECTURE

```
┌─────────────────────────────────────────────────────────────────┐
│                   VirtualCompositeDesign                         │
│                                                                  │
│  ┌──────────────────────┐    ┌──────────────────────────────┐   │
│  │   INPUT / MATERIAL   │    │        CLT ANALYSIS           │   │
│  │                      │    │                              │   │
│  │  materials.csv       │───▶│  clt.py                      │   │
│  │  (IM7/8552 props)    │    │   ├── Q_matrix()             │   │
│  │                      │    │   ├── ABD_matrix()           │   │
│  │  Layup definition    │    │   ├── mid_plane_solver()     │   │
│  │  [0/45/-45/90]s      │    │   ├── tsai_wu()             │   │
│  └──────────────────────┘    │   ├── hashin()   ← NEW       │   │
│                              │   ├── navier_deflection()   │   │
│                              │   └── evaluate_laminate()   │   │
│                              └──────────────┬───────────────┘   │
│                                             │                   │
│                              ┌──────────────▼───────────────┐   │
│                              │     SA OPTIMISER              │   │
│                              │                              │   │
│                              │  layup_optimizer_sa.py       │   │
│                              │   └── Candidate layups ×N   │   │
│                              └──────────────┬───────────────┘   │
│                                             │                   │
│  ┌──────────────────────┐    ┌──────────────▼───────────────┐   │
│  │   FEA LAYER (NEW)    │    │     COMPARISON PIPELINE       │   │
│  │                      │    │                              │   │
│  │  generate_inp.py     │───▶│  compare_clt_fea.py          │   │
│  │   └── baseline.inp   │    │   ├── evaluate_laminate()   │   │
│  │                      │    │   ├── run CalculiX (ccx)    │   │
│  │  ccx solver          │    │   ├── parse_ccx_results()   │   │
│  │   └── results.dat    │    │   ├── % error table         │   │
│  │                      │    │   └── bar chart (PNG)       │   │
│  │  parse_ccx_results   │    └──────────────┬───────────────┘   │
│  │   └── deflection     │                   │                   │
│  │   └── σₓₓ           │    ┌──────────────▼───────────────┐   │
│  └──────────────────────┘    │   VALIDATION NOTEBOOK         │   │
│                              │                              │   │
│                              │  validation.ipynb            │   │
│                              │   Cell 1: CLT baseline       │   │
│                              │   Cell 2: CalculiX run       │   │
│                              │   Cell 3: Parse + table      │   │
│                              │   Cell 4: Deflection plot    │   │
│                              │   Cell 5: Failure indices    │   │
│                              └──────────────────────────────┘   │
└─────────────────────────────────────────────────────────────────┘
```

### Layer Connections
- **CLT Layer → Comparison:** `evaluate_laminate()` produces ABD, strains, deflection, per-ply stresses
- **FEA Layer → Comparison:** `generate_inp.py` → `ccx` via subprocess → `parse_ccx_results.py` produces deflection + σₓₓ
- **Comparison → Notebook:** `compare_clt_fea.py` drives all notebook cells via importable functions
- **SA Optimiser → FEA spot-check:** top 3 candidate layups re-run through the FEA pipeline for cross-validation

---

## 4. WORK BREAKDOWN STRUCTURE (WBS)

```
PROJECT: VirtualCompositeDesign (Completion Sprint)
│
├── PHASE 1: FEA Setup (Days 1–3)
│   ├── MODULE 1.1: Environment + Code Review
│   │   ├── TASK 1.1.1: Review all existing code (clt.py, main.py, optimizer) — understand current API
│   │   ├── TASK 1.1.2: Install CalculiX via Homebrew + verify with `ccx -v`
│   │   └── TASK 1.1.3: Run `python src/main.py` — confirm existing pipeline still works
│   │
│   ├── MODULE 1.2: CalculiX Input Generator (Part 1)
│   │   ├── TASK 1.2.1: Write fea/generate_inp.py skeleton + argument interface
│   │   ├── TASK 1.2.2: Implement mesh generation (S8R elements, rectangular plate)
│   │   └── TASK 1.2.3: Implement *SHELL SECTION, COMPOSITE block from materials.csv values
│   │
│   └── MODULE 1.3: CalculiX Input Generator (Part 2)
│       ├── TASK 1.3.1: Implement simply-supported boundary conditions (4 edges)
│       ├── TASK 1.3.2: Implement uniform pressure load (matching CLT baseline)
│       ├── TASK 1.3.3: Implement *NODE PRINT + *EL PRINT output requests (U3, S)
│       └── TASK 1.3.4: Generate baseline.inp and verify CalculiX runs without errors
│
├── PHASE 2: FEA Parser + Comparison (Days 4–6)
│   ├── MODULE 2.1: Results Parser
│   │   ├── TASK 2.1.1: Write fea/parse_ccx_results.py — read CalculiX .dat file
│   │   ├── TASK 2.1.2: Extract centre node U3 deflection value
│   │   └── TASK 2.1.3: Extract σₓₓ at top and bottom surfaces of critical ply
│   │
│   ├── MODULE 2.2: Comparison Pipeline
│   │   ├── TASK 2.2.1: Replace TODO in compare_clt_fea.py — call evaluate_laminate()
│   │   ├── TASK 2.2.2: Add subprocess call to run CalculiX from compare_clt_fea.py
│   │   └── TASK 2.2.3: Compute % difference on deflection + σₓₓ, print comparison table
│   │
│   └── MODULE 2.3: Checkpoint 1 — FEA Pipeline Working
│       ├── TASK 2.3.1: Generate comparison bar chart (CLT vs FEA, deflection + σₓₓ)
│       ├── TASK 2.3.2: Save fea_comparison.csv and fea_comparison.png
│       └── TASK 2.3.3: Verify error < 5% on both metrics — document result
│
├── PHASE 3: Hashin Criteria (Days 7–8)
│   ├── MODULE 3.1: Hashin Implementation
│   │   ├── TASK 3.1.1: Add hashin(stress_local, strengths) to clt.py
│   │   ├── TASK 3.1.2: Implement Fibre Tension (FT) and Fibre Compression (FC) modes
│   │   ├── TASK 3.1.3: Implement Matrix Tension (MT) and Matrix Compression (MC) modes
│   │   └── TASK 3.1.4: Update evaluate_laminate() to return Hashin results alongside Tsai-Wu
│   │
│   └── MODULE 3.2: Integration Test
│       ├── TASK 3.2.1: Test hashin() against hand-calculated values for baseline layup
│       └── TASK 3.2.2: Confirm evaluate_laminate() output includes both criteria
│
├── PHASE 4: Validation Notebook (Days 9–10)
│   ├── MODULE 4.1: Notebook Cells 1–3
│   │   ├── TASK 4.1.1: Cell 1 — CLT baseline: display ABD matrix + per-ply stress table
│   │   ├── TASK 4.1.2: Cell 2 — CalculiX run: call compare_clt_fea.py via subprocess from notebook
│   │   └── TASK 4.1.3: Cell 3 — Parse + table: display FEA results + % error table
│   │
│   └── MODULE 4.2: Notebook Cells 4–5 + Checkpoint 2
│       ├── TASK 4.2.1: Cell 4 — Plot deflection + σₓₓ side-by-side with % error annotation
│       ├── TASK 4.2.2: Cell 5 — Plot Tsai-Wu vs Hashin failure indices per ply, side-by-side
│       └── TASK 4.2.3: Run full notebook top-to-bottom — confirm clean execution, no errors
│
├── PHASE 5: SA Spot-Check + Polish (Days 11–13)
│   ├── MODULE 5.1: FEA Spot-Check
│   │   ├── TASK 5.1.1: Run SA optimiser — save top 3 candidate layups
│   │   ├── TASK 5.1.2: Run each candidate through compare_clt_fea.py FEA pipeline
│   │   └── TASK 5.1.3: Confirm CLT ranking matches FEA ranking — document result
│   │
│   ├── MODULE 5.2: Documentation
│   │   ├── TASK 5.2.1: Update README.md — add accuracy table, usage instructions, validation results
│   │   ├── TASK 5.2.2: Update docs/Methodology_Notes — add CalculiX setup section + Hashin description
│   │   └── TASK 5.2.3: Finalize all TODO placeholders — no stubs remaining
│   │
│   └── MODULE 5.3: GitHub Launch
│       ├── TASK 5.3.1: Push all final code + notebook outputs to GitHub (public repo)
│       ├── TASK 5.3.2: Verify clean install in fresh venv (pip install -r requirements.txt)
│       └── TASK 5.3.3: Confirm notebook runs clean on fresh install
│
└── PHASE 6: Launch (Day 14)
    ├── TASK 6.1: Confirm GitHub repo is public
    ├── TASK 6.2: Capture 2–3 screenshots (notebook validation table, charts, failure index plot)
    ├── TASK 6.3: Draft LinkedIn post (validation results + methodology)
    └── TASK 6.4: Update resume Projects section with GitHub link + 1-line result summary
```

---

## 5. MASTER TIMELINE

**Total Duration:** ~2 weeks (14 days)
**Sprint Start:** March 27, 2026 (Friday) — Day 01
**Target Launch:** April 9, 2026 (Thursday)

```
PHASE 1  Mar 27–29   ▓▓▓  FEA SETUP: CalculiX install + generate_inp.py (Days 1–3)
PHASE 2  Mar 30–Apr 1 ▓▓▓  FEA PIPELINE: Parser + comparison script (Days 4–6)
PHASE 3  Apr 2–3     ▓▓   HASHIN: Failure criteria implementation (Days 7–8)
PHASE 4  Apr 4–5     ▓▓   NOTEBOOK: Complete validation.ipynb (Days 9–10)
PHASE 5  Apr 6–8     ▓▓▓  POLISH: SA spot-check + docs + GitHub (Days 11–13)
LAUNCH   Apr 9        ★   GITHUB PUBLIC + LINKEDIN POST
```

### Hard Deadlines

| Checkpoint | Date | What Must Be Done |
|-----------|------|-------------------|
| **Checkpoint 1** | April 1, 2026 | FEA pipeline running: generate_inp.py + ccx run + parse + comparison table with < 5% error |
| **Checkpoint 2** | April 5, 2026 | validation.ipynb runs clean top-to-bottom, all 5 cells producing output |
| **LAUNCH** | **April 9, 2026** | GitHub public, README with accuracy table, LinkedIn post published |

> **Policy:** If a checkpoint is missed, cut scope — not the deadline. Remove a feature before delaying launch.

---

## 6. WEEKLY EXECUTION PLAN

### WEEK 1 — March 27 – April 2 (Days 1–7)
**Objective:** FEA pipeline fully working. CalculiX installed, `.inp` file generated, solver runs, results parsed, CLT-vs-FEA comparison table produced with < 5% error.

**Deliverables to Complete:**
- CalculiX installed and verified
- `fea/generate_inp.py` — generates a valid `.inp` for the reference [0/45/-45/90]s layup
- `fea/baseline.inp` — generated file, committed to repo
- `fea/parse_ccx_results.py` — extracts U3 centre deflection + σₓₓ from `.dat` output
- `compare_clt_fea.py` — real pipeline (no TODO stub), outputs comparison table + bar chart

**Key Risks:**
- CalculiX `.inp` syntax errors → shell crashes silently. **Mitigation:** Start with simplest possible mesh (4-element plate), verify it runs before adding complexity.
- CalculiX `.dat` parsing edge cases. **Mitigation:** Print raw `.dat` output on Day 4 before writing parser.

**Expected Output:** `python compare_clt_fea.py` prints CLT vs FEA table and saves `fea_comparison.png`.

---

### WEEK 2 — April 3–9 (Days 8–14)
**Objective:** Hashin criteria integrated, validation notebook complete, SA spot-check done, project launched publicly.

**Deliverables to Complete:**
- `src/clt.py` — `hashin()` function with all 4 failure modes, `evaluate_laminate()` updated
- `notebooks/validation.ipynb` — all 5 cells run clean, outputs saved in notebook
- SA spot-check: top 3 candidates run through FEA, ranking verified
- `README.md` — updated with validation results table + accuracy summary
- GitHub: public repo, clean commit history
- LinkedIn post published

**Key Risks:**
- Hashin matrix compression mode is algebraically complex (Mohr-Coulomb friction term). **Mitigation:** Use simplified Hashin-Rotem form if full Hashin proves too complex in 1 day.
- Notebook subprocess calls to CalculiX may fail in Jupyter kernel. **Mitigation:** Test subprocess call in plain `.py` script first before embedding in notebook.

**Expected Output:** `jupyter nbconvert --to notebook --execute validation.ipynb` completes without errors.

---

## 7. DAILY EXECUTION PLAN (1 Hour/Day — 5–6 PM)

> **Format:** Day | Date | Task | Definition of Done

---

### PHASE 1: FEA Setup

| Day | Date | Task | Done When |
|-----|------|------|-----------|
| **D01** | Fri Mar 27 | **REVIEW SESSION.** Open and read `clt.py`, `main.py`, `layup_optimizer_sa.py`. Run `python src/main.py` to confirm existing pipeline works. Install CalculiX: `brew install calculix`, verify with `ccx -v`. Document current working state in a brief notes block at top of this roadmap. | `ccx -v` returns version string. `python src/main.py` produces output CSV/PNG without errors. |
| **D02** | Sat Mar 28 | Write `fea/generate_inp.py` Part 1: argument interface (layup, ply thickness, plate dims, load). Implement node + element mesh for a 4×4 grid rectangular plate using S8R shell elements. Write `*NODE` and `*ELEMENT` blocks to file. | Running `python fea/generate_inp.py` produces a `.inp` file with valid `*NODE` and `*ELEMENT` sections (verify by eye or with `grep '*NODE'`). |
| **D03** | Sun Mar 29 | Complete `fea/generate_inp.py` Part 2: add `*SHELL SECTION, COMPOSITE` block using IM7/8552 values from `materials.csv`. Add simply-supported BCs (`*BOUNDARY`), uniform pressure (`*DLOAD`), and output requests (`*NODE PRINT, U` + `*EL PRINT, S`). Run `ccx baseline` — verify solver exits cleanly (non-zero exit = error). | `ccx fea/baseline` runs without error, produces `baseline.dat`. File `fea/baseline.inp` committed. |
| **D04** | Mon Mar 30 | Write `fea/parse_ccx_results.py`. Open `baseline.dat`, read raw content. Extract: (1) U3 value at the centre node (highest magnitude), (2) σₓₓ (S11) at top and bottom of critical ply. Print both to stdout with units. | `python fea/parse_ccx_results.py baseline.dat` prints numerical values for deflection + σₓₓ. Values are non-zero and in physically reasonable range. |
| **D05** | Tue Mar 31 | Replace TODO in `compare_clt_fea.py`: (1) call `evaluate_laminate()` for baseline layup, (2) call `generate_inp.py` → run `ccx` via `subprocess.run()`, (3) call `parse_ccx_results.py`. Print raw CLT and FEA values side by side. | `python compare_clt_fea.py` runs end-to-end, prints CLT deflection and FEA deflection side by side without crashing. |
| **D06** | Wed Apr 1 | Add comparison metrics to `compare_clt_fea.py`: compute % error on deflection + σₓₓ, print formatted comparison table, save `fea_comparison.csv`. Generate side-by-side bar chart with % error annotation, save as `fea_comparison.png`. **Checkpoint 1:** Verify < 5% error on both metrics. | `python compare_clt_fea.py` outputs table + PNG. `fea_comparison.csv` committed. % error on deflection < 5%. |

---

### PHASE 2 + 3: Hashin Criteria

| Day | Date | Task | Done When |
|-----|------|------|-----------|
| **D07** | Thu Apr 2 | Implement `hashin(stress_local, strengths)` in `clt.py`. Add fibre failure modes: Fibre Tension (FT: σ₁ > 0, uses X_T + S₁₂) and Fibre Compression (FC: σ₁ < 0, uses X_C). Return `FI_FT` and `FI_FC` as floats. | Function exists, runs without error on dummy stress values. FI_FT > 1 for stress above fibre tension strength. |
| **D08** | Fri Apr 3 | Add matrix failure modes to `hashin()`: Matrix Tension (MT: σ₂ > 0, uses Y_T + S₁₂) and Matrix Compression (MC: σ₂ < 0, uses Hashin-Rotem simplified form with Y_C + S₁₂). Update `evaluate_laminate()` to include all 4 Hashin indices alongside Tsai-Wu. Test on baseline layup — hand-verify 1 ply. | `evaluate_laminate()` output dict contains `hashin_FT`, `hashin_FC`, `hashin_MT`, `hashin_MC` per ply. Spot-check one ply against manual calculation. |

---

### PHASE 4: Validation Notebook

| Day | Date | Task | Done When |
|-----|------|------|-----------|
| **D09** | Sat Apr 4 | Open `notebooks/validation.ipynb`. Write Cell 1: call `evaluate_laminate()`, display ABD matrix as DataFrame, display per-ply table (strains, Tsai-Wu, Hashin indices). Write Cell 2: call CalculiX via subprocess (same as `compare_clt_fea.py` flow). Write Cell 3: call `parse_ccx_results.py`, display FEA values + % error comparison table. | Cells 1–3 execute without error in Jupyter. Per-ply table and comparison table display inline. |
| **D10** | Sun Apr 5 | Write Cell 4: side-by-side bar chart of CLT vs FEA deflection + σₓₓ with % error annotation (Matplotlib). Write Cell 5: grouped bar chart of Tsai-Wu vs Hashin failure indices per ply. Run full notebook: `Kernel → Restart & Run All`. **Checkpoint 2.** | Full notebook runs top-to-bottom with no errors. All outputs (tables + charts) visible inline. Notebook committed with outputs. |

---

### PHASE 5: SA Spot-Check + Polish

| Day | Date | Task | Done When |
|-----|------|------|-----------|
| **D11** | Mon Apr 6 | Run SA optimiser (`python src/layup_optimizer_sa.py`). Save top 3 candidate layups. Run each through `compare_clt_fea.py` FEA pipeline. Print CLT rank vs FEA rank for all 3 candidates. | Top 3 candidates all run through FEA without errors. CLT and FEA rankings printed side by side. |
| **D12** | Tue Apr 7 | Verify SA candidate FEA results: confirm CLT ranking order matches FEA ranking (±1 position = acceptable). Document discrepancies if any. Update `docs/Methodology_Notes.md` — add CalculiX setup section, Hashin criteria description, SA spot-check results. | Methodology notes updated. If CLT rank ≠ FEA rank, document reason in notes. All stubs/TODOs removed from all files. |
| **D13** | Wed Apr 8 | Update `README.md`: (1) add CLT vs FEA accuracy table (deflection error, σₓₓ error), (2) add CalculiX installation instructions, (3) add usage section (`python compare_clt_fea.py`, `jupyter notebook validation.ipynb`), (4) add results summary. Test clean install in fresh venv. Tag commit `v1.0-complete`. | README shows accuracy table. Fresh venv install works. `v1.0-complete` tag pushed. |
| **D14** | Thu Apr 9 | **LAUNCH DAY.** (1) Confirm GitHub repo is public. (2) Run validation notebook one final time — all outputs visible. (3) Capture 2–3 screenshots (per-ply table, comparison bar chart, failure index plot). (4) Publish LinkedIn post. (5) Update resume Projects section with GitHub link + "CLT vs FEA error < 5%" as result bullet. | LinkedIn post live. GitHub repo public. Resume updated. |

---

## 8. TASK TRACKING SYSTEM

### Apple Reminders (Primary System)

**List Name:** `CLT Composite`

**Format for each reminder:**
```
Title: [D01] Review + CalculiX install
Date: Mar 27, 2026
Alert: 5:00 PM
Notes: Done when: ccx -v works, main.py runs clean
```

**All 14 daily reminders:**

| Reminder Title | Due Date | Time |
|----------------|----------|------|
| [D01] Code review + CalculiX install | Mar 27 | 5:00 PM |
| [D02] generate_inp.py: mesh + elements | Mar 28 | 5:00 PM |
| [D03] generate_inp.py: BCs + load + run ccx | Mar 29 | 5:00 PM |
| [D04] parse_ccx_results.py: deflection + σₓₓ | Mar 30 | 5:00 PM |
| [D05] compare_clt_fea.py: pipeline integration | Mar 31 | 5:00 PM |
| [D06] compare_clt_fea.py: table + chart ✅ Checkpoint 1 | Apr 1 | 5:00 PM |
| [D07] Hashin: fibre failure modes (FT + FC) | Apr 2 | 5:00 PM |
| [D08] Hashin: matrix modes (MT + MC) + evaluate_laminate update | Apr 3 | 5:00 PM |
| [D09] validation.ipynb: Cells 1–3 (CLT + FEA + table) | Apr 4 | 5:00 PM |
| [D10] validation.ipynb: Cells 4–5 + full run ✅ Checkpoint 2 | Apr 5 | 5:00 PM |
| [D11] SA spot-check: top 3 candidates through FEA | Apr 6 | 5:00 PM |
| [D12] SA verify + methodology notes update | Apr 7 | 5:00 PM |
| [D13] README update + clean install + tag v1.0 | Apr 8 | 5:00 PM |
| 🚀 [LAUNCH] GitHub public + LinkedIn + resume | Apr 9 | 5:00 PM |

---

## 9. TOOLS & STACK

| Tool | Role | Status |
|------|------|--------|
| Python 3.11 | Core language | Installed |
| NumPy + SciPy | Matrix algebra, linear solve | Installed |
| Matplotlib + Pandas | Plotting, data handling | Installed |
| CalculiX (`ccx`) | Free FEA solver (Abaqus replacement) | **Install Day 01** |
| Jupyter | Validation notebook | Installed |
| GitHub | Version control + public portfolio | Active |

---

## 10. DELIVERABLE STRUCTURE (GitHub)

```
VirtualCompositeDesign/
├── src/
│   ├── clt.py                 ✔ done  (+hashin → Day 7–8)
│   ├── main.py                ✔ done
│   ├── utils.py               ✔ done
│   └── layup_optimizer_sa.py  ✔ done
├── fea/                       ← BUILD THIS (Days 2–4)
│   ├── generate_inp.py        ❌ to create
│   ├── parse_ccx_results.py   ❌ to create
│   └── baseline.inp           ❌ to generate
├── compare_clt_fea.py         ⚠️ stub → Day 5–6
├── notebooks/
│   └── validation.ipynb       ⚠️ incomplete → Day 9–10
├── data/
│   ├── materials.csv          ✔ done
│   └── fea_comparison.csv     ❌ to generate
├── figures/                   ← OUTPUT (auto-generated)
│   └── fea_comparison.png
├── docs/
│   ├── Methodology_Notes.md   ✔ done (+CalculiX section → Day 12)
│   └── [CLT theory docs]      ✔ done
├── requirements.txt           ✔ done
└── README.md                  ⚠️ update → Day 13
```

---

## 11. QUALITY CONTROL SYSTEM

### Per-Day Check Before Ending Session
Before ending each 5–6 PM session, verify:
- [ ] Code runs without error (`python <file>` or notebook cell)
- [ ] Output is physically reasonable (non-zero values, correct sign on deflection)
- [ ] Changes committed to Git with a meaningful message
- [ ] Roadmap updated: current day marked complete

### Checkpoint Reviews
**Checkpoint 1 (Apr 1):** Run `python compare_clt_fea.py`. Confirm table printed, PNG saved, deflection error < 5%.

**Checkpoint 2 (Apr 5):** Run `Kernel → Restart & Run All` on `validation.ipynb`. Confirm all cells complete, all charts visible.

**Launch QC (Apr 8):** Fresh venv install, run full pipeline, README passes 30-second recruiter test.

---

## 12. RISKS & FAILURE POINTS

| Risk | Probability | Impact | Mitigation |
|------|-------------|--------|------------|
| CalculiX `.inp` syntax errors | High | High | Start with a 4-element plate (minimal mesh). Check CalculiX examples online. Verify before scaling. |
| `.dat` output format unexpected | Medium | Medium | Print raw `.dat` content first (Day 4), then write parser to match actual format. |
| Notebook subprocess to ccx fails in Jupyter | Medium | Medium | Test subprocess call as standalone `.py` script before embedding in notebook. |
| Hashin matrix compression algebraic complexity | Medium | Low | Use Hashin-Rotem simplified form (omits friction angle term). Document simplification. |
| CLT vs FEA error > 5% | Low | Medium | Check BCs in `.inp` — simply-supported requires 4-edge constraint, not pinned corners only. Verify load magnitude matches CLT baseline. |
| Git history messiness | Low | Low | One commit per day, imperative message format: "Add generate_inp.py mesh section" |

---

## 13. FINAL LAUNCH PLAN

### LinkedIn Post Template
```
Just completed a Python-based composite laminate analysis pipeline —
Classical Laminate Theory validated against CalculiX FEA with < 5%
error on centre deflection and σₓₓ.

The pipeline:
• CLT engine: ABD matrix, mid-plane solver, Tsai-Wu + Hashin failure criteria
• Simulated Annealing optimizer: finds optimal [0/±45/90] stacking sequences
• CalculiX FEA validation: automated .inp generation → solver run → result parse
• Jupyter notebook: end-to-end validation with visual comparison charts

Built with Python / NumPy / CalculiX — no proprietary tools required.
GitHub: [link]

#CompositesMaterials #StructuralAnalysis #AerospaceEngineering #Python #FEA
```

### Resume Bullet
```
• Built Python CLT + FEA validation pipeline for composite laminate optimization;
  CalculiX FEA cross-validated CLT predictions within 5% error on deflection and
  ply-level stress — Simulated Annealing optimiser confirmed via FEA spot-check
```

---

## 14. EXECUTION RULES

1. **One file at a time.** Don't split attention across multiple files per session. One session = one deliverable.
2. **Run before you commit.** Every commit must be a working state. No committing broken code.
3. **Physical sanity check every day.** Deflection should be downward (negative U3). σₓₓ should be tensile on bottom ply for simply-supported plate under pressure.
4. **CalculiX first, pretty code second.** Get the solver to run and produce output before refactoring anything.
5. **Document assumptions inline.** If you pick a mesh density or threshold, add a `# NOTE:` comment explaining why.
6. **50% done means 50% to go.** Treat this as a fresh sprint — don't assume old code is bug-free until you've run it.
7. **5–6 PM only.** Hard stop at 6 PM. If a task isn't done, continue the next day — don't bleed into FMEA time.
