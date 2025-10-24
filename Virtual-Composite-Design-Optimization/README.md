# Virtual Composite Structure Design & Optimization
## Overview
Parametric analysis of composite laminates using Classical Laminate Theory (CLT) with validation against FEA. Study sensitivity of deflection, stress, and failure indices to ply orientation and thickness.

## Objectives
- Implement CLT to compute [A], [B], [D], mid-plane strains/curvatures, and ply stresses.
- Validate analytical predictions against an FEA baseline.
- Run a small sweep over ply angles (0°, ±15°, ±30°, ±45°, 90°) and total thickness.
- Package results and figures into a short technical report.

## Methodology
Analytical module in Python (NumPy/SciPy) implements standard laminate theory for symmetric and unsymmetric stacks. FEA uses shell elements with matched geometry, BCs, and loads. Comparison metrics: tip deflection, peak ply stress, and failure index.

## Repository Structure
- /docs — methodology notes and report templates
- /src — CLT and comparison scripts
- /fea — solver inputs and results exports
- /data — materials and sweep CSVs
- /notebooks — validation notebooks
- /figures — plots saved by scripts

## How to Run
1. Create and activate a Python 3.10+ environment.
2. Install dependencies:
   ```bash
   pip install -r requirements.txt
   ```
3. Run a quick baseline:
   ```bash
   python src/main.py
   ```
4. Outputs:
   - CSVs in /data/sweeps/
   - Plots in /figures/
   - Console summary with CLT vs. baseline checks.

## Results & Validation
Start with a [0/90]s or [0/±45/90]s plate in bending. Target CLT–FEA deflection agreement ±8% after mesh convergence (<3% change between last two meshes).

## References
- Jones, R.M. (1999). *Mechanics of Composite Materials*.
- Daniel, I.M., & Ishai, O. (2006). *Engineering Mechanics of Composite Materials*.
- MIL-HDBK-17; ASTM D3039, D7264.
