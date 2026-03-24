# Methodology Notes
## Theory
- Classical Laminate Theory (CLT): ABD matrices, mid-plane strains \(\{\varepsilon_0\}\), curvatures \(\{\kappa\}\).
- Ply-level stress recovery via transformed reduced stiffness \(\bar{\mathbf{Q}}\).
- Failure indices: Tsai–Wu or Hashin (implementation hooks provided).

## Assumptions
- Linearly elastic orthotropic plies.
- Perfect bonding, small strains, Kirchhoff–Love kinematics for thin laminates.
- Loads within elastic regime.

## Validation Plan
- Baseline laminate with known properties.
- Mesh convergence to <3% change in deflection.
- Compare CLT vs FEA: deflection, \(\sigma_{xx}\) at critical ply, failure index.
