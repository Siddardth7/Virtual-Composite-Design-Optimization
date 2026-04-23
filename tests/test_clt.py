import sys
from pathlib import Path
import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from clt import (
    Q_matrix, Q_bar, abd_matrices, ply_interfaces,
    laminate_abd, solve_midplane, navier_center_deflection,
    hashin, tsai_wu, Ply, deg2rad, evaluate_laminate,
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


class TestEvaluateLaminate:
    def test_variable_thickness_z_interfaces(self):
        """evaluate_laminate z-interfaces: z_bot[0] must equal -total_t/2."""
        angles = [0, 90, 0]
        N = np.array([0.0, 0.0, 0.0])
        M = np.array([1.0, 0.0, 0.0])
        res = evaluate_laminate(E1, E2, G12, NU12, angles, T_PLY, N, M)
        total_t = 3 * T_PLY
        assert abs(res["plies"][0]["z_bot"] - (-total_t / 2)) < 1e-15
