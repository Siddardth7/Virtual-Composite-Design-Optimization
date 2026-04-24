"""
visualize.py
Laminate stack cross-section diagram and other composite visualizations.
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

_ANGLE_COLORS = {
    0:   "#2C7BB6",
    45:  "#D7191C",
    -45: "#FDAE61",
    90:  "#1A9641",
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
    total_t_mm = n * ply_t_m * 1e3

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
