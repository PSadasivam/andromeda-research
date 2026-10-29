"""Plot Milky Way-Andromeda separation tracks for two runs side by side.

Example:
    python scripts/plot_tracks.py --out ../deep_space_portal/Images/andromeda_futures.png
"""

import argparse
import json
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from andromeda.integrator import TIME_UNIT_GYR

ROOT = Path(__file__).resolve().parents[1]
PANELS = (
    ("sawala2025_fiducial_n10000", "2025 inputs (Sawala et al.)", "#c0392b"),
    ("wu2026_fiducial_n10000", "2026 inputs (Wu et al.)", "#2563eb"),
)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--n", type=int, default=200, help="tracks per panel")
    parser.add_argument("--dpi", type=int, default=150)
    args = parser.parse_args()

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4), sharey=True, dpi=args.dpi)
    for ax, (run, title, colour) in zip(axes, PANELS):
        data = np.load(ROOT / "runs" / run / "ensemble.npz")
        summary = json.loads((ROOT / "runs" / run / "summary.json").read_text(encoding="utf-8"))["summary"]
        tracks = data["mw_m31_separation_track"][: args.n]
        dt_gyr = float(summary["horizon_gyr"]) / (tracks.shape[1] - 1)
        t = np.arange(tracks.shape[1]) * dt_gyr

        for track in tracks:
            ax.plot(t, track, color=colour, alpha=0.12, lw=0.8)
        ends = np.array([np.flatnonzero(np.isfinite(tr))[-1] for tr in tracks])
        merged = ends < tracks.shape[1] - 1
        ax.scatter(t[ends[merged]], np.full(merged.sum(), 8.0), marker="|", s=40, color="#111", lw=0.8)

        frac = summary["pairs"]["MW-M31"]["fraction"]
        ax.set_title(f"{title}\n{frac:.0%} merge within {summary['horizon_gyr']:.1f} billion years",
                     fontsize=11, color="#111")
        ax.set_xlabel("Billions of years from now")
        ax.set_xlim(0, t[-1])
        ax.grid(alpha=0.25)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)

    axes[0].set_ylabel("Milky Way - Andromeda distance (kpc)")
    axes[0].set_ylim(0, 900)
    fig.text(0.5, -0.02, f"Each line is one simulated future ({args.n} of 10,000 per panel). "
             "Tick marks at the bottom show when a merger occurs.", ha="center", fontsize=9, color="#555")
    fig.tight_layout()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, bbox_inches="tight", facecolor="white", dpi=args.dpi)
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
