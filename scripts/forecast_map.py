"""Pre-DR4 forecast map: where a new M31 proper motion would raise, lower, or leave the forecast unchanged.

Reads a uniform-prior scan (config wu2026_m31pm_scan), writes forecast.json beside it, and plots
the pre-stated regions for two precision scenarios.

Example:
    python scripts/run_ensemble.py --config wu2026_m31pm_scan --n 100000 --workers 14 --name wu2026_m31pm_scan_n100000
    python scripts/forecast_map.py --out ../deep_space_portal/Images/andromeda_forecast_map.png
"""

import argparse
import json
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Ellipse, Patch

from andromeda import forecast as fc
from andromeda.configs import CONFIGS

ROOT = Path(__file__).resolve().parents[1]
CATEGORY_CODES = {"lower": 0, "unchanged": 1, "raise": 2}
CATEGORY_COLOURS = ("#f2b8ab", "#e3e3e3", "#a9c8f5")
GRID_STEP = 0.5


def _scan(run):
    data = np.load(ROOT / "runs" / run / "ensemble.npz")
    meta = json.loads((ROOT / "runs" / run / "summary.json").read_text(encoding="utf-8"))["metadata"]
    spec = meta["parameters"]["M31"]
    box = ((spec["pmra"][1], spec["pmra"][2]), (spec["pmdec"][1], spec["pmdec"][2]))
    merged = np.isfinite(data["merge_time_MW-M31"])
    return data["obs_M31_pmra"], data["obs_M31_pmdec"], merged, box, meta


def _range_within(frac, ra, dec, centre, sigma, n_sigma):
    rr, dd = np.meshgrid(ra, dec)
    inside = ((rr - centre[0]) / sigma[0]) ** 2 + ((dd - centre[1]) / sigma[1]) ** 2 <= n_sigma**2
    values = frac[inside & np.isfinite(frac)]
    return [float(values.min()), float(values.max())] if values.size else None


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--run", default="wu2026_m31pm_scan_n100000")
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--dpi", type=int, default=150)
    args = parser.parse_args()

    pmra, pmdec, merged, box, meta = _scan(args.run)
    ra = np.arange(np.ceil(box[0][0]), box[0][1] + 1e-9, GRID_STEP)
    dec = np.arange(np.ceil(box[1][0]), box[1][1] + 1e-9, GRID_STEP)
    speed, zero_pm = fc.transverse_speed_model()
    wu_centre = fc.MEASUREMENTS["Wu et al. 2025 (2026 forecast input)"]

    scenarios = {
        "target": {"sigma": fc.TARGET_SIGMA, "title": "If new measurements reach 2 \u03bcas/yr"},
        "current": {"sigma": fc.CURRENT_SIGMA, "title": "If precision stays at today's ~8 \u03bcas/yr"},
    }
    report = {
        "scan_run": args.run, "scan_samples": int(merged.size), "scan_box_uas_per_yr": box,
        "reference_run": fc.REFERENCE_RUN, "reference_fraction": fc.REFERENCE_FRACTION,
        "material_change": fc.MATERIAL_CHANGE, "clip_sigma": fc.CLIP_SIGMA,
        "zero_transverse_velocity_pm": zero_pm, "horizon_code_units": 10.0,
        "upstream_commit": meta["upstream_commit"], "scenarios": {},
    }

    validation = fc.predict_point(pmra, pmdec, merged, wu_centre, fc.CURRENT_SIGMA, box)
    combined_se = np.hypot(validation["se"], np.sqrt(fc.REFERENCE_FRACTION * (1 - fc.REFERENCE_FRACTION) / 10000))
    validation.update(reference=fc.REFERENCE_FRACTION,
                      within_two_combined_se=bool(abs(validation["fraction"] - fc.REFERENCE_FRACTION) < 2 * combined_se))
    report["validation_reweighted_fiducial"] = validation

    for key, sc in scenarios.items():
        frac, n_eff = fc.predict_grid(pmra, pmdec, merged, ra, dec, sc["sigma"], box)
        sc["frac"] = frac
        points = {}
        for name, centre in {**fc.MEASUREMENTS, "Zero transverse velocity": zero_pm}.items():
            pt = fc.predict_point(pmra, pmdec, merged, centre, sc["sigma"], box)
            pt["category"] = fc.classify(pt["fraction"]) if np.isfinite(pt["fraction"]) else None
            pt["transverse_speed_km_s"] = float(speed(*centre))
            points[name] = pt
        report["scenarios"][key] = {
            "sigma_uas_per_yr": sc["sigma"], "predictions": points,
            "min_n_eff_in_domain": float(np.nanmin(n_eff)),
            "range_within_current_1sigma": _range_within(frac, ra, dec, wu_centre, fc.CURRENT_SIGMA, 1),
            "range_within_current_2sigma": _range_within(frac, ra, dec, wu_centre, fc.CURRENT_SIGMA, 2),
        }

    (ROOT / "runs" / args.run / "forecast.json").write_text(json.dumps(report, indent=2), encoding="utf-8")

    fig, axes = plt.subplots(1, 2, figsize=(12, 5.6), sharey=True, dpi=args.dpi)
    rr, dd = np.meshgrid(ra, dec)
    vt = speed(rr, dd)
    cmap = ListedColormap(CATEGORY_COLOURS)
    cmap.set_bad("white")
    lo, hi = fc.REFERENCE_FRACTION - fc.MATERIAL_CHANGE, fc.REFERENCE_FRACTION + fc.MATERIAL_CHANGE
    for ax, sc in zip(axes, scenarios.values()):
        frac = sc["frac"]
        cats = np.where(frac <= lo, 0, np.where(frac >= hi, 2, 1)).astype(float)
        cats[~np.isfinite(frac)] = np.nan
        ax.pcolormesh(ra, dec, np.ma.masked_invalid(cats), cmap=cmap, vmin=-0.5, vmax=2.5, shading="nearest")
        if np.nanmin(frac) < fc.COIN_FLIP:
            cs = ax.contour(ra, dec, frac, levels=[fc.COIN_FLIP], colors="#7f1d1d", linewidths=1.6)
            ax.clabel(cs, fmt={fc.COIN_FLIP: "50%"}, fontsize=8)
        vc = ax.contour(ra, dec, vt, levels=[20, 50, 80, 110], colors="#555", linewidths=0.7, linestyles="--")
        ax.clabel(vc, fmt="%d km/s", fontsize=7)
        for n_sig in (1, 2):
            ax.add_patch(Ellipse(wu_centre, 2 * n_sig * fc.CURRENT_SIGMA[0], 2 * n_sig * fc.CURRENT_SIGMA[1],
                                 fill=False, ls=":", lw=1.0, ec="#1e3a8a"))
        ax.plot(*wu_centre, marker="*", ms=13, color="#2563eb", mec="white", ls="none")
        ax.plot(*fc.MEASUREMENTS["Salomon et al. 2021 (2025 forecast input)"], marker="s", ms=8,
                color="#c0392b", mec="white", ls="none")
        ax.plot(*zero_pm, marker="^", ms=9, color="#111", ls="none")
        ax.set_title(sc["title"], fontsize=11, color="#111")
        ax.set_xlabel("Andromeda proper motion east, \u03bc\u03b1* (\u03bcas/yr)")
        ax.set_xlim(*box[0])
        ax.set_ylim(*box[1])
        ax.set_aspect("equal")
    axes[0].set_ylabel("Andromeda proper motion north, \u03bc\u03b4 (\u03bcas/yr)")

    handles = [
        Patch(fc=CATEGORY_COLOURS[2], label=f"Raise: {hi:.1%} or more"),
        Patch(fc=CATEGORY_COLOURS[1], label=f"Unchanged: within 5 points of {fc.REFERENCE_FRACTION:.1%}"),
        Patch(fc=CATEGORY_COLOURS[0], label=f"Lower: {lo:.1%} or less"),
        Line2D([], [], marker="*", ms=12, color="#2563eb", ls="none", label="Wu et al. 2025 (today's input)"),
        Line2D([], [], marker="s", ms=8, color="#c0392b", ls="none", label="Salomon et al. 2021 (2025 input)"),
        Line2D([], [], marker="^", ms=8, color="#111", ls="none", label="Zero sideways velocity"),
        Line2D([], [], ls=":", color="#1e3a8a", label="Today's 1\u03c3 and 2\u03c3 range"),
        Line2D([], [], ls="--", color="#555", label="Sideways speed vs Milky Way"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=4, fontsize=8.5, frameon=False, bbox_to_anchor=(0.5, -0.08))
    fig.text(0.5, -0.12, "Colour = where a new central value would move the 90.1% forecast. "
             f"{merged.size:,} simulated futures, frozen 2026 model; blank areas are outside the scan.",
             ha="center", fontsize=8.5, color="#555")
    fig.tight_layout()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, bbox_inches="tight", facecolor="white", dpi=args.dpi)
    print(json.dumps(report, indent=2))
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
