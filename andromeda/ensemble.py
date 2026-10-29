"""Ensemble execution and summary statistics."""

import json
import platform
import subprocess
import time
from concurrent.futures import ProcessPoolExecutor
from importlib.metadata import version
from pathlib import Path

import numpy as np

from .configs import CONFIGS, SOURCES
from .initial_conditions import build_initial_conditions
from .integrator import TIME_UNIT_GYR, integrate

UPSTREAM_DIR = Path(__file__).resolve().parents[1] / "upstream" / "MW-M31"


def _integrate_chunk(args):
    ics, kwargs = args
    return integrate(ics, **kwargs)


def _slice_ics(ics, sl):
    return {k: (v[sl] if isinstance(v, np.ndarray) else v) for k, v in ics.items() if k != "observables"}


def run_ensemble(config, galaxies, n, workers=1, clip_sigma=2.0, **kwargs):
    started = time.time()
    ics = build_initial_conditions(CONFIGS[config], galaxies, n, clip_sigma)
    bounds = np.linspace(0, n, max(1, workers) + 1, dtype=int)
    chunks = [(_slice_ics(ics, slice(a, b)), dict(kwargs, track_samples=max(0, min(b, 200) - a)))
              for a, b in zip(bounds[:-1], bounds[1:]) if b > a]
    if workers > 1:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            parts = list(pool.map(_integrate_chunk, chunks))
    else:
        parts = [_integrate_chunk(c) for c in chunks]

    result = {k: v for k, v in parts[0].items() if not isinstance(v, np.ndarray)}
    for key, value in parts[0].items():
        if isinstance(value, np.ndarray):
            result[key] = np.concatenate([p[key] for p in parts])
    result["observables"] = {f"{g}_{k}": v for (g, k), v in ics["observables"].items()}
    result["metadata"] = _metadata(config, galaxies, n, clip_sigma, workers, time.time() - started)
    return result


def _metadata(config, galaxies, n, clip_sigma, workers, seconds):
    try:
        commit = subprocess.run(["git", "-C", str(UPSTREAM_DIR), "rev-parse", "HEAD"],
                                capture_output=True, text=True, check=True).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        commit = "unknown"
    return {
        "config": config, "config_source": SOURCES[config], "parameters": CONFIGS[config],
        "galaxies": list(galaxies), "n_samples": n, "clip_sigma": clip_sigma,
        "upstream_commit": commit, "workers": workers, "runtime_seconds": round(seconds, 1),
        "python": platform.python_version(),
        "packages": {p: version(p) for p in ("numpy", "scipy", "astropy", "colossus")},
        "time_unit_gyr": TIME_UNIT_GYR,
    }


def wilson_interval(k, n, z=1.959964):
    p = k / n
    centre = (p + z**2 / (2 * n)) / (1 + z**2 / n)
    half = z * np.sqrt(p * (1 - p) / n + z**2 / (4 * n**2)) / (1 + z**2 / n)
    return centre - half, centre + half


def _quantiles(values):
    values = values[np.isfinite(values)]
    if values.size == 0:
        return None
    q16, q50, q84 = np.percentile(values, [16, 50, 84])
    return {"median": q50, "minus": q50 - q16, "plus": q84 - q50, "n": int(values.size)}


def summarize(result):
    n = result["metadata"]["n_samples"]
    summary = {"n_samples": n, "horizon_code_units": result["t_end"],
               "horizon_gyr": result["t_end"] * TIME_UNIT_GYR, "pairs": {}}
    for key in [k[len("merge_time_"):] for k in result if k.startswith("merge_time_")]:
        times = result[f"merge_time_{key}"]
        k = int(np.isfinite(times).sum())
        lo, hi = wilson_interval(k, n)
        summary["pairs"][key] = {
            "merged": k, "fraction": k / n, "fraction_95ci": [lo, hi],
            "merge_time_code_units": _quantiles(times),
            "merge_time_gyr": _quantiles(times * TIME_UNIT_GYR),
        }
    summary["first_pericentre_time_code_units"] = _quantiles(result["first_pericentre_time"])
    summary["first_pericentre_distance_kpc"] = _quantiles(result["first_pericentre_distance"])
    if "absorbed_time_M33" in result:
        m33, mw = result["absorbed_time_M33"], result["merge_time_MW-M31"]
        before = np.isfinite(m33) & (~np.isfinite(mw) | (m33 < mw))
        summary["m33_absorbed_before_mw_m31"] = float(before.mean())
        summary["m33_absorbed_after_mw_m31"] = float((np.isfinite(m33) & ~before).mean())
    return summary


def save(result, summary, out_dir):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    arrays = {k: v for k, v in result.items() if isinstance(v, np.ndarray)}
    arrays.update({f"obs_{k}": v for k, v in result["observables"].items()})
    np.savez_compressed(out_dir / "ensemble.npz", **arrays)
    (out_dir / "summary.json").write_text(
        json.dumps({"metadata": result["metadata"], "summary": summary}, indent=2, default=float),
        encoding="utf-8")
