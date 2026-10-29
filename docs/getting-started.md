# Getting Started - Andromeda Research

## Prerequisites

- Python 3.11 or newer (published runs used 3.13).
- Git, to fetch the upstream model submodule.
- Multi-core CPU recommended: 10,000 samples take about 12 minutes on 14 workers.

## 1. Environment Setup

```powershell
git clone --recurse-submodules <repository-url> andromeda-research
cd andromeda-research
python -m venv .venv
.venv\Scripts\activate
pip install -r requirements.txt
pip install -e .
```

If you cloned without `--recurse-submodules`, run `git submodule update --init`.

## 2. Run the Tests

```powershell
python -m pytest
```

The tests execute the upstream notebook's own functions and require our port to match them exactly. Sky coordinates are pinned in `andromeda/configs.py`, so no network access is needed.

## 3. Run an Ensemble

```powershell
python scripts/run_ensemble.py --config wu2026 --n 10000 --workers 14
python scripts/run_ensemble.py --config sawala2025 --n 10000 --workers 14
```

Options:

| Option | Default | Meaning |
|--------|---------|---------|
| `--config` | required | `sawala2025` or `wu2026` fiducial priors |
| `--galaxies` | `MW,M31,LMC,M33` | Galaxies to include; MW and M31 are required |
| `--n` | 10000 | Monte Carlo samples |
| `--workers` | 1 | Parallel processes |
| `--friction` | `proportional` | `proportional`, `hierarchical`, or `none` |
| `--softening` | 20 | Gravitational softening (kpc) |
| `--merger-threshold` | 20 | Merger separation (kpc) |

Results go to `runs/<name>/`: `summary.json` (statistics and full metadata) and `ensemble.npz` (per-sample arrays, not committed).

## 4. Reading the Summary

- `fraction` and `fraction_95ci`: merger fraction within the horizon, with a Wilson 95% interval for Monte Carlo sampling error only.
- `merge_time_*`: median and 16th-84th percentile offsets, among merging samples only, in code units and physical Gyr.
- `m33_absorbed_before_mw_m31` / `_after_`: whether M33 is absorbed before or after the Milky Way-Andromeda merger.

Sampling intervals do not include observational, systematic, or model uncertainty.
