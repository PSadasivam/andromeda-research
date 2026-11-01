# Andromeda Research

Reproducible forecasts of the Milky Way-Andromeda encounter, built to publish a dated forecast before Gaia Data Release 4 and update it openly afterwards.

**ORCID**: [https://orcid.org/0009-0008-0069-693X](https://orcid.org/0009-0008-0069-693X)

## Overview

Published forecasts disagree: Sawala et al. (2025) found roughly even odds of a merger within 10 Gyr, while Wu et al. (2026), using largely the same model with updated proper motions, found about 90%. This repository reproduces both from the public model and makes every input, seed, and version explicit.

The research plan, story, and milestones live in the portal repository: `deep_space_portal/docs/andromeda-research-plan.md`. This repository holds the numerical work; the portal shows precomputed, versioned results.

## Model

A port of the public Sawala et al. (2025) semi-analytic model ([TillSawala/MW-M31](https://github.com/TillSawala/MW-M31), pinned at commit `e04820e` as a submodule):

- Four galaxies (Milky Way, Andromeda, LMC, M33) as softened point masses with NFW halo properties.
- Chandrasekhar-type dynamical friction, shared between host and satellite in proportion to mass.
- Leapfrog integration, step 0.001 code units, horizon 10 code units; merger flag at 20 kpc separation.
- Monte Carlo sampling of 20 observables from Gaussians truncated at plus or minus 2 sigma.

The port integrates all samples at once (vectorized) and is tested against the upstream notebook's own functions. See [docs/reproduction-log.md](docs/reproduction-log.md) for differences, findings, and validation.

**Pre-DR4 forecast:** [docs/pre-dr4-forecast.md](docs/pre-dr4-forecast.md) states, before Gaia DR4, what the frozen model expects and which new Andromeda proper motions would raise, lower, or leave it unchanged. The frozen code is git tag `pre-dr4-forecast-v1`.

**Time units:** upstream uses G = 43007 (kpc, km/s, 1e10 Msun), so one code time unit is 0.978 Gyr, although it is labelled "Gyr" upstream. Results are reported in both code units (for comparison with the papers) and physical Gyr.

## Repository Structure

```
andromeda-research/
├── README.md
├── LICENSE                     # GPL-3.0 (inherited from the upstream model)
├── pyproject.toml
├── requirements.txt            # Pinned versions used for published runs
├── andromeda/
│   ├── configs.py              # Sawala (2025) and Wu (2026) fiducial priors; M31 proper-motion scan; pinned sky coordinates
│   ├── initial_conditions.py   # Matched-seed Monte Carlo sampling and Galactocentric transform
│   ├── integrator.py           # Vectorized orbit integrator, dynamical friction, merger logic
│   ├── ensemble.py             # Parallel runs, summaries, Wilson intervals, metadata
│   └── forecast.py             # Pre-DR4 threshold, update predictions by reweighting
├── scripts/
│   ├── run_ensemble.py         # Command-line entry point
│   ├── plot_tracks.py          # Separation-track figure
│   └── forecast_map.py         # Pre-DR4 proper-motion map and forecast.json
├── tests/
│   ├── test_upstream_parity.py # Executes upstream notebook cells and compares results
│   └── test_forecast.py        # Scan sampling, reweighting, threshold
├── runs/                       # summary.json and logs per run (arrays are regenerated, not committed)
├── docs/
│   ├── getting-started.md
│   ├── reproduction-log.md
│   └── pre-dr4-forecast.md     # Stage 2: dated forecast, threshold, decision rule
└── upstream/MW-M31/            # Submodule: Sawala et al. public notebook
```

## Sources

- Sawala et al. (2025), *No certainty of a Milky Way-Andromeda collision*, Nature Astronomy 9, 1206. [doi:10.1038/s41550-025-02563-1](https://doi.org/10.1038/s41550-025-02563-1)
- Wu et al. (2026), *The Fate of the Milky Way-Andromeda System: To Merge or Not?*, [arXiv:2603.22863](https://arxiv.org/abs/2603.22863)
- Wu et al. (2025), *Revisiting the proper motions of M31 and M33 using massive supergiant stars with Gaia DR3*, A&A 701, A265. [doi:10.1051/0004-6361/202555477](https://doi.org/10.1051/0004-6361/202555477)

## License

GPL-3.0, because the model is derived from the GPL-3.0 upstream code. See [LICENSE](LICENSE).
