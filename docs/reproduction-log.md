# Reproduction Log

Stage 1 of the research plan: reproduce the published fiducial forecasts before building anything new.

- **Date:** 2026-10-28.
- **Upstream:** [TillSawala/MW-M31](https://github.com/TillSawala/MW-M31) at commit `e04820ec366c57a9f03ac3af14bcf62b5467eb2e`.
- **Environment:** Python 3.13.14; NumPy 2.5.4, SciPy 1.18.1, Astropy 8.0.2, Colossus 1.4.0 (`requirements.txt`). Upstream used Python 3.10, NumPy 1.26.4, SciPy 1.13.0, Astropy 6.0.1, Colossus 1.3.5.
- **Settings:** four galaxies (MW, M31, LMC, M33); proportional dynamical friction; softening 20 kpc; merger threshold 20 kpc; step 0.001; horizon 10 code units; priors truncated at plus or minus 2 sigma; 10,000 samples with upstream seed order.

## Results

Times are in code units (1 unit = 0.978 Gyr), matching how both papers report them. Intervals on fractions are Wilson 95% intervals for Monte Carlo sampling only; times are medians with 16th-84th percentile offsets among merging samples.

### Wu et al. (2026) fiducial

| Quantity | This reproduction | Wu et al. (2026) |
|----------|-------------------|------------------|
| MW-M31 merger fraction within horizon | 90.1% [89.5, 90.7] | 90% |
| MW-M31 merger time | 6.53 -1.53/+1.27 | 6.5 -1.5/+1.3 |
| LMC-MW merger time | 1.27 -0.37/+0.61 | 1.3 -0.4/+0.6 |
| M33 absorbed before / after MW-M31 merger | 20.7% / 5.6% | 21% / 5% |
| First pericentre time | 5.34 -0.60/+1.00 | 5.3 -0.6/+0.9 |
| First pericentre distance | 26.2 -6.2/+26.8 kpc | 28.9 -8.9/+33.8 kpc |

Runtime: 82 s on 14 workers. Run: `runs/wu2026_fiducial_n10000/`.

### Sawala et al. (2025) fiducial

| Quantity | This reproduction | Sawala et al. (2024 preprint) |
|----------|-------------------|-------------------------------|
| MW-M31 merger fraction, 10,000 samples | 52.3% [51.4, 53.3] | "just over 50%" |
| MW-M31 merger fraction, first 2,500 samples | 53.5% | 54% (Fig. 2, 2,500 samples) |
| MW-M31 merger time | 7.52 -1.06/+1.35 | 7.6 (median) |
| LMC-MW merger time | 1.27 -0.37/+0.61 | 1.3 (median) |
| M33-M31 merger fraction / median time | 89.5% / 3.32 | about 86% / 3.3 |
| M33 absorbed before MW-M31 merger | 86.0% | about 86% |

Runtime: 101 s on 14 workers. Run: `runs/sawala2025_fiducial_n10000/`. The preprint values are used until the final Nature Astronomy text has been checked against them (Stage 3).

## Exit-Gate Assessment

**Passed.** All headline quantities agree within two combined sampling standard errors; for a 2,500-sample published fraction near 50%, that tolerance is about 2.8 percentage points.

- The 0.5-point residual on the matched 2,500-sample Sawala comparison is plausibly a library-version effect: the seed order is identical, but random draws for the same seed can differ across SciPy and NumPy versions. Confirming this requires a run in the upstream's original environment (Python 3.10, NumPy 1.26, SciPy 1.13); it is a hypothesis, not yet tested.
- The first-pericentre distance differs by 2.7 kpc at the median. Wu et al. do not define their pericentre detector, and ours (below) is designed to ignore jumps caused by satellite mergers. Treat pericentre distances as definition-dependent until Stage 3.

## Port Differences From Upstream

None of these change results relative to the upstream functions, as verified by `tests/test_upstream_parity.py`, except where noted.

1. **Pinned sky coordinates.** Upstream resolves galaxy positions online by name (`SkyCoord.from_name`). Coordinates are pinned in `andromeda/configs.py` from the same resolver on 2026-10-28. During setup the resolver timed out for the SMC, which shows why runs should not depend on it.
2. **Vectorized integration.** All samples are integrated simultaneously instead of one at a time, giving the same arithmetic per sample.
3. **Analytic NFW properties.** Scale radius, r200, and Vmax are computed analytically and tested against Colossus to 1e-6.
4. **Absorbed galaxies.** A zero-mass guard in the friction share prevents NaN values spreading through the vectorized gravity sum. The guard was needed only because of vectorization; upstream skips absorbed galaxies before reaching this code.
5. **Satellite fates.** Upstream marks a pair "merged" when either galaxy disappears, so an M33-M31 "merger" is also recorded when M31 merges into the Milky Way. We keep that statistic for comparison and add per-galaxy absorption times, which match the M33 before/after percentages in Wu et al. (2026).
6. **First pericentre.** The first MW-M31 pericentre is the running minimum separation, confirmed once the separation rises 20 kpc above it. A naive "first increase" detector misfires when the LMC merges into the Milky Way and shifts its centre by a few kpc. Samples that merge before rebounding are assigned a pericentre equal to the merger threshold, an upper bound.
7. **SMC omitted.** Neither fiducial configuration includes the SMC.

## Model Findings for Act III

These are properties of the upstream model, recorded for the Topic 9 stress tests. They are not errors in the published papers' reproduction; both papers inherit them.

1. **Time units.** With G = 43007 in kpc, km/s, and 1e10 Msun, one code time unit is 0.978 Gyr, but upstream labels it "Gyr". The "10 Gyr" horizon is 9.78 Gyr, and quoted merger times are about 2.3% high in physical units: the Wu median of 6.5 is 6.4 Gyr. Small, but every public number must say which unit it uses.
2. **Friction reverses sign at close range.** The Coulomb logarithm is ln(r/epsilon), with epsilon = 2.2 r_s - 14 kpc for the satellite. For a Milky Way halo of 1e12 Msun at concentration 10, r_s = 21.1 kpc, so epsilon = 32.4 kpc. Between the 20 kpc merger flag and about 32 kpc the logarithm is negative and "friction" accelerates instead of decelerating, exactly where typical first pericentres fall (about 26-29 kpc). The size of the effect on merger fractions is untested; the 20 versus 30 kpc threshold comparisons in Sawala et al. bound it only partly. Testing it is a priority Topic 9 check.
3. **Halo properties frozen at merger.** After a merger, the remnant keeps the absorbing halo's scale radius and concentration while its mass grows, which affects later friction on M33.
4. **Merger-loop aliasing.** Upstream iterates over a list while removing from it (`exists_const = exists` is an alias, not a copy). This can skip pair checks in the rare step where two mergers happen together. Our port checks all pairs, matching upstream on all tested samples.
5. **Fixed solar parameters.** The Sun's Galactocentric distance and velocity are fixed (Astropy v4.0 frame: R0 = 8.122 kpc, v_sun = (12.9, 245.6, 7.78) km/s), not sampled. Topic 9 will test whether sampling them matters.

## Validation

- `test_sawala_config_matches_upstream_fiducial`: our Sawala priors equal upstream `DefineParameters('fiducial')`.
- `test_initial_conditions_match_upstream`: positions, velocities, masses, and concentrations match upstream `CreateICs`.
- `test_nfw_properties_match_colossus`: analytic halo properties match Colossus.
- `test_orbits_and_mergers_match_upstream`: four full integrations, covering MW-M31, LMC, and M33 mergers, match upstream `IntegrateOrbits` and `calculate_mergers` merger times exactly and minimum separations to 1e-6.
- Full ensembles ran with NumPy runtime warnings promoted to errors.

## Next

Stage 2 is in [pre-dr4-forecast.md](pre-dr4-forecast.md): the materiality threshold, the decision rule, and the M31 proper-motion map in the frozen Wu (2026) configuration.
