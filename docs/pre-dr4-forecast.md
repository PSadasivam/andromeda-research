# Pre-DR4 Forecast: Milky Way-Andromeda Merger

Stage 2 of the research plan. This note commits, before Gaia Data Release 4 (ESA date: 2 December 2026), to what the frozen model expects and what new measurements of Andromeda's proper motion would change it.

- **Author:** Prabhu Sadasivam ([ORCID 0009-0008-0069-693X](https://orcid.org/0009-0008-0069-693X)).
- **Threshold and decision rule fixed:** 2026-11-01, committed before the proper-motion scan was run.
- **Frozen code:** the git tag on the commit that contains this note's results.
- **Our contribution:** the advance commitment, open reproducibility, and the public update. The proper-motion probability map itself was first published by Wu et al. (2026, Fig. 5); this is a reproduction of their approach, not a new model.

## 1. The Reference Forecast

The frozen configuration is the Wu et al. (2026) fiducial four-body model (MW, M31, LMC, M33), reproduced in [reproduction-log.md](reproduction-log.md):

- **Merger within the horizon:** 90.1% (9,010 of 10,000 futures; Wilson 95% sampling interval 89.5-90.7%).
- **Horizon:** 10 code units = 9.78 Gyr. Merger flag at 20 kpc, softening 20 kpc, proportional dynamical friction, step 0.001 code units.
- **Inputs:** `andromeda/configs.py`, `WU2026`. M31 proper motion (45.9 +/- 8.1, -20.5 +/- 6.6) uas/yr, heliocentric; priors truncated at +/- 2 sigma; upstream seed order.
- **Run:** `runs/wu2026_fiducial_n10000/summary.json`.

The 90.1% is conditional on this model and these inputs. It is not a model-independent probability, and it must not be averaged with the 52.3% Sawala et al. (2025) result.

## 2. What Counts as a Material Change

Fixed before any scan result was seen:

| Category | Updated merger fraction |
|----------|-------------------------|
| **Raise** | 95.1% or higher (at least 5 percentage points above 90.1%) |
| **Unchanged** | between 85.1% and 95.1% |
| **Lower** | 85.1% or lower (at least 5 percentage points below 90.1%) |

A result below 50% is additionally flagged as **"coin flip or worse"**, the 2025 headline.

Why 5 percentage points: it is about 17 times the Monte Carlo standard error of a 10,000-sample run near 90% (0.3 points), larger than the headline MW-M31 residuals in the Stage 1 reproduction (0.1 and 0.5 points), and small enough that the categories still discriminate, since 90.1% leaves only 9.9 points of headroom. If an updated value lies within two combined sampling standard errors of a boundary, it is reported as borderline as well as categorized.

## 3. The Decision Rule

1. **Trigger.** The first M31 proper motion derived from Gaia DR4 data, posted on arXiv or published in a refereed journal, with stated uncertainties. Later DR4-based values are reported as separate, labelled updates. No update is made from DR4 release-day data alone.
2. **Substitution.** Replace only the M31 proper-motion mean and sigma in `WU2026`. If statistical and systematic errors are given separately, combine them in quadrature per component. Convert to heliocentric (pm_ra_cosdec, pm_dec) in uas/yr if needed, documenting the conversion.
3. **Nothing else changes.** Same code at the tag, same priors for every other input (including M33 and the M31 concentration), same seeds, 10,000 samples, +/- 2 sigma truncation, horizon, merger flag, softening, and friction. The frozen pipeline samples the two components independently; a reported correlation is ignored in the first update and used only in a separate, labelled update.
4. **Report.** Publish the updated fraction beside 90.1%, its category from Section 2, and the prediction this note made for that value (Section 4), with any disagreement between prediction and re-run explained.
5. **M33, separately.** If a DR4-based M33 proper motion appears, run M33 alone (M31 unchanged) and then both together, each labelled.
6. **Model changes are separate.** Topic 9 fixes, mass priors, halo profiles, or any other change are later, labelled updates, never folded into the first.

## 4. Pre-Stated Regions

Computed by `scripts/forecast_map.py` from the 100,000-sample scan (`runs/wu2026_m31pm_scan_n100000/`, 591 s on 14 workers); full numbers in its `forecast.json`. Figure: `andromeda_forecast_map.png` in the portal repository.

**Check passed.** Reweighting the scan to today's measurement gives 90.0% +/- 0.3% (effective sample size 10,410), against the reference 90.1%.

**Predicted updated fraction** for three possible new central values (sampling standard error in brackets):

| New central value (heliocentric, uas/yr) | Sideways speed vs MW | New precision 2 uas/yr | Today's precision (8.1, 6.6) |
|---|---|---|---|
| Today's value, Wu et al. 2025 (45.9, -20.5) | 26 km/s | 95.4% (0.8): **raise**, borderline | 90.0% (0.3): unchanged |
| The 2025 input, Salomon et al. 2021 (48.9, -36.9) | 67 km/s | 69.3% (1.7): **lower** | 55.9% (0.5): **lower** |
| Zero sideways velocity (38.8, -21.4) | 0 km/s | 91.0% (1.0): unchanged | 88.7% (0.3): unchanged |

Sideways speeds are in the model's Galactocentric frame at central distance and radial velocity.

**Stated in advance:**

1. **Precision alone moves the number.** If DR4-based measurements confirm today's central value at about 2 uas/yr, the forecast rises to about 95%, at the raise boundary. Today's 90.1% averages over a wide uncertainty range that includes low-probability motions; a sharper measurement at the same centre removes them.
2. **Most of today's 1 sigma range leaves the forecast unchanged** at 2 uas/yr precision. The raise region is only a small patch around and just below today's central value. Across today's 1 sigma and 2 sigma ranges, the 2 uas/yr forecast spans 86.9-95.8% and 69.9-95.8%.
3. **A return to "coin flip or worse" needs a large shift.** At 2 uas/yr the 50% contour lies entirely outside today's 2 sigma range, at sideways speeds of roughly 80-110 km/s. A DR4 value back at the 2025 input lowers the forecast to about 69% at 2 uas/yr, or 56% at today's precision, not to the 2025 result of 52%: the 2026 configuration also changed M31's concentration and M33's motion, and those changes stay frozen.
4. **Head-on is not the maximum.** Zero sideways velocity gives about 91%, below the 95% at today's central value: the LMC's pull on the Milky Way makes a slightly oblique M31 orbit more likely to merge than a head-on one (Wu et al. 2026, Sec. III.2).
5. **If precision does not improve,** no central value in the scanned range raises the forecast, and only values within about 1 sigma of today's leave it unchanged.

**Comparison with Wu et al. (2026).** They report 73.3-100% and 64.7-100% across the 1 sigma and 2 sigma regions, from fixed-proper-motion grid cells. Our ranges average over a 2 uas/yr prior instead, which smooths the extremes; the ranges are related but not the same estimator. Near category boundaries the map's sampling error (0.7-1.7 points at 2 uas/yr) matters, which is why borderline results are reported as such.

## 5. Method for the Map

- **Scan.** Config `wu2026_m31pm_scan`: the frozen Wu (2026) configuration with the M31 proper-motion prior replaced by a uniform box of +/- 5 sigma per component, as in Wu et al. (2026, Sec. III.2). Every other input keeps its fiducial prior and seed, so the scan differs from the reference only in M31's proper motion. 100,000 samples.
- **Prediction for any new measurement.** For a new central value and uncertainty, the expected updated fraction is the scan's merger fraction importance-weighted by the new truncated Gaussian prior (`andromeda/forecast.py`, `predict_point`). Because the scan's prior is uniform, the weights are simply the new prior's density. This predicts what the decision-rule re-run will give, up to sampling error, provided the new prior lies inside the scanned box.
- **Two precision scenarios.** Regions are drawn for 2 uas/yr per component, the precision Wu et al. (2026) estimate is needed to settle the question, and for today's precision (8.1, 6.6 uas/yr). The category depends on precision as well as on the central value, because the forecast averages over the whole uncertainty range.
- **Check.** Reweighting the scan to today's measurement must recover the reference 90.1% within two combined sampling standard errors.

## 6. Limitations

- **One model.** Rigid NFW halos, point-mass orbits, a 20 kpc merger flag, and Chandrasekhar friction applied to near-equal galaxies. The friction term reverses sign between 20 and about 32 kpc ([reproduction-log.md](reproduction-log.md), Model Findings 2), near typical first pericentres. These are Stage 3 stress tests and do not enter the first update.
- **Conditional, not absolute.** All fractions are conditional on the model, priors, and horizon. "No merger within 9.78 Gyr" is not "no merger ever".
- **Frame.** The Sun's Galactocentric position and velocity are fixed (Astropy v4.0 frame). Sideways speeds quoted here use that frame at central distance and radial velocity.
- **Gaia DR4 does not publish an M31 proper motion.** Deriving one requires careful source selection and systematics control; the update waits for it.
