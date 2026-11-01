"""Pre-DR4 forecast: the frozen reference, the materiality threshold, and update predictions.

The decision rule itself is stated in docs/pre-dr4-forecast.md. Everything here was fixed
before the proper-motion scan was run.
"""

import numpy as np

from .configs import CONFIGS
from .initial_conditions import to_galactocentric

# Frozen reference: runs/wu2026_fiducial_n10000, 9,010 of 10,000 futures merge within 10 code units.
REFERENCE_RUN = "wu2026_fiducial_n10000"
REFERENCE_FRACTION = 0.901
MATERIAL_CHANGE = 0.05
COIN_FLIP = 0.5
CLIP_SIGMA = 2.0

# Precision scenarios for the pre-stated regions, uas/yr per component.
CURRENT_SIGMA = tuple(CONFIGS["wu2026"]["M31"][k][1] for k in ("pmra", "pmdec"))
TARGET_SIGMA = (2.0, 2.0)  # Wu et al. (2026): the precision needed to settle the question

# Heliocentric (pm_ra_cosdec, pm_dec) in uas/yr, as used by each published forecast.
MEASUREMENTS = {
    "Wu et al. 2025 (2026 forecast input)": (45.9, -20.5),
    "Salomon et al. 2021 (2025 forecast input)": (48.9, -36.9),
}


def classify(p, reference=REFERENCE_FRACTION, threshold=MATERIAL_CHANGE):
    delta = round(p - reference, 9)  # so a boundary value like 0.951 is not lost to float error
    if delta >= threshold:
        return "raise"
    if delta <= -threshold:
        return "lower"
    return "unchanged"


def _gaussian_weights(samples, centres, sigma, clip):
    z = (np.asarray(centres)[:, None] - samples[None, :]) / sigma
    return np.where(np.abs(z) <= clip, np.exp(-0.5 * z * z), 0.0)


def predict_grid(pmra, pmdec, merged, centres_ra, centres_dec, sigma, box, clip=CLIP_SIGMA):
    """Merger fraction expected if M31's proper-motion prior became N(centre, sigma) truncated at clip.

    Importance-reweights a scan whose M31 proper motions are uniform over `box`
    ((ra_lo, ra_hi), (dec_lo, dec_hi)); every other input keeps its fiducial prior, so the
    uniform density cancels. Returns fraction and effective sample size, shaped (dec, ra);
    centres whose truncated prior leaves the box are NaN.
    """
    w_ra = _gaussian_weights(pmra, centres_ra, sigma[0], clip)
    w_dec = _gaussian_weights(pmdec, centres_dec, sigma[1], clip)
    m = np.asarray(merged, dtype=float)
    den = w_dec @ w_ra.T
    with np.errstate(invalid="ignore", divide="ignore"):
        frac = ((w_dec * m) @ w_ra.T) / den
        n_eff = den**2 / ((w_dec**2) @ (w_ra**2).T)
    (ra_lo, ra_hi), (dec_lo, dec_hi) = box
    ok_ra = (np.asarray(centres_ra) - clip * sigma[0] >= ra_lo) & (np.asarray(centres_ra) + clip * sigma[0] <= ra_hi)
    ok_dec = (np.asarray(centres_dec) - clip * sigma[1] >= dec_lo) & (np.asarray(centres_dec) + clip * sigma[1] <= dec_hi)
    inside = ok_dec[:, None] & ok_ra[None, :]
    return np.where(inside, frac, np.nan), np.where(inside, n_eff, np.nan)


def predict_point(pmra, pmdec, merged, centre, sigma, box, clip=CLIP_SIGMA):
    frac, n_eff = predict_grid(pmra, pmdec, merged, [centre[0]], [centre[1]], sigma, box, clip)
    p, n = float(frac[0, 0]), float(n_eff[0, 0])
    return {"fraction": p, "n_eff": n, "se": float(np.sqrt(p * (1 - p) / n)) if np.isfinite(p) else np.nan}


def transverse_speed_model(galaxy="M31", config="wu2026"):
    """M31's Galactocentric velocity is affine in its proper motion at fixed distance and vr.

    Returns a function giving the speed perpendicular to the MW-M31 line (km/s) for arrays of
    proper motions, and the proper motion at which that speed vanishes.
    """
    p = CONFIGS[config][galaxy]
    dmod, vr = np.full(3, p["dmod"][0]), np.full(3, p["vr"][0])
    pos, vel = to_galactocentric(galaxy, dmod, np.array([0.0, 1.0, 0.0]), np.array([0.0, 0.0, 1.0]), vr)
    r_hat = pos[0] / np.linalg.norm(pos[0])
    b = vel[0]
    a = np.stack([vel[1] - b, vel[2] - b], axis=1)
    perp = np.eye(3) - np.outer(r_hat, r_hat)
    zero = np.linalg.lstsq(perp @ a, -perp @ b, rcond=None)[0]

    def speed(pmra, pmdec):
        v = np.stack([np.asarray(pmra), np.asarray(pmdec)], axis=-1) @ a.T + b
        return np.linalg.norm(v @ perp.T, axis=-1)

    return speed, tuple(float(x) for x in zero)
