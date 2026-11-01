"""Monte Carlo initial conditions in the Galactocentric frame (kpc, km/s)."""

import astropy.units as u
import numpy as np
from astropy.coordinates import Distance, Galactocentric, SkyCoord, galactocentric_frame_defaults
from scipy.stats import truncnorm

from .configs import GALAXY_ORDER, SKY_COORDS_DEG

# Upstream draws parameter k after np.random.seed(k); keeping this order keeps draws matched across configurations.
SEED_ORDER = (
    ("MW", "m"), ("MW", "c"),
    ("M31", "dmod"), ("M31", "pmra"), ("M31", "pmdec"), ("M31", "vr"), ("M31", "m"), ("M31", "c"),
    ("LMC", "dmod"), ("LMC", "pmra"), ("LMC", "pmdec"), ("LMC", "vr"), ("LMC", "m"), ("LMC", "c"),
    ("M33", "dmod"), ("M33", "pmra"), ("M33", "pmdec"), ("M33", "vr"), ("M33", "m"), ("M33", "c"),
)


def sample_observables(params, n, clip_sigma=2.0):
    draws = {}
    for seed, (galaxy, key) in enumerate(SEED_ORDER):
        spec = params[galaxy][key]
        np.random.seed(seed)
        if spec[0] == "uniform":
            draws[galaxy, key] = np.random.uniform(spec[1], spec[2], n)
        else:
            mean, sigma = spec
            draws[galaxy, key] = truncnorm(-clip_sigma, clip_sigma, loc=mean, scale=sigma).rvs(n)
    return draws


def to_galactocentric(galaxy, dmod, pmra, pmdec, vr):
    ra, dec = SKY_COORDS_DEG[galaxy]
    n = len(dmod)
    sky = SkyCoord(
        ra=np.full(n, ra) * u.deg, dec=np.full(n, dec) * u.deg,
        distance=Distance(distmod=np.asarray(dmod) * u.mag),
        pm_ra_cosdec=np.asarray(pmra) * u.uas / u.yr, pm_dec=np.asarray(pmdec) * u.uas / u.yr,
        radial_velocity=np.asarray(vr) * u.km / u.s, frame="icrs",
    )
    with galactocentric_frame_defaults.set("v4.0"):
        gal = sky.transform_to(Galactocentric())
    pos = np.stack([gal.x.to_value(u.kpc), gal.y.to_value(u.kpc), gal.z.to_value(u.kpc)], axis=-1)
    vel = np.stack([gal.v_x.to_value(u.km / u.s), gal.v_y.to_value(u.km / u.s),
                    gal.v_z.to_value(u.km / u.s)], axis=-1)
    return pos, vel


def build_initial_conditions(params, galaxies, n, clip_sigma=2.0):
    names = [g for g in GALAXY_ORDER if g in galaxies]
    if names[:2] != ["MW", "M31"] or len(names) != len(set(galaxies)):
        raise ValueError(f"galaxies must include MW and M31 and come from {GALAXY_ORDER}")

    draws = sample_observables(params, n, clip_sigma)
    positions = np.zeros((n, len(names), 3))
    velocities = np.zeros_like(positions)
    masses = np.empty((n, len(names)))
    concentrations = np.empty_like(masses)

    for k, name in enumerate(names):
        masses[:, k] = draws[name, "m"]
        concentrations[:, k] = draws[name, "c"]
        if name != "MW":
            positions[:, k], velocities[:, k] = to_galactocentric(
                name, draws[name, "dmod"], draws[name, "pmra"], draws[name, "pmdec"], draws[name, "vr"])

    return {"names": names, "positions": positions, "velocities": velocities,
            "masses": masses, "concentrations": concentrations, "observables": draws}
