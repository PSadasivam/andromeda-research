"""Parity tests: our vectorized port against the upstream notebook's own functions."""

import json
import math
from pathlib import Path

import astropy.units as u
import numpy as np
import pytest
from astropy.coordinates import Distance, Galactocentric, SkyCoord, galactocentric_frame_defaults
from colossus.cosmology import cosmology
from colossus.halo import profile_nfw
from scipy.stats import truncnorm

from andromeda.configs import CONFIGS, SKY_COORDS_DEG
from andromeda.initial_conditions import build_initial_conditions
from andromeda.integrator import integrate, nfw_halo_properties

NOTEBOOK = Path(__file__).resolve().parents[1] / "upstream" / "MW-M31" / "galaxies.ipynb"
UPSTREAM_CELLS = (10, 12, 13, 15, 16, 18, 19, 20, 22, 23)
GALAXIES = ["MW", "M31", "LMC", "M33"]


class PinnedSkyCoord(SkyCoord):
    @classmethod
    def from_name(cls, name, *args, **kwargs):
        ra, dec = SKY_COORDS_DEG[name]
        return SkyCoord(ra=ra * u.deg, dec=dec * u.deg)


@pytest.fixture(scope="module")
def upstream():
    if not NOTEBOOK.exists():
        pytest.skip("upstream notebook not cloned; see docs/getting-started.md")
    galactocentric_frame_defaults.set("v4.0")
    cells = json.loads(NOTEBOOK.read_text(encoding="utf-8"))["cells"]
    t_end, dt = 10, 0.001
    n_t = int(np.ceil(t_end / dt))
    ns = {
        "np": np, "math": math, "truncnorm": truncnorm, "u": u, "Galactocentric": Galactocentric,
        "SkyCoord": PinnedSkyCoord, "Distance": Distance, "profile_nfw": profile_nfw,
        "cosmo": cosmology.setCosmology("planck18"), "tqdm": lambda x: x, "G": 43007,
        "dt": dt, "tEnd": t_end, "Nt": n_t, "t_all": np.arange(n_t + 1) * dt,
        "softening": 20, "merger_threshold": 20, "save_observables": False,
    }
    for idx in UPSTREAM_CELLS:
        exec("".join(cells[idx]["source"]), ns)
    return ns


def test_sawala_config_matches_upstream_fiducial(upstream):
    values = upstream["DefineParameters"]("fiducial", GALAXIES)
    cfg = CONFIGS["sawala2025"]
    # Upstream returns MW (m, c), then m, dmod, pmra, pmdec, vr, c for M31, M33, LMC, each as (mean, sigma).
    expected = [*cfg["MW"]["m"], *cfg["MW"]["c"]]
    for galaxy in ("M31", "M33", "LMC"):
        for key in ("m", "dmod", "pmra", "pmdec", "vr", "c"):
            expected.extend(cfg[galaxy][key])
    assert values[:40] == expected


def test_initial_conditions_match_upstream(upstream):
    n = 6
    params = upstream["DefineParameters"]("fiducial", GALAXIES)
    ref = upstream["CreateICs"](GALAXIES, params, 2.0, n, "")
    ours = build_initial_conditions(CONFIGS["sawala2025"], GALAXIES, n)
    np.testing.assert_allclose(ours["masses"], ref["masses"], rtol=1e-12)
    np.testing.assert_allclose(ours["concentrations"], ref["concentrations"], rtol=1e-12)
    np.testing.assert_allclose(ours["positions"], ref["positions"], rtol=1e-9, atol=1e-9)
    np.testing.assert_allclose(ours["velocities"], ref["velocities"], rtol=1e-9, atol=1e-9)


def test_nfw_properties_match_colossus():
    cosmo = cosmology.setCosmology("planck18")
    masses = np.array([[15.0, 100.0, 130.0, 300.0]])
    concs = np.array([[7.0, 10.0, 15.0, 12.0]])
    r_s, r_200, vmax = nfw_halo_properties(masses, concs)
    for k in range(masses.shape[1]):
        prof = profile_nfw.NFWProfile(M=masses[0, k] * 1e10 * cosmo.h, c=concs[0, k], z=0.0, mdef="200c")
        assert r_s[0, k] == pytest.approx(prof.par["rs"] / cosmo.h, rel=1e-9)
        assert r_200[0, k] == pytest.approx(prof.par["rs"] / cosmo.h * concs[0, k], rel=1e-9)
        assert vmax[0, k] == pytest.approx(prof.Vmax()[0], rel=1e-6)


@pytest.mark.slow
def test_orbits_and_mergers_match_upstream(upstream):
    n = 4
    params = upstream["DefineParameters"]("fiducial", GALAXIES)
    ref_ics = upstream["CreateICs"](GALAXIES, params, 2.0, n, "")
    orbits, _, existence = upstream["IntegrateOrbits"](ref_ics, "proportional")
    dmin, t_m, _, dmin_lmc, t_m_lmc, _, dmin_m33, t_m_m33, _ = upstream["calculate_mergers"](
        GALAXIES, orbits, existence, n)

    ours = integrate(build_initial_conditions(CONFIGS["sawala2025"], GALAXIES, n))
    for ref_t, key in ((t_m, "MW-M31"), (t_m_lmc, "LMC-MW"), (t_m_m33, "M33-M31")):
        expected = np.where(ref_t > 10.5, np.nan, ref_t)
        np.testing.assert_allclose(ours[f"merge_time_{key}"], expected, atol=1e-9, equal_nan=True)
    live = np.isnan(ours["merge_time_MW-M31"])
    np.testing.assert_allclose(ours["min_separation_MW-M31"][live], dmin[live], rtol=1e-6)
