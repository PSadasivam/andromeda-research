"""Tests for the pre-DR4 forecast tools."""

import numpy as np
import pytest

from andromeda import forecast as fc
from andromeda.configs import CONFIGS
from andromeda.initial_conditions import sample_observables

BOX = ((0.0, 80.0), (-50.0, 10.0))


def test_scan_changes_only_m31_proper_motion():
    fiducial = sample_observables(CONFIGS["wu2026"], 200)
    scan = sample_observables(CONFIGS["wu2026_m31pm_scan"], 200)
    for key, values in fiducial.items():
        if key in (("M31", "pmra"), ("M31", "pmdec")):
            mean, sigma = CONFIGS["wu2026"]["M31"][key[1]]
            assert scan[key].min() >= mean - 5 * sigma and scan[key].max() <= mean + 5 * sigma
            assert scan[key].std() > 2 * sigma
        else:
            np.testing.assert_array_equal(scan[key], values)


def test_reweighting_recovers_known_fraction():
    rng = np.random.default_rng(1)
    n = 200_000
    pmra, pmdec = rng.uniform(*BOX[0], n), rng.uniform(*BOX[1], n)
    merged = pmra < 40.0
    centre_on_edge = fc.predict_point(pmra, pmdec, merged, (40.0, -20.0), (8.0, 6.0), BOX)
    assert abs(centre_on_edge["fraction"] - 0.5) < 3 * centre_on_edge["se"]
    far_side = fc.predict_point(pmra, pmdec, merged, (60.0, -20.0), (8.0, 6.0), BOX)
    assert far_side["fraction"] == 0.0


def test_centres_whose_prior_leaves_the_box_are_undefined():
    rng = np.random.default_rng(2)
    pmra, pmdec = rng.uniform(*BOX[0], 1000), rng.uniform(*BOX[1], 1000)
    frac, _ = fc.predict_grid(pmra, pmdec, np.ones(1000, bool), [5.0, 40.0], [-20.0], (8.0, 6.0), BOX)
    assert np.isnan(frac[0, 0]) and frac[0, 1] == 1.0


@pytest.mark.parametrize("p, expected", [(0.951, "raise"), (0.95, "unchanged"), (0.852, "unchanged"),
                                         (0.851, "lower"), (0.30, "lower"), (1.0, "raise")])
def test_classify_uses_five_point_threshold(p, expected):
    assert fc.classify(p) == expected


def test_zero_transverse_velocity_point():
    speed, zero = fc.transverse_speed_model()
    assert speed(*zero) < 1e-6
    assert speed(*fc.MEASUREMENTS["Wu et al. 2025 (2026 forecast input)"]) > speed(*zero)
