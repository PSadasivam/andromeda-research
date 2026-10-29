"""Observational priors as (mean, sigma) pairs.

Units: M200 in 1e10 Msun, distance modulus in mag, proper motions in uas/yr,
heliocentric line-of-sight velocity in km/s.
"""

from copy import deepcopy

# ICRS centres from the CDS Sesame resolver (as queried upstream), resolved 2026-10-28 and pinned.
SKY_COORDS_DEG = {
    "M31": (10.68470833, 41.26875),
    "M33": (23.46206906, 30.66017511),
    "LMC": (80.89416667, -69.75611111),
}

# Upstream array order; the forecast always needs MW and M31 at indices 0 and 1.
GALAXY_ORDER = ("MW", "M31", "LMC", "M33")

SAWALA2025 = {
    "MW": {"m": (100.0, 20.0), "c": (10.0, 2.0)},
    "M31": {
        "m": (130.0, 40.0), "c": (10.0, 2.0), "dmod": (24.407, 0.032),
        "pmra": (48.9, 10.5), "pmdec": (-36.9, 8.1), "vr": (-301.0, 1.0),
    },
    "M33": {
        "m": (30.0, 10.0), "c": (10.0, 2.0), "dmod": (24.67, 0.07),
        "pmra": (31.0, 19.0), "pmdec": (-29.0, 16.0), "vr": (-179.2, 1.8),
    },
    "LMC": {
        "m": (15.0, 5.0), "c": (10.0, 2.0), "dmod": (18.477, 0.026),
        "pmra": (1910.0, 20.0), "pmdec": (229.0, 47.0), "vr": (262.2, 3.4),
    },
}

WU2026 = deepcopy(SAWALA2025)
WU2026["M31"].update(c=(15.0, 2.5), pmra=(45.9, 8.1), pmdec=(-20.5, 6.6))
WU2026["M33"].update(pmra=(30.6, 5.7), pmdec=(19.0, 5.7))

CONFIGS = {"sawala2025": SAWALA2025, "wu2026": WU2026}

SOURCES = {
    "sawala2025": "Sawala et al. (2025), Nature Astronomy 9, 1206; upstream notebook "
                  "commit e04820e, DefineParameters('fiducial')",
    "wu2026": "Wu et al. (2026), arXiv:2603.22863v1, Table 1",
}
