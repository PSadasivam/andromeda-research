"""Vectorized port of the Sawala et al. (2025) orbit integrator (upstream commit e04820e).

Integrates all Monte Carlo samples simultaneously. Units follow upstream: kpc, km/s,
1e10 Msun, G = 43007, so one time unit is 1 kpc/(km/s) = 0.9778 Gyr (labelled "Gyr" upstream).
"""

import numpy as np
from colossus.cosmology import cosmology
from colossus.utils import constants
from scipy.special import erf

G = 43007.0
TIME_UNIT_GYR = 0.977792221
NFW_XMAX = 2.162581587064612
PERICENTRE_REBOUND_KPC = 20.0


def g_nfw(x):
    return np.log1p(x) - x / (1.0 + x)


def nfw_halo_properties(mass_1e10, concentration):
    """Scale radius (kpc), r200 (kpc), and Vmax (km/s) for M200c at z = 0 in Planck18."""
    cosmo = cosmology.setCosmology("planck18")
    h = cosmo.h
    m_h = mass_1e10 * 1e10 * h
    r200_h = (3.0 * m_h / (4.0 * np.pi * 200.0 * cosmo.rho_c(0.0))) ** (1.0 / 3.0)
    rs_h = r200_h / concentration
    rmax_h = NFW_XMAX * rs_h
    vmax = np.sqrt(constants.G * m_h * g_nfw(NFW_XMAX) / g_nfw(concentration) / rmax_h)
    r_s = rs_h / h
    return r_s, r_s * concentration, vmax


def accelerations(pos, mass, softening):
    dx = pos[:, None, :, :] - pos[:, :, None, :]
    inv_r3 = (np.sum(dx**2, axis=-1) + softening**2) ** -1.5
    return G * np.einsum("sijk,sij,sj->sik", dx, inv_r3, mass)


def dynamical_friction(pos, vel, mass, conc, r_s, r_200, vmax, scheme):
    fric = np.zeros_like(pos)
    if scheme == "none":
        return fric
    n_body = mass.shape[1]
    for i in range(n_body):
        for j in range(n_body):
            if i == j:
                continue
            active = (mass[:, i] > 0) & (mass[:, j] > mass[:, i])
            r = np.maximum(10.0, np.linalg.norm(pos[:, i] - pos[:, j], axis=-1))
            active &= r < r_200[:, j]
            if not active.any():
                continue
            dv = vel[:, i] - vel[:, j]
            v = np.linalg.norm(dv, axis=-1)
            rho = mass[:, j] / (4.0 * np.pi * g_nfw(conc[:, j]) * r * (r + r_s[:, j]) ** 2)
            x = r / r_s[:, j]
            sigma_v = vmax[:, j] * 1.4393 * x**0.354 / (1.0 + 1.1756 * x**0.725)
            y = v / (np.sqrt(2.0) * sigma_v)
            integral = erf(y) - 2.0 * y / np.sqrt(np.pi) * np.exp(-y * y)
            epsilon = np.where(r_s[:, i] > 8.0, 2.2 * r_s[:, i] - 14.0, 0.45 * r_s[:, i])
            fricmag = (4.0 * np.pi * G**2 * np.log(r / epsilon) * mass[:, i] * rho * integral
                       / np.maximum(v, 5.0) ** 3)
            fricmag = np.where(active, fricmag, 0.0)[:, None] * dv
            if scheme == "proportional":
                total = mass[:, i] + mass[:, j]
                share_i = np.divide(mass[:, j], total, out=np.zeros_like(total), where=total > 0)
                share_j = np.divide(mass[:, i], total, out=np.zeros_like(total), where=total > 0)
                fric[:, i] -= fricmag * share_i[:, None]
                fric[:, j] += fricmag * share_j[:, None]
            elif scheme == "hierarchical":
                fric[:, i] -= fricmag
            else:
                raise ValueError(f"unknown dynamical friction scheme: {scheme}")
    return fric


def integrate(ics, t_end=10.0, dt=0.001, softening=20.0, merger_threshold=20.0,
              friction="proportional", track_every=10, track_samples=200):
    """Integrate all samples; returns per-sample merger and pericentre summaries."""
    names = ics["names"]
    pos = ics["positions"].copy()
    vel = ics["velocities"].copy()
    mass = ics["masses"].copy()
    conc = ics["concentrations"]
    n_s, n_body = mass.shape

    r_s, r_200, vmax = nfw_halo_properties(mass, conc)

    acc = accelerations(pos, mass, softening)
    fric = dynamical_friction(pos, vel, mass, conc, r_s, r_200, vmax, friction)
    vel -= (np.sum(mass[:, :, None] * vel, axis=1) / np.sum(mass, axis=1)[:, None])[:, None, :]

    alive = np.ones((n_s, n_body), dtype=bool)
    absorbed_step = np.full((n_s, n_body), -1)
    merge_step = {pair: np.full(n_s, -1) for pair in _tracked_pairs(names)}
    min_sep = {(a, b): np.linalg.norm(pos[:, a] - pos[:, b], axis=-1) for a, b in merge_step}

    d_start = np.linalg.norm(pos[:, 1] - pos[:, 0], axis=-1)
    run_min = d_start.copy()
    run_min_step = np.full(n_s, -1)
    peri_step = np.full(n_s, -1)
    peri_dist = np.full(n_s, np.nan)

    n_steps = int(np.ceil(t_end / dt))
    n_tracked = min(track_samples, n_s)
    track = np.full((n_tracked, n_steps // track_every + 1), np.nan, dtype=np.float32)
    track[:, 0] = d_start[:n_tracked]

    for step in range(n_steps):
        vel += (acc + fric) * dt / 2.0
        pos += vel * dt
        acc = accelerations(pos, mass, softening)
        fric = dynamical_friction(pos, vel, mass, conc, r_s, r_200, vmax, friction)
        vel += (acc + fric) * dt / 2.0

        for l in range(n_body):
            for m in range(l):
                close = (alive[:, l] & alive[:, m]
                         & (np.linalg.norm(pos[:, l] - pos[:, m], axis=-1) < merger_threshold))
                if not close.any():
                    continue
                keep_l = mass[:, l] > mass[:, m]
                keep = np.where(keep_l, l, m)[close]
                remove = np.where(keep_l, m, l)[close]
                idx = np.flatnonzero(close)
                m_k, m_r = mass[idx, keep], mass[idx, remove]
                total = m_k + m_r
                pos[idx, keep] = (m_k[:, None] * pos[idx, keep] + m_r[:, None] * pos[idx, remove]) / total[:, None]
                vel[idx, keep] = (m_k[:, None] * vel[idx, keep] + m_r[:, None] * vel[idx, remove]) / total[:, None]
                mass[idx, keep] = total
                mass[idx, remove] = 0.0
                alive[idx, remove] = False
                absorbed_step[idx, remove] = step

        for (a, b), steps in merge_step.items():
            gone = ~(alive[:, a] & alive[:, b])
            newly = gone & (steps < 0)
            steps[newly] = step
            sep = np.linalg.norm(pos[:, a] - pos[:, b], axis=-1)
            live = ~gone
            min_sep[a, b][live] = np.minimum(min_sep[a, b][live], sep[live])

        mw_m31_live = alive[:, 0] & alive[:, 1]
        d_now = np.linalg.norm(pos[:, 1] - pos[:, 0], axis=-1)
        pending = peri_step < 0
        lower = pending & mw_m31_live & (d_now < run_min)
        run_min[lower] = d_now[lower]
        run_min_step[lower] = step
        # A rebound threshold ignores the few-kpc jumps when a satellite merges into its host.
        rebound = pending & mw_m31_live & (d_now > run_min + PERICENTRE_REBOUND_KPC)
        peri_step[rebound] = run_min_step[rebound]
        peri_dist[rebound] = run_min[rebound]
        merged_first = pending & ~mw_m31_live
        peri_step[merged_first] = step
        peri_dist[merged_first] = merger_threshold

        if (step + 1) % track_every == 0:
            track[:, (step + 1) // track_every] = np.where(mw_m31_live[:n_tracked], d_now[:n_tracked], np.nan)

    # Upstream reports merger time as t_all[i] for a merger detected in loop step i.
    to_time = lambda s: np.where(s >= 0, s * dt, np.nan)
    result = {
        "names": names, "t_end": t_end, "dt": dt, "softening": softening,
        "merger_threshold": merger_threshold, "friction": friction,
        "first_pericentre_time": to_time(peri_step), "first_pericentre_distance": peri_dist,
        "track_dt": dt * track_every, "mw_m31_separation_track": track,
    }
    for (a, b), steps in merge_step.items():
        key = f"{names[a]}-{names[b]}"
        result[f"merge_time_{key}"] = to_time(steps)
        result[f"min_separation_{key}"] = min_sep[a, b]
    # Pair "merge" times follow upstream (either galaxy gone); absorption times record each galaxy's own fate.
    for k, name in enumerate(names):
        result[f"absorbed_time_{name}"] = to_time(absorbed_step[:, k])
    return result


def _tracked_pairs(names):
    pairs = [(0, 1)]
    if "LMC" in names:
        pairs.append((names.index("LMC"), 0))
    if "M33" in names:
        pairs.append((names.index("M33"), 1))
    return pairs
