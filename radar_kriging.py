"""Radar-gauge merging: pairwise variogram, ordinary kriging and KED.

Idea in one sentence: the radar error is mostly *multiplicative* and
*spatially correlated*, so we interpolate the error itself between gauges
and use it to correct the radar where there is no gauge.

Steps
-----
1. ``pairwise_variogram``: one variogram point per station pair (no distance
   binning). With only ~40 gauges, binning would hide the short-range
   structure that matters most for convective rain.
2. ``fit_all_models``: weighted least squares (weights = sqrt of shared
   hours, Cressie 1985) for exponential, spherical and Gaussian models;
   the lowest weighted error wins. Parameters stuck on a bound are flagged.
3. ``ordinary_kriging_ratio``: krige z = ln(S/R) and return R * exp(z_hat).
   Working on the log-ratio turns the multiplicative error into an additive
   one and stabilises its variance across intensities.
4. ``kriging_external_drift``: krige ln(S) using ln(R) as drift, with two
   safeguards (no drift extrapolation, cap on |weights|) and fallback to
   the raw radar when the system is unreliable.
5. ``nrmse``, ``pearson_r``, ``relative_difference_pct``: verification.

Translated from the MATLAB scripts in this repository
(``variogram_fit_seasonal.m``, ``ked_correction.m``).
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import least_squares


# ---------------------------------------------------------------- variogram
MODELS = ("exponential", "spherical", "gaussian")
LOWER = np.array([0.0, 0.05, 500.0])        # C0, CH, AH (m)
UPPER = np.array([0.5, 4.0, 200_000.0])


def variogram_model(name: str, params, h):
    """gamma(h) for nugget C0, partial sill CH and range parameter AH."""
    c0, ch, ah = params
    h = np.asarray(h, dtype=float)
    r = h / ah
    if name == "exponential":
        return c0 + ch * (1 - np.exp(-r))
    if name == "spherical":
        return c0 + ch * np.where(h <= ah, 1.5 * r - 0.5 * r**3, 1.0)
    if name == "gaussian":
        return c0 + ch * (1 - np.exp(-(r**2)))
    raise ValueError(f"unknown model {name!r}")


def pairwise_variogram(z: np.ndarray, xy: np.ndarray, time_mask=None, min_shared: int = 10):
    """Empirical variogram, one point per station pair.

    Parameters
    ----------
    z : (n_hours, n_stations) values on a common hourly grid, NaN where missing
        (e.g. not wet-wet).
    xy : (n_stations, 2) coordinates in metres.
    time_mask : optional boolean (n_hours,) to restrict to a season.
    min_shared : pairs sharing fewer valid hours are skipped.

    Returns
    -------
    h, gamma, n : arrays over the retained pairs.
    """
    z = np.asarray(z, dtype=float)
    n_st = z.shape[1]
    if time_mask is not None:
        z = z[np.asarray(time_mask, dtype=bool)]
    ok = ~np.isnan(z)
    h, g, n = [], [], []
    for i in range(n_st - 1):
        for j in range(i + 1, n_st):
            m = ok[:, i] & ok[:, j]
            k = int(m.sum())
            if k < min_shared:
                continue
            d = z[m, i] - z[m, j]
            h.append(np.hypot(*(xy[i] - xy[j])))
            g.append(np.sum(d**2) / (2 * k))
            n.append(k)
    return np.array(h), np.array(g), np.array(n)


def initial_guess(h, g):
    """Starting values read off the empirical cloud (Chiles & Delfiner, 2012).

    C0: mean gamma below the 10th percentile of distances; total sill: 85th
    percentile of gamma; AH: distance where gamma reaches C0 + 0.63 CH on a
    temporary 50-bin profile up to the 95th distance percentile.
    """
    pct = lambda a, q: np.percentile(a, q, method="hazen")  # MATLAB prctile
    short = h <= pct(h, 10)
    c0 = max(np.mean(g[short]), 0.005)
    sill = max(pct(g, 85), c0 + 0.05)
    ch = sill - c0
    target = c0 + 0.63 * ch
    edges = np.linspace(0, pct(h, 95), 51)
    bh, bg = [], []
    for a, b in zip(edges[:-1], edges[1:]):
        m = (h >= a) & (h < b)
        if m.any():
            bh.append(h[m].mean())
            bg.append(g[m].mean())
    ah = max(bh[int(np.argmin(np.abs(np.array(bg) - target)))], 500.0)
    p0 = np.array([c0, ch, ah])
    return np.clip(p0, LOWER + 1e-6, UPPER - 1e-6)


def fit_variogram(name, h, g, n, p0=None, lower=LOWER, upper=UPPER):
    """Weighted least squares fit with weights sqrt(n). Returns (params, wsse)."""
    p0 = initial_guess(h, g) if p0 is None else np.asarray(p0, dtype=float)
    w = np.sqrt(n)
    res = least_squares(
        lambda p: w * (variogram_model(name, p, h) - g),
        p0, bounds=(lower, upper), xtol=1e-10, ftol=1e-10, gtol=1e-10, max_nfev=5000,
    )
    return res.x, float(np.sum(res.fun**2))


def fit_all_models(h, g, n, upper=UPPER):
    """Fit the three models; returns dict and the name of the lowest-WSSE model.

    A fitted parameter that ends up on its bound (e.g. the range in a season
    with too few pairs) is reported in ``at_bound`` so it can be treated as
    unreliable rather than silently used.
    """
    p0 = initial_guess(h, g)
    out = {}
    for name in MODELS:
        p, wsse = fit_variogram(name, h, g, n, p0, upper=upper)
        at_bound = bool(np.any(np.isclose(p, LOWER, rtol=1e-3) | np.isclose(p, upper, rtol=1e-3)))
        out[name] = {"C0": p[0], "CH": p[1], "AH": p[2], "wsse": wsse, "at_bound": at_bound}
    best = min(out, key=lambda k: out[k]["wsse"])
    return out, best


# ---------------------------------------------------------------- kriging
def _cov_matrix(model, params, xy_nb):
    c0, ch, _ = params
    sill = c0 + ch
    d = np.hypot(xy_nb[:, None, 0] - xy_nb[None, :, 0], xy_nb[:, None, 1] - xy_nb[None, :, 1])
    c = sill - variogram_model(model, params, d)
    c[d == 0] = sill
    return c


def _neighbours(values_t, xy_train, xy_target, radius):
    ok = ~np.isnan(values_t)
    if np.isfinite(radius):
        ok &= np.hypot(*(xy_train - xy_target).T) <= radius
    return np.flatnonzero(ok)


def ordinary_kriging_ratio(lnr_t, xy_train, xy_target, radar_target, params,
                           model="exponential", radius=np.inf, min_stations=3):
    """OK of ln(S/R) at one hour. Returns the corrected rainfall or NaN.

    lnr_t : (n_train,) ln(S/R) of the training stations at this hour (NaN if
        not wet-wet).
    """
    idx = _neighbours(lnr_t, xy_train, xy_target, radius)
    n = idx.size
    if n < min_stations:
        return np.nan
    c0, ch, _ = params
    a = np.zeros((n + 1, n + 1))
    a[:n, :n] = _cov_matrix(model, params, xy_train[idx])
    a[n, :n] = a[:n, n] = 1.0
    b = np.ones(n + 1)
    hi = np.hypot(*(xy_train[idx] - xy_target).T)
    b[:n] = (c0 + ch) - variogram_model(model, params, hi)
    if 1.0 / np.linalg.cond(a, 1) < 1e-12:
        return np.nan
    lam = np.linalg.solve(a, b)
    return float(radar_target * np.exp(lam[:n] @ lnr_t[idx]))


def kriging_external_drift(lns_t, lnr_train_t, xy_train, xy_target, radar_target, params,
                           model="exponential", radius=np.inf, min_stations=3,
                           safeguard_extrapolation=True, max_abs_weight=10.0,
                           return_diagnostics=False):
    """KED of ln(S) with ln(R) as external drift at one hour.

    Safeguards (both optional, as in the MATLAB scripts):

    1. ``safeguard_extrapolation``: reject the hour when ln(R) at the target is
       outside the range of ln(R) at the neighbours (drift extrapolation).
    2. ``max_abs_weight``: reject when max |lambda| exceeds this value (set to
       None to disable). As in the original code, the maximum is taken over
       the full solution vector, including the two Lagrange multipliers.
    """
    diag = {"extrapolated": False, "max_abs_lambda": np.nan}
    nan = (np.nan, diag) if return_diagnostics else np.nan
    idx = _neighbours(lns_t, xy_train, xy_target, radius)
    n = idx.size
    if n < min_stations:
        return nan
    c0, ch, _ = params
    drift = lnr_train_t[idx]
    t_drift = np.log(radar_target)
    diag["extrapolated"] = bool(t_drift < drift.min() or t_drift > drift.max())
    if safeguard_extrapolation and diag["extrapolated"]:
        return nan
    a = np.zeros((n + 2, n + 2))
    a[:n, :n] = _cov_matrix(model, params, xy_train[idx])
    a[n, :n] = a[:n, n] = 1.0
    a[n + 1, :n] = a[:n, n + 1] = drift
    b = np.zeros(n + 2)
    hi = np.hypot(*(xy_train[idx] - xy_target).T)
    b[:n] = (c0 + ch) - variogram_model(model, params, hi)
    b[n], b[n + 1] = 1.0, t_drift
    if 1.0 / np.linalg.cond(a, 1) < 1e-12:
        return nan
    lam = np.linalg.solve(a, b)
    diag["max_abs_lambda"] = float(np.max(np.abs(lam)))
    if max_abs_weight is not None and diag["max_abs_lambda"] > max_abs_weight:
        return nan
    val = float(np.exp(lam[:n] @ lns_t[idx]))
    return (val, diag) if return_diagnostics else val


# ---------------------------------------------------------------- metrics
def _clean(pred, obs):
    pred, obs = np.asarray(pred, float), np.asarray(obs, float)
    m = ~np.isnan(pred) & ~np.isnan(obs)
    return pred[m], obs[m]


def nrmse(pred, obs):
    """RMSE normalised by the mean observation."""
    p, o = _clean(pred, obs)
    return float(np.sqrt(np.mean((p - o) ** 2)) / np.mean(o))


def pearson_r(pred, obs):
    p, o = _clean(pred, obs)
    return float(np.corrcoef(p, o)[0, 1])


def relative_difference_pct(pred, obs):
    """Relative difference of totals, (sum(pred) - sum(obs)) / sum(obs) * 100."""
    p, o = _clean(pred, obs)
    return float((p.sum() - o.sum()) / o.sum() * 100)
