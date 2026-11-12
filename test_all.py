import numpy as np
import pytest

from gauge_qc import run_hourly_qc, repeated_value_flag
from radar_kriging import (variogram_model, pairwise_variogram, fit_all_models,
                     ordinary_kriging_ratio, kriging_external_drift, nrmse, relative_difference_pct)


def test_repeated_value_flag():
    s = np.array([0, 1.2, 1.2, 1.2, 1.2, 0, 0, 0, 0, 3.0, 3.0])
    f = repeated_value_flag(s, 4)
    assert f[1:5].all() and not f[5:9].any() and not f[9:].any()


@pytest.mark.parametrize("model", ["exponential", "spherical", "gaussian"])
def test_model_shape(model):
    p = (0.1, 0.5, 10_000)
    assert variogram_model(model, p, 0.0) == pytest.approx(0.1)
    assert variogram_model(model, p, 1e7) == pytest.approx(0.6, rel=1e-3)
    assert np.all(np.diff(variogram_model(model, p, np.linspace(0, 5e4, 50))) >= -1e-12)


def test_variogram_recovers_parameters():
    rng = np.random.default_rng(0)
    xy = rng.uniform(0, 50_000, (40, 2))
    d = np.hypot(xy[:, None, 0] - xy[None, :, 0], xy[:, None, 1] - xy[None, :, 1])
    cov = 0.5 * np.exp(-d / 10_000) + 0.05 * np.eye(40)
    z = rng.standard_normal((3000, 40)) @ np.linalg.cholesky(cov).T
    h, g, n = pairwise_variogram(z, xy)
    fits, best = fit_all_models(h, g, n)
    e = fits["exponential"]
    assert e["C0"] == pytest.approx(0.05, abs=0.03)
    assert e["CH"] == pytest.approx(0.5, rel=0.25)
    assert e["AH"] == pytest.approx(10_000, rel=0.35)


def test_kriging_exact_for_constant_field():
    xy = np.array([[0, 0], [1000, 0], [0, 1000], [1000, 1000.0]])
    p = (0.0, 1.0, 5000)
    lnr = np.full(4, np.log(1.5))
    assert ordinary_kriging_ratio(lnr, xy, np.array([500, 500.0]), 2.0, p) == pytest.approx(3.0)
    lns = np.log(np.array([1.0, 2.0, 3.0, 4.0]))
    lnr_d = lns.copy()                      # gauge == radar at neighbours
    val = kriging_external_drift(lns, lnr_d, xy, np.array([500, 500.0]), 2.5, p)
    assert val == pytest.approx(2.5, rel=1e-6)


def test_too_few_neighbours_returns_nan():
    xy = np.array([[0, 0], [1000, 0], [0, 1000.0]])
    lnr = np.array([0.1, np.nan, 0.2])
    assert np.isnan(ordinary_kriging_ratio(lnr, xy, np.array([10, 10.0]), 1.0, (0, 1, 5000)))


def test_qc_flags():
    rng = np.random.default_rng(1)
    xy = rng.uniform(0, 20_000, (10, 2))
    rain = np.zeros((200, 10))
    rain[50, :] = 2.0
    rain[50, 0] = 2.1
    rain[120, 4] = 30.0                     # isolated, all neighbours dry
    rain[10:15, 2] = 1.4                    # stuck
    rain[60, 6] = 500.0                     # > 200 mm/h
    r = run_hourly_qc(rain, xy)
    assert r.flag_repeat[10:15, 2].all()
    assert r.flag_high[60, 6]
    assert r.flag_all_neighbours_zero[120, 4]
    assert not r.valid3[120, 4] and r.valid3[50, 0]


def test_metrics():
    obs = np.array([1.0, 2.0, 3.0])
    assert nrmse(obs, obs) == 0
    assert relative_difference_pct(obs * 1.1, obs) == pytest.approx(10.0)
