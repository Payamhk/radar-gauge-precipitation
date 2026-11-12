"""End-to-end demo on synthetic data (no real data needed).

1. Simulate hourly gauge rainfall at 30 stations with a few planted errors
   and run the QC.
2. Simulate a multiplicative radar error field, fit the variogram of
   ln(S/R), and correct the radar at held-out stations with OK and KED.
"""
import numpy as np

from gauge_qc import QCConfig, run_hourly_qc
from radar_kriging import (pairwise_variogram, fit_all_models,
                           ordinary_kriging_ratio, kriging_external_drift, nrmse, pearson_r,
                           relative_difference_pct)

rng = np.random.default_rng(42)
n_st, n_h = 30, 2000
xy = rng.uniform(0, 40_000, size=(n_st, 2))


def gaussian_field(corr_len, n_samples):
    d = np.hypot(xy[:, None, 0] - xy[None, :, 0], xy[:, None, 1] - xy[None, :, 1])
    cov = np.exp(-d / corr_len) + 1e-8 * np.eye(n_st)
    return rng.standard_normal((n_samples, n_st)) @ np.linalg.cholesky(cov).T


# ---------- 1. gauge QC ----------
wet = gaussian_field(15_000, n_h) > 0.8
rain = np.where(wet, np.exp(0.8 * gaussian_field(15_000, n_h)), 0.0).round(1)
rain[100:106, 3] = 2.4          # stuck sensor
rain[500, 7] = 350.0            # impossible value
rain[800, 11] = 25.0            # isolated spike, dry neighbours
res = run_hourly_qc(rain, xy, QCConfig())
print("QC summary:", res.summary())
assert res.flag_repeat[100:106, 3].all() and res.flag_high[500, 7]

# ---------- 2. radar-gauge merging ----------
truth = np.where(wet, np.exp(0.8 * gaussian_field(15_000, n_h)) + 0.1, np.nan)
log_err = 0.5 * gaussian_field(12_000, n_h) - 0.2          # spatially correlated
radar = truth * np.exp(log_err)
test = np.arange(n_st) % 5 == 0                              # 6 held-out stations
train = ~test

lnr = np.log(truth / radar)
h, g, n = pairwise_variogram(lnr[:, train], xy[train])
fits, best = fit_all_models(h, g, n)
p = [fits[best][k] for k in ("C0", "CH", "AH")]
print(f"best variogram: {best}  C0={p[0]:.3f} CH={p[1]:.3f} AH={p[2]/1000:.1f} km")

obs, raw, ok_pred, ked_pred = [], [], [], []
for t in range(n_h):
    for j in np.flatnonzero(test):
        if np.isnan(truth[t, j]):
            continue
        o = ordinary_kriging_ratio(lnr[t, train], xy[train], xy[j], radar[t, j], p, best)
        k = kriging_external_drift(np.log(truth[t, train]), np.log(radar[t, train]), xy[train],
                                   xy[j], radar[t, j], p, best)
        obs.append(truth[t, j]); raw.append(radar[t, j])
        ok_pred.append(radar[t, j] if np.isnan(o) else o)      # fallback to raw radar
        ked_pred.append(radar[t, j] if np.isnan(k) else k)

for name, pred in (("radar", raw), ("OK", ok_pred), ("KED", ked_pred)):
    print(f"{name:6s} NRMSE={nrmse(pred, obs):.3f}  r={pearson_r(pred, obs):.3f}  "
          f"RD={relative_difference_pct(pred, obs):+.1f}%")
