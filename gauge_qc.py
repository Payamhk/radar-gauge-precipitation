"""Hourly rain-gauge quality control with a neighbour consistency check.

Translated from ``qc_hourly_consistency_check.m`` in this repository.

Progressive filtering on an hourly matrix ``rain`` (hours x stations):

1. ``flag_repeat``: runs of >= ``repeat_run_length`` equal, non-zero hourly
   values (stuck sensor).                      Valid1 = Total - flag_repeat
2. ``flag_high``: values above ``hourly_max_mm`` (default 200 mm/h).
                                               Valid2 = Valid1 - flag_high
3. Spatial consistency on Valid2. The target is compared with an inverse-
   distance-squared (IDW) estimate from its 4 nearest stations, using only
   neighbours that are themselves Valid2 at that hour.

   * fewer than ``min_valid_neighbours`` valid neighbours -> NotEvaluated
     (neutral status, not an error)
   * otherwise the target value is put in one of four WMO intensity classes:

     ========  ======================  =================================
     class     target rain P (mm/h)    metric
     ========  ======================  =================================
     NoRain    P == 0                  P_IDW, flagged if >= 3 mm
     Light     0 < P < 2.5             SARD = |2 (P - P_IDW) / (P + P_IDW)|
     Moderate  2.5 <= P < 10           SARD
     Heavy     P >= 10                 SARD
     ========  ======================  =================================

   * complementary rule: target >= 1 mm while every valid neighbour is 0.

Two final masks are returned:

* ``valid3``: Light/Moderate/Heavy are *not* flagged by SARD (original table)
* ``valid3_q``: Light/Moderate/Heavy flagged when SARD >= the class-specific
  quantile (default Q95) of its own pool.

Records with P > 0 but P_IDW == 0 are excluded from the SARD pools, because
SARD is then always exactly 2 regardless of P; those cases are handled by the
complementary all-neighbours-zero rule instead.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

CLASSES = ("NoRain", "Light", "Moderate", "Heavy")


@dataclass
class QCConfig:
    hourly_max_mm: float = 200.0
    repeat_run_length: int = 4
    min_valid_neighbours: int = 3
    n_neighbours: int = 4
    light_upper_mm: float = 2.5
    moderate_upper_mm: float = 10.0
    norain_fixed_threshold_mm: float = 3.0
    all_neighbours_zero_target_mm: float = 1.0
    sard_quantile: float = 0.95


@dataclass
class QCResult:
    flag_repeat: np.ndarray
    flag_high: np.ndarray
    valid1: np.ndarray
    valid2: np.ndarray
    not_evaluated: np.ndarray
    eligible: np.ndarray
    p_idw: np.ndarray
    n_valid_neighbours: np.ndarray
    rain_class: np.ndarray
    metric: np.ndarray
    flag_all_neighbours_zero: np.ndarray
    flag_not_aligned: np.ndarray
    flag_not_aligned_q: np.ndarray
    valid3: np.ndarray
    valid3_q: np.ndarray
    sard_thresholds: dict = field(default_factory=dict)

    def summary(self) -> dict:
        """Counts of flagged records, comparable to the MATLAB console log."""
        return {
            "repeat": int(self.flag_repeat.sum()),
            "above_max": int(self.flag_high.sum()),
            "not_evaluated": int(self.not_evaluated.sum()),
            "all_neighbours_zero": int(self.flag_all_neighbours_zero.sum()),
            "not_aligned": int(self.flag_not_aligned.sum()),
            "not_aligned_q": int(self.flag_not_aligned_q.sum()),
            "valid3": int(self.valid3.sum()),
            "valid3_q": int(self.valid3_q.sum()),
            **{f"sard_q_{k}": v for k, v in self.sard_thresholds.items()},
        }


def repeated_value_flag(series: np.ndarray, run_length: int) -> np.ndarray:
    """True where a non-zero value belongs to a run of >= run_length equal values."""
    s = np.asarray(series, dtype=float)
    same = np.r_[False, np.diff(s) == 0] & (s != 0)
    run_id = np.cumsum(~same)
    run_len = np.bincount(run_id)[run_id]
    return (run_len >= run_length) & (s != 0)


def nearest_neighbours(xy: np.ndarray, k: int) -> tuple[np.ndarray, np.ndarray]:
    """Indices and distances of the k nearest other stations, for each station."""
    xy = np.asarray(xy, dtype=float)
    d = np.hypot(xy[:, None, 0] - xy[None, :, 0], xy[:, None, 1] - xy[None, :, 1])
    np.fill_diagonal(d, np.inf)
    idx = np.argsort(d, axis=1)[:, :k]
    return idx, np.take_along_axis(d, idx, axis=1)


def run_hourly_qc(rain: np.ndarray, xy: np.ndarray, cfg: QCConfig | None = None) -> QCResult:
    """Run the three-step QC on an hourly matrix.

    Parameters
    ----------
    rain : (n_hours, n_stations) hourly precipitation in mm. Use 0 for hours
        that were artificially filled; they are treated as ordinary values,
        as in the MATLAB version.
    xy : (n_stations, 2) projected station coordinates in metres.
    """
    cfg = cfg or QCConfig()
    rain = np.asarray(rain, dtype=float)
    n_h, n_s = rain.shape

    flag_repeat = np.column_stack(
        [repeated_value_flag(rain[:, i], cfg.repeat_run_length) for i in range(n_s)]
    )
    valid1 = ~flag_repeat
    flag_high = valid1 & (rain > cfg.hourly_max_mm)
    valid2 = valid1 & ~flag_high

    nb_idx, nb_dist = nearest_neighbours(xy, cfg.n_neighbours)

    p_idw = np.full((n_h, n_s), np.nan)
    n_valid = np.zeros((n_h, n_s), dtype=int)
    eligible = np.zeros((n_h, n_s), dtype=bool)
    not_eval = np.zeros((n_h, n_s), dtype=bool)
    rain_class = np.full((n_h, n_s), "", dtype=object)
    metric = np.full((n_h, n_s), np.nan)
    flag_anz = np.zeros((n_h, n_s), dtype=bool)

    for i in range(n_s):
        cols = nb_idx[i]
        w = 1.0 / nb_dist[i] ** 2
        nb_rain = rain[:, cols]
        nb_ok = valid2[:, cols]
        nv = nb_ok.sum(axis=1)
        w_ok = w[None, :] * nb_ok
        sw = w_ok.sum(axis=1)
        pidw = (w_ok * nb_rain).sum(axis=1) / np.maximum(sw, np.finfo(float).eps)

        p = rain[:, i]
        elig = valid2[:, i] & (nv >= cfg.min_valid_neighbours)
        n_valid[:, i] = nv
        p_idw[:, i] = pidw
        eligible[:, i] = elig
        not_eval[:, i] = valid2[:, i] & ~elig

        m_no = elig & (p == 0)
        m_li = elig & (p > 0) & (p < cfg.light_upper_mm)
        m_mo = elig & (p >= cfg.light_upper_mm) & (p < cfg.moderate_upper_mm)
        m_he = elig & (p >= cfg.moderate_upper_mm)
        for name, m in zip(CLASSES, (m_no, m_li, m_mo, m_he)):
            rain_class[m, i] = name

        metric[m_no, i] = pidw[m_no]
        sard_ok = (m_li | m_mo | m_he) & (pidw != 0)
        metric[sard_ok, i] = np.abs(
            2 * (p[sard_ok] - pidw[sard_ok]) / np.maximum(p[sard_ok] + pidw[sard_ok], np.finfo(float).eps)
        )

        nb_nonzero = (nb_ok & (nb_rain > 0)).sum(axis=1)
        flag_anz[:, i] = elig & (p >= cfg.all_neighbours_zero_target_mm) & (nb_nonzero == 0)

    masks = {c: rain_class == c for c in CLASSES}

    # Original table: only the fixed No-Rain rule (+ complementary rule)
    flag_na = np.zeros((n_h, n_s), dtype=bool)
    flag_na[masks["NoRain"]] = metric[masks["NoRain"]] >= cfg.norain_fixed_threshold_mm
    flag_na = (flag_na & eligible) | flag_anz

    # Quantile table: SARD >= class quantile for Light/Moderate/Heavy
    thr = {}
    flag_q = np.zeros((n_h, n_s), dtype=bool)
    flag_q[masks["NoRain"]] = metric[masks["NoRain"]] >= cfg.norain_fixed_threshold_mm
    for c in CLASSES[1:]:
        pool = metric[masks[c]]
        pool = pool[~np.isnan(pool)]
        # 'hazen' reproduces MATLAB's quantile() definition
        thr[c] = float(np.quantile(pool, cfg.sard_quantile, method="hazen")) if pool.size else np.nan
        vals = metric[masks[c]]
        flag_q[masks[c]] = np.where(np.isnan(vals), False, vals >= thr[c])
    flag_q = (flag_q & eligible) | flag_anz

    valid3 = valid2 & ~not_eval & ~flag_na
    valid3_q = valid2 & ~not_eval & ~flag_q

    return QCResult(
        flag_repeat=flag_repeat, flag_high=flag_high, valid1=valid1, valid2=valid2,
        not_evaluated=not_eval, eligible=eligible, p_idw=p_idw, n_valid_neighbours=n_valid,
        rain_class=rain_class, metric=metric, flag_all_neighbours_zero=flag_anz,
        flag_not_aligned=flag_na, flag_not_aligned_q=flag_q, valid3=valid3, valid3_q=valid3_q,
        sard_thresholds=thr,
    )
