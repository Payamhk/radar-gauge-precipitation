# radar-gauge-precipitation

Code from my PhD project at the University of Brescia, carried out in collaboration with MeteoSwiss: *Accuracy Enhancement of Radar Rainfall Data and Predictive Flood Modelling*. It covers two steps:

1. **Quality control of hourly rain-gauge data**: removing bad gauge records before they are used.
2. **Correcting radar precipitation with gauges**: radar–gauge merging by ordinary kriging and kriging with external drift.

Study area: Franciacorta, Northern Italy. Data: 44 Consorzio gauges, 11 ARPA Lombardia gauges and the MeteoSwiss CombiPrecip product, hourly, March 2025 – March 2026. Data are not included; `demo.py` runs on synthetic data.

## Files

| File | What it does |
|---|---|
| `gauge_qc.py` | Hourly gauge QC (Python) |
| `radar_kriging.py` | Variogram fitting, ordinary kriging of ln(S/R), KED, verification scores (Python) |
| `demo.py` | End-to-end example on synthetic data |
| `test_all.py` | Unit tests (`pytest`) |
| `qc_hourly_consistency_check.m` | Original MATLAB version of the QC |
| `variogram_fit_seasonal.m` | Original MATLAB seasonal variogram fitting |
| `ked_correction.m` | Original MATLAB KED correction with validation |

```bash
pip install -r requirements.txt
python demo.py
pytest -q
```

## Key ideas

**1. Gauge QC by rain intensity, not one threshold** (`gauge_qc.py`)
- Each gauge is compared with an inverse-distance estimate from its 4 nearest neighbours. Neighbours are used only if they passed the earlier checks.
- The comparison depends on the rain class (WMO): no rain, light, moderate, heavy. For dry hours, the gauge is flagged if the neighbours give ≥ 3 mm. For wet hours, the symmetric relative difference SARD = |2(P − P_IDW)/(P + P_IDW)| is computed, and the flag threshold is the 95th percentile *of that class*. Errors of 1 mm mean different things in light and heavy rain.
- One special case is handled separately. If all neighbours are dry, SARD is always exactly 2 and carries no information. So the gauge is flagged if it reports ≥ 1 mm while every neighbour reports zero.
- If fewer than 3 neighbours are valid, the hour is marked *not evaluated* rather than wrong.

**2. Correct the radar error, not the rainfall** (`radar_kriging.py`)
- The radar error is mostly multiplicative, so the variable interpolated between gauges is z = ln(S/R) (gauge over radar).
- The log turns a multiplicative error into an additive one and keeps its variance stable across intensities. I checked this against S−R, S/R and ln(S) before choosing it.
- The corrected radar is R × exp(ẑ).

**3. Variogram from station pairs, without binning**
- With only ~40 gauges, distance bins average away the short-range structure, which matters most for convective rain.
- Each station pair is therefore one point: γ = mean of squared differences over the hours both gauges were wet.
- Pairs are weighted by √(shared hours) in a weighted least-squares fit (Cressie, 1985). Exponential, spherical and Gaussian models are compared, and the lowest weighted error wins.
- Fitted parameters that stick to their bounds (e.g. winter, with few events) are flagged as unreliable instead of being used silently.

**4. KED with safeguards**
- Kriging with external drift interpolates ln(S) using the radar ln(R) as a trend.
- It can extrapolate badly, so an hour is rejected if:
  - the radar value at the target is outside the neighbours' radar range, or
  - any kriging weight is larger than 10 in absolute value.
- Rejected hours fall back to the raw radar.

**5. Fair validation**
- Variograms are fitted only on training gauges. Correction is tested on gauges that were never used: 6 Consorzio and 11 ARPA.
- Scores are NRMSE, Pearson r and the relative difference of totals.

On the synthetic demo, ordinary kriging of ln(S/R) reduces NRMSE from 0.73 (raw radar) to 0.65, and r rises from 0.76 to 0.82.

## Related publication

Heidarian, P., Pilotti, M., Khoshkonesh, A., Nazari, R., Nikoo, M. R., Gabella, M., & Okhravi, S. *Adaptive selection of spatial rainfall interpolation methods across temporal scales: a multi-event analysis in Franciacorta Region, Northern Italy.*

## License

MIT
