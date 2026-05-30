# WDM wavelet whitening normalization — fix plan

## Symptom

Running `test_mojito_noise.sh` (`noise_wavelet_mcmc --CD1L` on a MOJITO file), the data
is not whitened: GLASS's own `data/whitened_data.dat` had **variance ~4e7** (should be
~1). Equivalently, the data wavelet power `|w|^2` is ~1e8x the noise model.

## Diagnosis (confirmed numerically, 2026-05-29)

Two compounding normalization bugs in the wavelet data path:

1. **Hard-coded cadence.** `LISA_CADENCE = 5` (`utils/src/glass_lisa.h:34`) drives
   `wdm->NF` and the WDM cadence (`utils/src/glass_wavelet.c:214`);
   `initialize_wavelet(data->wdm, data->T)` (`utils/src/glass_data.c:277`) only gets
   `T`, never the file's `dt`. CD1L files are 2.5 s (or 0.25 s), so the data length is
   2x (or 20x) the WDM grid `NF*NT`, breaking the FFT->WDM normalization.
   Worth ~80x here.

2. **Unnormalized data FFT.** `glass_forward_real_fft_outplace`
   (`utils/src/glass_math.c:808`) is raw `kiss_fftr` (no `dt`, no `1/N`). The model
   side (`generate_full_dynamic_covariance_matrix`) populates `noise->C` with the raw
   physical PSD `S(f)` and never applies the `/ND` conversion. Net: data `|w|^2` too
   big by `~NF*NT/dt^2` (~5e5 here).

Quick-test ladder (variance):

| stage | whitened variance |
|---|---|
| original (2.5 s data) | ~4e7 |
| + cadence hack (read every other point -> 5 s) | ~5e5 |
| + FFT-norm scaling (`dt/sqrt(Nout)`) | ~1.4 |

The residual **~1.4 ~ sqrt(2)** is a one-sided/two-sided (real-FFT) convention factor
(see below).

## The normalization chain

Real series `x_i`, `N` samples, cadence `dt`, `T = N*dt`. Model encodes the physical
**one-sided** PSD `S(f)` (`XYZnoise_FF`, units 1/Hz), `<x^2> = integral_0^fNy S df`.

- **Forward FFT (kiss, unnormalized):** `X_k = sum_i x_i exp(-2pi i k i / N)`.
  Physical FFT is `Dt * X_k`. One-sided periodogram: `S(f_k) = (2 dt / N) |X_k|^2`, so
  `<|X_k|^2> = (N / 2dt) * S(f_k)` for `0 < k < N/2`. The **factor 2** is the one-sided
  fold; DC (`k=0`) and Nyquist (`k=N/2`) have no fold -> no factor 2.
- **Wavelet transform:** windows `X_k` with the Meyer window, inverse-FFTs each layer
  (÷`Nt`). Unit test validates `<|w_nm|^2> = <|X_k|^2> / ND`, `ND = NF*NT = N`.
  Therefore the **correctly normalized per-pixel noise variance is**
  `<|w_nm|^2> = S(f) / (2 dt)`.

### Why the sqrt(2)

`S(f)` in the model is **one-sided** (carries the x2 fold), but a wavelet pixel is built
from an inverse **complex** FFT of windowed bins — a **two-sided** object. Mixing a
one-sided PSD with a two-sided transform leaves a factor 2 in power = sqrt(2) in
whitened amplitude. The code already does this bookkeeping at the band edges in
`wavelet_transform_freq`:

```c
if (m==0)       wdmdata[n]   = DX_t[2*n] * M_SQRT2;   // DC layer
else if (m==Nf) wdmdata[n+1] = DX_t[2*n] * M_SQRT2;   // Nyquist layer
...
if ((j==0) && (m==0 || m==Nf)) weight /= 2.0;         // "less d.o.f. in boundary layers"
```

The interior layers and the one-sided `S(f)` model aren't reconciled with that same
convention -> bulk pixels off by this factor.

## Fix (step 2 — the real fix)

Pick ONE convention and enforce it on both sides:

1. **Normalize the forward FFT by `dt`** in the wavelet data path so it is the physical
   FFT that `S(f)` is defined against (alternatively carry the factor analytically in
   the model — but normalizing the FFT is cleaner).
2. **Build the wavelet noise model as `S(f)/(2 dt)` per pixel** — i.e. route the model
   through `stationary_dft_psd_to_wdm_psd` (which encodes `<|FFT|^2>/ND` and the same
   edge-layer `M_SQRT2` handling) instead of dumping raw `S(f)` into `noise->C`. The
   production builder `generate_full_dynamic_covariance_matrix` should call this
   conversion.

Both the `NF*NT/dt^2` and the factor-2 then fall out by construction; no magic
constants.

## Fix (step 0 — cadence, removes the decimation hack)

Derive the WDM cadence from `tdi->delta` (the file's `dt`) instead of the compile-time
`LISA_CADENCE`. Set `NF`/cadence from the data so any cadence works; then the every-
other-point decimation is unnecessary.

## Temporary hacks to remove

Both live in `LISA_Read_HDF5_CD1L_TDI` (`utils/src/glass_lisa.c`), tagged:
- `TODO(remove-cadence-hack)` — decimate by 2 (read every other point), `dt *= 2`.
- `TODO(remove-fft-norm-hack)` — scale data by `dt/sqrt(Nout)`.

## Verification

- Re-run `test_mojito_noise.sh`; confirm whitened variance ~= 1.0 (not 1.4).
- Print the **per-layer** `<|w|^2> / model` ratio and confirm it sits flat at 1.0
  across the band (not just the aggregate variance) — this distinguishes a clean
  factor-2 from a frequency-dependent residual.
- Add a whitening sanity assertion (var ~ 1 within tolerance) so this can't regress
  silently.
- Note: `full_noise_model.dat` is raw physical PSD `S(f)` (= `instrument_noise_model.dat`
  replicated over time, no WDM normalization). Do NOT whiten against it; use the
  internal model / `whitened_data.dat`.

## Suggested sequencing

1. Branch/worktree off `mojito-noise`.
2. Implement step 0 (cadence from file) + step 2 (FFT norm + model conversion) together.
3. Delete both hacks.
4. Verify per-layer ratio flat at 1.0; add assertion.
5. Cross-check the 0.25 s file too (should now also whiten).
