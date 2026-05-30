# LISA Models — Analytic TDI XYZ Noise Models

A summary of the analytical residual-noise models implemented in `lisamodels`, with
emphasis on the **TDI X/Y/Z** spectra (and the A/E/T and η spectra derived from them).

Source files:
- `lisamodels/noise_models.py` — noise sources, input PSDs, TDI assembly (`NoiseSet`).
- `lisamodels/noise_transfer_functions.py` — the `H1…H9` TDI transfer functions.

Primary reference: **Quang Nam et al., *Time-delay interferometry noise transfer
functions for LISA*, Phys. Rev. D 108, 082004 (2023)** (`Dam_2023`). The AET
expressions in the code cite Eqs. 58–61 of that paper. Conventions follow
`Bayle_2023` (PRD 107, 083019) and the LISA TN *Analytical formulation for
performance model* (LISA-LCST-INST-TN).

> **Scope warning (from the code):** the package currently models **only the
> MOJITO-light secondary noise set** — test-mass acceleration noise and OMS noise
> in the science (SCI) interferometer. See "MOJITO / data generation" below.

---

## 1. Conventions and units

- Angular frequency: `ω = 2πf`.
- `c` is the speed of light, `ν₀ = central_freq` the laser central frequency
  (the simulations use `ν₀ = 2.816e14 Hz`).
- A link `ij` carries light from spacecraft `j` to spacecraft `i`; `d_ij` is its
  one-way light-travel time (LTT) in seconds.
- All input PSDs and all TDI PSDs/CSDs are returned in **total frequency units**,
  `Hz²/Hz`.
- Spectra are evaluated on a grid of (epoch `t`, frequency `f`); returned arrays are
  shape `(#epochs, #frequencies)`. The epoch dependence enters only through the
  delays (orbits), which are spline-interpolated by the `Delays` class.

Two LTT combinations drive the unequal-arm transfer functions:
- average 2-way delay: `d̄_ij = (d_ij + d_ji) / 2`
- 2-way delay difference: `Δd_ij = d_ij − d_ji`

In the equal-arm approximation a single `d̄` = average of all six links is used.

---

## 2. Input (single-link) noise PSDs

Each noise source produces a "link-level" input PSD `S(f)` that is then propagated
through a TDI transfer function. Both implemented sources are assumed **independent
across the six MOSAs** and independent of each other. Per-MOSA values may be a single
float (same for all six) or a dict keyed by MOSA (`12, 23, 31, 13, 32, 21`).

### 2.1 Test-mass acceleration noise — `TestMassAccelNoise` (→ uses `H1`)

Acceleration ASD `A` in `m·s⁻²/√Hz`. Input PSD in acceleration:

```
S_acc(f) = A² · [1 + (f_knee / f)²] · [1 + (f / f_break)⁴]
```

Optional `shape = "lowfreq-relax"` adds a pessimistic low-frequency `f⁻⁴` tail
(extends below the nominal 1e-4 Hz band):

```
S_acc(f) = … × [1 + (f_relax / f)⁴]      # when shape == "lowfreq-relax"
```

Conversion acceleration → total frequency `Hz²/Hz`:

```
S(f) = S_acc(f) · ( ν₀ / (2π f c) )²
```

Parameters (keywords pulled from the instrument config):
`central_freq`, `testmass_asds`, `testmass_fknees`, `testmass_fbreak`,
`testmass_frelax`, `testmass_shape` (`"original"` | `"lowfreq-relax"`).

### 2.2 OMS science-channel noise — `OMSNoiseSci` (→ uses `H2`)

Displacement ASD `A` in `m/√Hz`. Input PSD in displacement:

```
S_disp(f) = A² · [1 + (f_knee / f)⁴]
```

Conversion displacement → total frequency `Hz²/Hz`:

```
S(f) = S_disp(f) · ( 2π f ν₀ / c )²
```

Parameters: `central_freq`, `oms_sci_carrier_asds`, `oms_fknees`.

---

## 3. From input PSD to TDI X/Y/Z spectra

The TDI spectrum for one source is `(transfer function) × (input PSD)`. The transfer
function depends on the TDI generation (`1` or `2`) and on whether arms are treated as
equal or unequal.

| Source              | Equal-arm TF        | Unequal-arm PSD TF          | Unequal-arm CSD TF          |
|---------------------|---------------------|-----------------------------|-----------------------------|
| Test-mass accel.    | `h1`                | `h1_xx_unequal_armlength`   | `h1_xy_unequal_armlength`   |
| OMS (SCI)           | `h2`                | `h2_xx_unequal_armlength`   | `h2_xy_unequal_armlength`   |

### 3.1 Equal-arm case (`equal_armlength=True`)

A single averaged delay `d̄` is used. PSDs of `XX/YY/ZZ` and CSD of `XY/YZ/ZX` are

```
S_TDI = H(d̄, f, tdi, gen) · S_input(f)
```

with `H = h1` (test mass) or `h2` (OMS). The dispatcher `h1`/`h2` selects the right
closed form for the requested combination (`XX`, `AA`, `TT`, or `XY`). For X/Y/Z
CSDs all three off-diagonal pairs share the same `XY` transfer function.

### 3.2 Unequal-arm case (default, `equal_armlength=False`)

For a Michelson channel starting on S/C `i` (X→1, Y→2, Z→3), with `j`, `k` the other
two spacecraft, the PSD depends on the two 2-way delays `d̄_ij`, `d̄_ik`:

```
S_XX = h?_xx_unequal_armlength(f, d̄_ij, d̄_ik, gen) · S_input(f)
```

The CSD between channels `i` and `j` depends on `d̄_ij, d̄_ik, d̄_jk` and the delay
difference `Δd_ij`:

```
S_ij = h?_xy_unequal_armlength(f, d̄_ij, d̄_ik, d̄_jk, Δd_ij, gen) · S_input(f)
```

(e.g. `ZX` ↔ `i=3, j=1, k=2`). These CSDs are complex.

### 3.3 A/E/T spectra (derived from X/Y/Z)

A, E, T spectra are **not** computed from their own transfer functions; they are
assembled from the X/Y/Z PSDs and CSDs using `Dam_2023` Eqs. 58–61
(`GenericNoiseSource.compute_psd_tdi` / `compute_csd_tdi`):

```
S_AA = ( S_ZZ + S_XX − 2 Re S_ZX ) / 2
S_EE = ( S_XX + 4 S_YY + S_ZZ − 2 Re( 2 S_XY − S_ZX + 2 S_YZ ) ) / 6
S_TT = ( S_XX + S_YY + S_ZZ + 2 Re( S_XY + S_ZX + S_YZ ) ) / 3

Cov(A,E) = ( −S_XX + S_ZZ + 2 S_XY − 2 S_ZY − 2i Im S_XZ ) / √12
Cov(E,T) = ( S_XX − 2 S_YY + S_ZZ + 2 Re S_ZX + S_XY − 2 S_YX + S_ZY − 2 S_YZ ) / √18
Cov(T,A) = ( −S_XX + S_ZZ + S_YZ − S_YX − 2i Im S_ZX ) / √6
```

The 3×3 AET covariance matrix is also available directly via the orthogonal
`get_AET_covariance_matrix`, which rotates the XYZ covariance with the standard
`XYZ→AET` matrix:

```
        | −1/√2     0     1/√2 |
M_AET = |  1/√6  −2/√6   1/√6 |
        |  1/√3   1/√3   1/√3 |
```

### 3.4 η (single-link) spectra

The `compute_psd_csd_eta` / `compose_csd_eta` / `get_eta_covariance_matrix` API gives
the 6×6 covariance of the intermediate η variables:

- **Test mass:** non-zero only for same link and reverse link:
  ```
  CSD[η_ij, η_ij] = S^acc_ij + S^acc_ji
  CSD[η_ij, η_ji] = S^acc_ij · e^{+iωL_ji} + S^acc_ji · e^{−iωL_ij}
  ```
- **OMS (SCI):** diagonal only — `CSD[η_ij, η_kl] = δ_ik δ_jl · S^disp_ij`.

---

## 4. Transfer-function catalog (`noise_transfer_functions.py`)

### 4.1 Common factors

Equal-arm forms are built from two reusable factors (gen-1 base, gen-2 multiplier):

```
C_XX = 4 sin²(ω d̄)                       # ×4 sin²(2ω d̄) for generation 2
C_XY = −4 sin(ω d̄) sin(2ω d̄)            # ×4 sin²(2ω d̄) for generation 2
```

### 4.2 H1 and H2 — the wired-in models

These two families are the ones actually used by `NoiseSet` (test mass → H1, OMS → H2).

**H1 (equal arm):**
```
H1_XX = 4 C_XX · [3 + cos(2ω d̄)]
H1_AA = 4 C_XX · [3 + 2cos(ω d̄) + cos(2ω d̄)]      (also used for EE)
H1_TT = 32 C_XX · sin⁴(ω d̄ / 2)
H1_XY = 4 C_XY
```

**H1 (unequal arm):**
```
H1_XX = 8 [ sin²(ω d̄_ij)(3 + cos 2ω d̄_ik) + sin²(ω d̄_ik)(3 + cos 2ω d̄_ij) ]
        × 4 sin²(ω(d̄_ij + d̄_ik))                       # gen-2 factor

H1_ij = −32 cos(ω d̄_ij) sin(ω d̄_ik) sin(ω d̄_jk)
        · e^{−iω(d̄_ik − d̄_jk + Δd_ij/2)}
        × 4 sin(ω(d̄_ij+d̄_ik)) sin(ω(d̄_ij+d̄_jk)) e^{−iω(d̄_ik−d̄_jk)}   # gen-2 factor
```

**H2 (equal arm):**
```
H2_XX = 4 C_XX
H2_AA = 2 C_XX · (2 + cos ω d̄)                       (also used for EE)
H2_TT = 4 C_XX · (1 − cos ω d̄)
H2_XY = C_XY
```

**H2 (unequal arm):**
```
H2_XX = 8 [ sin²(ω d̄_ij) + sin²(ω d̄_ik) ]
        × 4 sin²(ω(d̄_ij + d̄_ik))                       # gen-2 factor

H2_ij = −8 cos(ω d̄_ij) sin(ω d̄_ik) sin(ω d̄_jk)
        · e^{−iω(d̄_ik − d̄_jk + Δd_ij/2)}
        × 4 sin(ω(d̄_ij+d̄_ik)) sin(ω(d̄_ij+d̄_jk)) e^{−iω(d̄_ik−d̄_jk)}   # gen-2 factor
```

### 4.3 H3–H9 — implemented but not yet wired into `NoiseSet`

These closed forms exist in the module and are documented/tested, but no noise
source currently maps to them (they cover correlated/anti-correlated noises not in
the MOJITO-light set). They are equal-arm only and provide `XX`, `AA`, `TT`, `XY`.
Each is a multiple of `C_XX`/`C_XY`:

| TF  | Representative noise (per docs)                  | `XX`                         | `AA`                            | `TT`                         |
|-----|--------------------------------------------------|------------------------------|---------------------------------|------------------------------|
| H3  | OMS noise in TMI channel                         | `C_XX (3 + cos 2ω d̄)`        | `C_XX (3 + 2cos + cos 2ω d̄)`    | `8 C_XX sin⁴(ω d̄/2)`         |
| H4  | Anti-correlated TMI / adjacent-MOSA jitter       | `8 C_XX (2 + cos 2ω d̄)`      | `4 C_XX (1 + 2cos ω d̄)²`        | `64 C_XX sin⁴(ω d̄/2)`        |
| H5  | Fully-correlated adjacent TM accel.              | `8 C_XX`                     | `12 C_XX`                       | `0`                          |
| H6  | Anti-correlated OMS in adjacent SCI/REF          | `6 C_XX`                     | `C_XX (5 + 4cos ω d̄)`           | `−8 C_XX (cos ω d̄ − 1)`      |
| H7  | Correlated OMS in adjacent SCI/REF/TMI           | `2 C_XX`                     | `3 C_XX`                        | `0`                          |
| H8  | Anti-correlated adjacent TMI                     | `2 C_XX (2 + cos 2ω d̄)`      | `C_XX (1 + 2cos ω d̄)²`          | `16 C_XX sin⁴(ω d̄/2)`        |
| H9  | Uncorrelated TMI backlink                        | —                            | `C_XX (3 + 2cos + cos 2ω d̄)`    | `8 C_XX sin⁴(ω d̄/2)`         |

Off-diagonal (`XY`) forms for H3–H8 (H9 has none):
`H3_XY = C_XY`, `H4_XY = 4 C_XX(1 − 4cos ω d̄)`, `H5_XY = −4 C_XX`,
`H6_XY = C_XX(1 − 4cos ω d̄)`, `H7_XY = −C_XX`, `H8_XY = C_XX(1 − 4cos ω d̄)`.

---

## 5. Public API (`NoiseSet`)

Construct from a `lisainstrument` `Instrument`, an HDF5 simulation file (written by
`Instrument.export_hdf5`, lisainstrument ≥ 2.0.0), or a config dict. Delays come from
a `Delays` object, a LISA Orbits HDF5 file, or `"fromInstrument"` (default).

```python
import numpy as np
import lisamodels.noise_models as nm

ns = nm.NoiseSet(instrument="./data/simulation.h5",
                 delays="fromInstrument",
                 noises=["testmass", "oms"])

freq = 10.0**(-5 + np.arange(6))
psd_TT = ns.compose_psd(t=ns.delays.t0, freq=freq, tdi="TT")     # PSD
csd_ZX = ns.compose_csd(t=ns.delays.t0, freq=freq, tdi="ZX")     # CSD
cov    = ns.get_AET_covariance_matrix(t=ns.delays.t0, freq=freq) # 3×3 cov
```

- PSD combinations: `XX, YY, ZZ, AA, EE, TT`. CSD pairs: `XY, YZ, ZX, AE, ET, TA`.
- Keyword args: `equal_armlength` (default `False`), `tdi_generation` (default `2`),
  `noises_list` (subset of active sources).
- Covariance helpers: `get_XYZ_covariance_matrix`, `get_AET_covariance_matrix`,
  `get_eta_covariance_matrix` → arrays of shape `(#t, #f, n, n)`.

---

## 6. MOJITO / data generation references

The repository ties these analytic models to **MOJITO** simulations produced with
`lisainstrument` + `pytdi`. (No reference to "CD1L" / a "CD-1L" dataset appears
anywhere in this repo — the only data-generation lineage present is MOJITO.)

- **"MOJITO simple" / "MOJITO-light"** = the minimal secondary-noise set: test-mass
  acceleration + OMS-SCI only, with all other instrument noises disabled. `NoiseSet`'s
  docstring explicitly states the package only supports this set.
- **`examples/Test_01_Each_All_Noises/1_datgen_sim_tdi.ipynb`** — generates the
  reference simulations. Six cases via `case_id`:

  |              | TMA | OMS | ALL (TM+OMS) |
  |--------------|-----|-----|--------------|
  | equal arms   |  1  |  2  |  3           |
  | unequal arms |  4  |  5  |  6           |

  **Case 6** (unequal arms, TM+OMS) is annotated *"i.e. MOJITO simple"*. The pipeline:
  simulate with `Instrument`, build TDI via `pytdi`
  (`compute_etas` + `compute_factorized_michelson` for X/Y/Z, `rot=0,1,2`), estimate
  PSD/CSD with Welch/CSD, and overlay the `lisamodels` analytic model.
- **`examples/Test_02_Model_Over_A_Year/study_mojito_simple_yearlong.ipynb`** —
  evaluates the "MOJITO simple" model over a full year (epoch dependence through the
  trailing-constellation ESA orbits, `orbits_unequalarms.h5`); builds
  `NoiseSet(instrument=instru, noises=["testmass", "oms"])` and plots X/Y/Z and A/E/T.
- **`examples/Test_03_Tdi_Nperseg/`** — same datagen with Welch `nperseg`/kwargs study.

Instrument noise parameters used by the example simulations (and thus the values the
analytic model is validated against):

| Parameter            | Value                                  |
|----------------------|----------------------------------------|
| `central_freq`       | `2.816e14` Hz                          |
| `testmass_asds`      | `3.0e-15` m·s⁻²/√Hz                     |
| `testmass_fknees`    | `0.4e-3` Hz                            |
| `testmass_fbreak`    | `8e-3` Hz                              |
| `testmass_shape`     | `"original"` (`testmass_frelax = 0`)   |
| `oms_asds`           | `15.0e-12` m/√Hz (link `12` only)      |
| `oms_fknees`         | `2e-3` Hz                              |
| `noises_f_min_hz`    | `5.0e-5` Hz                            |
| sampling             | `dt = 0.25 s`, `physics_upsampling = 4`|

The OMS configuration in the datagen notebook also documents the relation between the
displacement PSD and the lisainstrument filter transfer function (a `filter_approx`
high/low-frequency split), used as a cross-check at the edges of the FFT band.

---

## 7. Quick reference — which TF for which spectrum

| Want                                  | Call                                                   |
|---------------------------------------|--------------------------------------------------------|
| X/Y/Z PSD, equal arms                 | `h1`/`h2`(`d̄`, f, `"XX"`, gen)                          |
| X/Y/Z PSD, unequal arms               | `h1/h2_xx_unequal_armlength`(f, d̄_ij, d̄_ik, gen)       |
| X/Y/Z CSD, equal arms                 | `h1`/`h2`(`d̄`, f, `"XY"`, gen)                          |
| X/Y/Z CSD, unequal arms               | `h1/h2_xy_unequal_armlength`(f, d̄_ij, d̄_ik, d̄_jk, Δd_ij, gen) |
| A/E/T PSD & CSD                       | derived from X/Y/Z (Dam_2023 Eq. 58–61)                |
| η covariance                          | `compute_psd_csd_eta` per source                       |

All TDI PSDs/CSDs = transfer function × input PSD (Section 2), in `Hz²/Hz`.
