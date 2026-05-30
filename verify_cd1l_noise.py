"""Verify GLASS analytic noise model against CD1L (MOJITO) ground-truth noise_estimates.

Ground truth: noise_estimates/XYZ (LISAModels covariance), total-freq units (Hz^2/Hz).
GLASS works on time series divided by laser_frequency nu0 -> fractional-freq units,
so we compare GLASS model to (ground truth / nu0^2).

Covers both the PSD (XX) and the CSD (XY), in two arm treatments:
  - equal-arm:  single averaged 2-way delay d̄ (CSD is purely real).
  - unequal-arm: per-arm 2-way delays d̄_ij + delay difference Δd_ij; this is the
    only treatment that produces the (small) imaginary part of the CSD.
TDI generation 2 throughout. Transfer functions follow NOISE_SPEC.md (Dam_2023).
"""
import sys
import h5py
import numpy as np

FN = sys.argv[1] if len(sys.argv) > 1 else \
    "NOISE_731d_0.25s_L1_source0_0_20251204T192628758262Z.h5"
CLIGHT = 299792458.0
PI2 = 2 * np.pi

f = h5py.File(FN, "r")
nu0 = float(f.attrs["laser_frequency"])

# Per-link one-way LTTs (epoch 0). Build the three two-way arm averages
# d̄_ij = (d_ij + d_ji)/2 and the 2-way delay difference Δd_ij = d_ij - d_ji.
dl = {ij: f["ltts/ltt_%d" % ij][:2000].mean()
      for ij in (12, 13, 21, 23, 31, 32)}
dbar = {(1, 2): (dl[12] + dl[21]) / 2,
        (1, 3): (dl[13] + dl[31]) / 2,
        (2, 3): (dl[23] + dl[32]) / 2}
Dd12 = dl[12] - dl[21]                       # Δd for the X-Y pair (i=1, j=2)
dmean = np.mean([dl[ij] for ij in dl])       # single equal-arm delay
fstar = 1.0 / (PI2 * dmean)

freq = np.logspace(np.log10(1e-5), np.log10(1.0), 1000)  # match file sampling
gt_XX = np.real(f["noise_estimates/XYZ"][0, :, 0, 0]) / nu0**2
gt_XY = f["noise_estimates/XYZ"][0, :, 0, 1] / nu0**2


def levels(model, fq):
    """Single-link input PSDs (fractional-freq units): test-mass Spm, OMS Sop."""
    fonc = PI2 * fq / CLIGHT
    Sa_a = 9.00e-30 * (1 + (0.4e-3 / fq)**2) * (1 + (fq / 8e-3)**4)
    if model == "mojito":                    # + lowfreq-relax frelax=1e-4
        Sa_a = Sa_a * (1 + (1e-4 / fq)**4)
    Spm = Sa_a * (PI2 * fq)**-4 * fonc**2
    Sop = 2.25e-22 * (1 + (2e-3 / fq)**4) * fonc**2
    return Spm, Sop


def glass_eq(model, fq, kind):
    """Equal-arm gen-2 TDI spectrum (XX real, XY real)."""
    Spm, Sop = levels(model, fq)
    x = PI2 * fq * dmean                      # = ω d̄
    g2 = 4 * np.sin(2 * x)**2                  # generation-2 factor
    if kind == "XX":
        H1 = 16 * np.sin(x)**2 * (3 + np.cos(2 * x))   # 4 C_XX (3 + cos 2ωd̄)
        H2 = 16 * np.sin(x)**2                          # 4 C_XX
    else:                                     # XY
        Cxy = -4 * np.sin(x) * np.sin(2 * x)
        H1, H2 = 4 * Cxy, Cxy
    return (H1 * Spm + H2 * Sop) * g2


def glass_uneq(model, fq, kind):
    """Unequal-arm gen-2 TDI spectrum. XX for channel X (arms 12,13);
    XY for the X-Y pair (i=1, j=2, k=3). XY is complex."""
    Spm, Sop = levels(model, fq)
    w = PI2 * fq
    if kind == "XX":
        a, b = w * dbar[(1, 2)], w * dbar[(1, 3)]
        H1 = 8 * (np.sin(a)**2 * (3 + np.cos(2 * b))
                  + np.sin(b)**2 * (3 + np.cos(2 * a)))
        H2 = 8 * (np.sin(a)**2 + np.sin(b)**2)
        g2 = 4 * np.sin(a + b)**2
        return (H1 * Spm + H2 * Sop) * g2
    # XY: d̄_ij=arm12, d̄_ik=arm13, d̄_jk=arm23, Δd_ij=Δd12
    dij, dik, djk = dbar[(1, 2)], dbar[(1, 3)], dbar[(2, 3)]
    a, b, c = w * dij, w * dik, w * djk
    phase = np.exp(-1j * (w * (dik - djk) + w * Dd12 / 2))
    g2 = 4 * np.sin(a + b) * np.sin(a + c) * np.exp(-1j * w * (dik - djk))
    base = np.cos(a) * np.sin(b) * np.sin(c) * phase * g2
    return (-32 * base) * Spm + (-8 * base) * Sop


lines = ["nu0=%.4e  mean_ltt=%.5f s  fstar=%.5f Hz" % (nu0, dmean, fstar),
         "arms (d̄ s): 12=%.5f 13=%.5f 23=%.5f   Δd12=%.3e s"
         % (dbar[(1, 2)], dbar[(1, 3)], dbar[(2, 3)], Dd12),
         "ratios = GLASS_model / noise_estimates (gen-2):",
         "         XX(PSD)        XY Re(CSD)     XY Im(CSD)",
         " f[mHz]  eq    uneq   eq     uneq   gt[abs]   uneq/gt"]
for model in ("scirdv1", "mojito"):
    lines.append("-- model: %s" % model)
    for ftest in (0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50):
        i = int(np.argmin(np.abs(freq - ftest * 1e-3)))
        xx_eq = glass_eq(model, freq, "XX")[i] / gt_XX[i]
        xx_un = glass_uneq(model, freq, "XX")[i] / gt_XX[i]
        xy_eq = glass_eq(model, freq, "XY")[i].real / gt_XY[i].real
        xy_un = glass_uneq(model, freq, "XY")
        re_un = xy_un[i].real / gt_XY[i].real
        im_un = (xy_un[i].imag / gt_XY[i].imag) if gt_XY[i].imag else 0.0
        lines.append(" %6.2f  %5.3f %5.3f  %5.3f %5.3f  %9.2e %6.3f"
                     % (freq[i] * 1e3, xx_eq, xx_un, xy_eq, re_un,
                        gt_XY[i].imag, im_un))

with open("/tmp/v.txt", "w") as fh:
    fh.write("\n".join(lines) + "\n")
for l in lines:
    print(l)
