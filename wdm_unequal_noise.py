"""Build the WDM (wavelet-domain) PSD of the *unequal-arm* TDI noise model and plot
it in the style of wavelet_noise_plotting.ipynb.

The shipped data/full_noise_model.dat is equal-arm + stationary (XX=YY=ZZ, no time
dependence) and only covers 0.39-8 mHz. Here we keep its time-pixel grid but extend
the WDM frequency layers up to FMAX, normalize with GLASS's phase-PSD convention
(frac-freq PSD / (2pi f)^2, scaled by a calibrated constant K), and evaluate every
channel with the arm delays at each time pixel (orbit breathing). Output:
data/full_noise_model_unequal.dat (8-col: t f XX YY ZZ XY YZ ZX) and a 3-panel PNG.
"""
import numpy as np
import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import colors

CL = 299792458.0
P = 2 * np.pi
H5 = "NOISE_731d_0.25s_L1_source0_0_20251204T192628758262Z.h5"

# --- WDM grid: keep GLASS time pixels, extend frequency layers up to FMAX ------
DTPIX = 7680.0                                         # WDM pixel time = 1/(2 df)
df = 1.0 / (2 * DTPIX)
FMAX = 3.4e-2                                          # past the ~30 mHz gen-2 null
M0 = 6                                                  # first layer in GLASS files
t_pix = np.unique(np.loadtxt("data/full_noise_model.dat", usecols=(0,)))
Nt = t_pix.size
m = np.arange(M0, int(np.ceil(FMAX / df)) + 1)
f_layers = m * df; Nf = f_layers.size
w = P * f_layers                                        # (Nf,)

# --- per-pixel arm delays from the LTTs (orbit breathing) ----------------------
h = h5py.File(H5, "r")
DT = 0.25                                               # LTT sample spacing [s]
Nsamp = h["ltts/ltt_12"].shape[0]
idx = np.clip((t_pix / DT).astype(int), 0, Nsamp - 1)
dl = {ij: h["ltts/ltt_%d" % ij][:][idx] for ij in (12, 13, 21, 23, 31, 32)}
dbar = {(1, 2): (dl[12] + dl[21]) / 2,
        (1, 3): (dl[13] + dl[31]) / 2,
        (2, 3): (dl[23] + dl[32]) / 2}
Dd = {(1, 2): dl[12] - dl[21], (1, 3): dl[13] - dl[31], (2, 3): dl[23] - dl[32]}
dmean = np.mean([dl[ij] for ij in dl], axis=0)         # equal-arm delay per pixel
dref = float(dmean.mean())                             # scalar reference armlength


def levels(fq):
    fonc = P * fq / CL
    Saa = 9.00e-30 * (1 + (0.4e-3 / fq)**2) * (1 + (fq / 8e-3)**4)
    Spm = Saa * (P * fq)**-4 * fonc**2
    Sop = 2.25e-22 * (1 + (2e-3 / fq)**4) * fonc**2
    return Spm, Sop


Spm, Sop = levels(f_layers)                            # (Nf,)


def un_XX(di, dk):                                     # unequal frac-freq XX, arms i,k
    a, b = np.outer(di, w), np.outer(dk, w)
    H1 = 8 * (np.sin(a)**2 * (3 + np.cos(2 * b)) + np.sin(b)**2 * (3 + np.cos(2 * a)))
    H2 = 8 * (np.sin(a)**2 + np.sin(b)**2)
    return (H1 * Spm + H2 * Sop) * 4 * np.sin(a + b)**2


def un_CSD(dij, dik, djk, Dij):                        # unequal complex frac-freq CSD
    a, b, c = np.outer(dij, w), np.outer(dik, w), np.outer(djk, w)
    ph = np.exp(-1j * (np.outer(dik - djk, w) + np.outer(Dij, w) / 2))
    g2 = 4 * np.sin(a + b) * np.sin(a + c) * np.exp(-1j * np.outer(dik - djk, w))
    base = np.cos(a) * np.sin(b) * np.sin(c) * ph * g2
    return (-32 * base) * Spm + (-8 * base) * Sop


# --- calibrate GLASS normalization K against the shipped equal-arm file --------
cal = np.loadtxt("data/instrument_noise_model.dat")    # f XX YY ZZ XY YZ ZX
fc = cal[:, 0]
Spc, Soc = levels(fc); xc = P * fc * dref
eq_c = (16 * np.sin(xc)**2 * (3 + np.cos(2 * xc)) * Spc
        + 16 * np.sin(xc)**2 * Soc) * 4 * np.sin(2 * xc)**2
K = float(np.median(cal[:, 1] * (P * fc)**2 / eq_c))   # phase-PSD scale
norm = K / w**2                                         # frac-freq -> GLASS units, (Nf,)

# --- channels: X on S/C1 (arms 12,13), Y on S/C2 (12,23), Z on S/C3 (13,23) ----
XX = norm * un_XX(dbar[(1, 2)], dbar[(1, 3)])
YY = norm * un_XX(dbar[(1, 2)], dbar[(2, 3)])
ZZ = norm * un_XX(dbar[(1, 3)], dbar[(2, 3)])
XY = norm * un_CSD(dbar[(1, 2)], dbar[(1, 3)], dbar[(2, 3)], Dd[(1, 2)])
YZ = norm * un_CSD(dbar[(2, 3)], dbar[(1, 2)], dbar[(1, 3)], Dd[(2, 3)])
ZX = norm * un_CSD(dbar[(1, 3)], dbar[(2, 3)], dbar[(1, 2)], -Dd[(1, 3)])

# --- write same 8-col format (CSD columns: real part, like GLASS) --------------
T = np.repeat(t_pix, Nf); F = np.tile(f_layers, Nt)
out = np.column_stack([T, F, XX.ravel(), YY.ravel(), ZZ.ravel(),
                       XY.real.ravel(), YZ.real.ravel(), ZX.real.ravel()])
np.savetxt("data/full_noise_model_unequal.dat", out, fmt="%.6e")


def spec(ax, z, title, cmap, norm):
    tp = np.append(t_pix, 2 * t_pix[-1] - t_pix[-2]) / 86400.0          # days
    fp = np.append(f_layers, 2 * f_layers[-1] - f_layers[-2]) * 1e3     # mHz
    mesh = ax.pcolormesh(tp, fp, z.T, cmap=cmap, norm=norm,
                         shading="flat", rasterized=True)
    ax.set_xlabel("t [days]"); ax.set_ylabel("f [mHz]"); ax.set_title(title)
    plt.colorbar(mesh, ax=ax)


# --- the gen-2 null wanders in frequency as the arms breathe -------------------
# |chan|^2 carries an overall 4 sin^2(omega(d_ij+d_ik)) factor -> exact null at
# f_null = 1/(2(d_ij+d_ik)), a different (moving) frequency for each channel.
fnull = {"X": 1 / (2 * (dbar[(1, 2)] + dbar[(1, 3)])),
         "Y": 1 / (2 * (dbar[(1, 2)] + dbar[(2, 3)])),
         "Z": 1 / (2 * (dbar[(1, 3)] + dbar[(2, 3)]))}
allnull = np.concatenate(list(fnull.values()))
lo_f, hi_f = allnull.min() - 4 * df, allnull.max() + 4 * df   # zoom band [Hz]

tp = np.append(t_pix, 2 * t_pix[-1] - t_pix[-2]) / 86400.0
fp = np.append(f_layers, 2 * f_layers[-1] - f_layers[-2]) * 1e3
td = t_pix / 86400.0
band = np.ma.masked_outside(XX, *np.nanpercentile(
    XX[:, (f_layers > lo_f) & (f_layers < hi_f)], [1, 99.9]))
norm0 = colors.LogNorm(vmin=band.compressed().min(), vmax=band.compressed().max())

fig, ax = plt.subplots(1, 3, figsize=(16, 4.6), sharey=True)
for a, (ch, z) in zip(ax, [("X", XX), ("Y", YY), ("Z", ZZ)]):
    m = a.pcolormesh(tp, fp, z.T, cmap="viridis", norm=norm0,
                     shading="flat", rasterized=True)
    a.plot(td, fnull[ch] * 1e3, "w--", lw=1.0,
           label=r"$f_{\rm null}=1/[2(\bar d_{ij}+\bar d_{ik})]$")
    a.set_ylim(lo_f * 1e3, hi_f * 1e3)
    a.set_xlabel("t [days]"); a.set_title(r"$|%s|^2$  gen-2 null" % ch)
    a.legend(loc="upper right", fontsize=8, framealpha=0.7)
ax[0].set_ylabel("f [mHz]")
fig.colorbar(m, ax=ax, label=r"WDM PSD", fraction=0.025, pad=0.01)
fig.savefig("wdm_unequal_noise.png", dpi=130, bbox_inches="tight")
print("null band: %.3f-%.3f mHz   X-null swing = %.4f mHz over %.0f days"
      % (lo_f * 1e3, hi_f * 1e3,
         (fnull["X"].max() - fnull["X"].min()) * 1e3, t_pix[-1] / 86400))
print("wrote data/full_noise_model_unequal.dat and wdm_unequal_noise.png")
print("Nt=%d Nf=%d  f=[%.3f,%.3f] mHz  span=%.1f days  K=%.4e"
      % (Nt, Nf, f_layers[0] * 1e3, f_layers[-1] * 1e3, t_pix[-1] / 86400, K))
print("gen-2 TDI null at f=1/(4 d̄)=%.3f mHz" % (1 / (4 * dref) * 1e3))
print("channel asym (Z-Y)/X: max = %.3e   Im(XY)/XX: max = %.3e"
      % (np.max(np.abs((ZZ - YY) / XX)), np.max(np.abs(XY.imag / XX))))
