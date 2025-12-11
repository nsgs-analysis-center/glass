import numpy as np
import scipy.signal
import scipy.fft
import matplotlib.pyplot as plt

N_gen = 8192*4

fs = 1.0
dt = 1/fs

# Generate white noise and filter it to give it an interesting PSD
#w = np.random.normal(0, 1, (N_gen,))
#x = scipy.signal.filtfilt(np.array([1, -2, -3, -4, 5]), np.array([1]), w) + 2

# Generate data with a theoretically expected PSD
farr = np.fft.rfftfreq(N_gen,dt)
fref = 1e-1 
#expected_PSD = 1e-10*np.ones(farr.shape[0]-1)
expected_PSD = 1e-10*(farr[farr>0] / fref)**(-2/3)
expected_PSD = np.hstack((np.array([0.0]),expected_PSD)) # fix DC bin
Nf = farr.shape[0]

periodograms = []
for _ in range(100):
    X = np.sqrt(expected_PSD/2)*np.random.randn(Nf) + 1j*np.sqrt(expected_PSD/2)*np.random.randn(Nf)
    x = np.real(np.fft.irfft(X)) * np.sqrt(N_gen/(2*dt))
    f, P = scipy.signal.periodogram(x, fs, 'boxcar', nfft=N_gen, detrend=False)
    periodograms.append(P)

# Find PSD using built-in function and manually:
# Built-in function
f, P = scipy.signal.periodogram(x, fs, 'boxcar', nfft=N_gen, detrend=False)
plt.loglog(f,P,label="periodogram")
plt.loglog(farr[farr>0],expected_PSD[farr>0],label="expected PSD")
plt.loglog(f,np.mean(periodograms,axis=0),label="average of many periodograms")
plt.ylim([1e-14,1e-7])
plt.legend()
plt.show()

# whitening test
bins = np.linspace(-6,6,100)
wht = np.real(X[1:]) / np.sqrt(expected_PSD[1:]/2)
counts,bins,_ = plt.hist(wht, bins=bins, density=True)
plt.plot(bins, np.max(counts)*np.exp(-bins**2/2), ls='dashed', color='black')
plt.yscale('log')
plt.show()
print(f"stdev: {np.std(wht)}")

# Manually
N = len(x)
X = np.fft.rfft(x) * np.sqrt(2*dt/N)
P_fft_one_sided = np.abs(X*np.conj(X))
P_fft_one_sided[0] = P_fft_one_sided[0] / 2
P_fft_one_sided[-1] = P_fft_one_sided[-1] / 2
print(np.allclose(P, P_fft_one_sided))  # True, P matches P_fft_one_sided
'''
N = len(x)
X = scipy.fft.fft(x) / np.sqrt(N) * np.sqrt(2*dt)
P_fft = np.abs(X*np.conj(X))
P_fft_one_sided = P_fft[0:int(N/2)+1]  # P_fft_one_sided is identical to P
P_fft_one_sided[0] = P_fft_one_sided[0] / 2
P_fft_one_sided[-1] = P_fft_one_sided[-1] / 2
print(np.allclose(P, P_fft_one_sided))  # True, P matches P_fft_one_sided
'''

# Now undo the manual operation of P_fft_one_sided to get back to a time series, insert your PSD here

'''
N_P = len(P_fft_one_sided)  # Length of PSD
N = 2*(N_P - 1)

# Because P includes both DC and Nyquist (N/2+1), P_fft must have 2*(N_P-1) elements
P_fft_one_sided[0] = P_fft_one_sided[0] * 2
P_fft_one_sided[-1] = P_fft_one_sided[-1] * 2
P_fft_new = np.zeros((N,), dtype=complex)
P_fft_new[0:int(N/2)+1] = P_fft_one_sided
P_fft_new[int(N/2)+1:] = P_fft_one_sided[-2:0:-1]
'''
P_fft_one_sided[0] *= 2
P_fft_one_sided[-1] *= 2
X_new = np.sqrt(P_fft_one_sided)

# Create random phases for all FFT terms other than DC and Nyquist
phases = np.hstack((np.array([0]),np.random.uniform(0, 2*np.pi, (int(N/2)-1,)), np.array([0])))
X_new = X_new * np.exp(2j*phases)

# This is the new time series with a given PSD
X_new = X_new * np.sqrt(N) / np.sqrt(2*dt)
x_new = np.real(np.fft.irfft(X_new))

# Verify that P matches P_new
f_new, P_new = scipy.signal.periodogram(x_new, fs, 'boxcar', nfft=N_gen, detrend=False)

import matplotlib.pyplot as plt
plt.loglog(f_new, P_new)
plt.loglog(f, P)
plt.loglog(farr, expected_PSD)
plt.ylim([1e-14,1e-7])
plt.show()
print(np.allclose(P, P_new))  # True, P matches P_new
