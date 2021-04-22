# THIS PROGRAM FINDS THE FIR COEFFICIENTS FOR A BANDPASS FILTER USING THE
# KAISER WINDOW AND THE DESIGN SPECIFICATIONS

import numpy as np
import matplotlib.pyplot as plt

L = 114 # Filter number
Fs = 48 # Sampling freqency (kHz)
const = 2 * np.pi /  Fs # Constant used to get the normalized digital freqencies
delta = 0.15 # The permissible filter amplitude deviation from unity

# Bandpass filter specifications
 
# Passband: F_p2 < F_p1
F_p1 = 3 + 0.6 * (L - 107)
F_p2 = 3 + 0.6 * (L - 109)

delF = 0.3 # Transition band

# Stopband F_s2 < F_p21; F_p1 < F_s1
F_s1 = F_p1 + 0.3
F_s2 = F_p2 - 0.3

# Normalized digital filter specifications (radians/sec)
omega_p1 = const * F_p1
omega_p2 = const * F_p2

omega_c = (omega_p1 + omega_p2)/2;
omega_l = (omega_p1 - omega_p2)/2;

omega_s1 = const * F_s1;
omega_s2 = const * F_s2;
delomega = 2 * np.pi * delF / Fs


# The Kaiser window design
A = -20 * np.log10(delta)
N = np.ceil((A - 8)/(4.57 * delomega))

N = 100
n = np.arange(-N, N)
hlp = np.sin(n * omega_l)/(n * np.pi)
hlp[N] = omega_l / np.pi # This handles the case when n[i] = 0

# The Bandpass filter
hbp = 2 * hlp * np.cos(n * omega_c)

# The bandpass filter plot
omega = np.arange(-np.pi/2, np.pi/2, np.pi/200)

Hbp = abs(np.polyval(hbp, np.exp(-1j * omega)))

plt.figure()
plt.plot(omega / np.pi, Hbp)
plt.xlabel('$\omega/\pi$')
plt.ylabel('$|H_{bp}(\omega)|$')
plt.grid()
plt.show()
plt.savefig("../figs/fir_kaiser_window.eps")
fir_coeff = hbp