# IIR  FILTER DESIGN USING THE CHEBYSCHEV APPROXIMATION
# The order of the filter is N = 4

from funcs import *
import numpy as np
import matplotlib.pyplot as plt
#If using termux
import subprocess
import shlex
#end if

N = 4
L = 114 # Filter number
Fs = 48 # Sampling freqency (kHz)
const = 2 * np.pi /  Fs # Constant used to get the normalized digital freqencies
 
# Passband: F_p2 < F_p1
F_p1 = 3 + 0.6 * (L - 107)
F_p2 = 3 + 0.6 * (L - 109)
omega_p1 = const * F_p1
omega_p2 = const * F_p2
Omega_p1 = np.tan(omega_p1/2)
Omega_p2 = np.tan(omega_p2/2)
Omega_0 = np.sqrt(Omega_p1 * Omega_p2)
B = Omega_p1 - Omega_p2

epsilon = 0.4

# The analog lowpass filter
p, G_lp = lp_stable_cheb(epsilon, N)

Omega_L = np.arange(-2, 2 + 0.01, 0.01)
H_analog_lp = G_lp * abs(1/np.polyval(p, 1j * Omega_L))

plt.figure()
plt.plot(Omega_L, H_analog_lp)
plt.xlabel('$\Omega$')
plt.ylabel('$|H_{a,LP}(j\Omega)|$')
plt.title("Analog Lowpass filter")
plt.grid()

# The analog bandpass filter
num, den, G_bp = lpbp(p,Omega_0,B,Omega_p1)

Omega = np.arange(-0.65, 0.65 + 0.01, 0.01)
H_analog_bp = G_bp * abs(np.polyval(num,1j * Omega) / np.polyval(den,1j * Omega))

plt.figure()
plt.plot(Omega, H_analog_bp)
plt.xlabel('$\Omega$')
plt.ylabel('$|H_{a,BP}(j\Omega)|$')
plt.title("Analog Bandpass filter")
plt.grid()
# plt.savefig('../../figs/iir/analog_bp.eps')
# plt.savefig('../../figs/iir/analog_bp.pdf')
subprocess.run(shlex.split("termux-open ../../figs/iir/analog_bp.pdf"))

# The digital bandpass filter
dignum, digden, G = bilin(den, omega_p1)

omega = np.arange(-2 * np.pi / 5, 2 * np.pi / 5 + (np.pi / 1000), np.pi / 1000)
H_dig_bp = G * abs(np.polyval(dignum, np.exp(-1j*omega))/np.polyval(digden, np.exp(-1j*omega)))

plt.figure()
plt.plot(omega/np.pi, H_dig_bp)
plt.xlabel('$\omega/\pi$')
plt.ylabel('$|H_{d,BP}(\omega)|$')
plt.title("Digital Bandpass filter")
plt.grid()
# plt.savefig('../../figs/iir/digital_bp.eps')
# plt.savefig('../../figs/iir/digital_bp.pdf')
subprocess.run(shlex.split("termux-open ../../figs/iir/digital_bp.pdf"))

# plt.show()

iir_num = G * dignum
iir_den = digden

np.savetxt('./iir_num.dat', iir_num)
np.savetxt('./iir_den.dat', iir_den)
np.savetxt('./dignum.dat', dignum)
np.savetxt('./digden.dat', digden)
np.savetxt('./G.dat', np.array([G]))