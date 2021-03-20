import soundfile as sf
import numpy as np
from scipy import signal
import matplotlib.pyplot as plt

#If using termux
import subprocess
import shlex

def to_complex(field):
	field = str(field)[2:]
	field = field[0 : len(field) - 1]
	return complex(field.replace('+-', '-').replace('i', 'j'))

#Reading the soundfile 
input_signal, fs = sf.read('../data/Sound_Noise.wav')
sampl_freq = fs
order = 4
cutoff_freq = 4000
Wn = 2 * cutoff_freq / sampl_freq
n = int(len(input_signal))
n = int(2 ** np.floor(np.log2(n)))
input_signal = input_signal[0 : n]

#Passing butterworth filter
b, a = signal.butter(order, Wn, 'low')

# Computing H(z)
w = 2 * np.pi * np.arange(n)/n

Y1 = np.loadtxt('../data/Y.dat', converters={0: to_complex}, dtype = np.complex128, delimiter = '\n') # Load the DFT computed using fft.c

output_signal_1 = np.loadtxt('../data/y.dat', converters={0: to_complex}, dtype = np.complex128, delimiter = '\n').real
output_signal_2 = signal.filtfilt(b, a, input_signal)
Y2 = np.fft.fft(output_signal_2)
Y2 = np.fft.fftshift(Y2)

sf.write('../data/Sound_With_ReducedNoise_1.wav', output_signal_1, sampl_freq)
sf.write('../data/Sound_With_ReducedNoise_2.wav', output_signal_2, sampl_freq)

#Plotting
w[0] = -np.pi
for i in range(1, n):
	w[i] = w[i - 1] + 2 * np.pi/n

# Time domain
t = np.arange(0, n / sampl_freq, 1/sampl_freq)
f2 = plt.figure(figsize = (10, 10))
plt.subplot(2, 1, 1)
plt.plot(t, output_signal_1)
plt.xlabel("t (in sec)")
plt.ylabel("y(t)")
plt.grid()
plt.title("With own C routine")

plt.subplot(2, 1, 2)
plt.plot(t, output_signal_2)
plt.xlabel("t (in sec)")
plt.ylabel("y(t)")
plt.title("With library function")
plt.grid()
# plt.savefig('../figs/ee18btech11039_1.eps')

# Frequency domain
f1 = plt.figure(figsize = (10, 10))
plt.subplot(2, 1, 1)
plt.plot(w, abs(Y1))
plt.xlabel("w (in rad)")
plt.ylabel("|Y(w)|")
plt.grid()
plt.title("With own C routine")

plt.subplot(2, 1, 2)
plt.plot(w, abs(Y2))
plt.xlabel("w (in rad)")
plt.ylabel("|Y(w)|")
plt.grid()
plt.title("With library function")
# plt.savefig('../figs/ee18btech11039_2.eps')
plt.show()

#If using termux
#subprocess.run(shlex.split("termux-open ../figs/ee18btech11039_1.pdf"))
#subprocess.run(shlex.split("termux-open ../figs/ee18btech11039_2.pdf"))
