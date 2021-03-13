import soundfile as sf
import numpy as np
from scipy import signal
import matplotlib.pyplot as plt

#If using termux
import subprocess
import shlex

#Reading the soundfile 
input_signal, fs = sf.read('./Sound_Noise.wav')
sampl_freq = fs
order = 4
cutoff_freq = 4000
Wn = 2 * cutoff_freq / sampl_freq
n = len(input_signal)
print(n)
s = 2 ** np.ceil(np.log2(n))
while n < s: # Making the file size as 2 ^ k.
	np.append(input_signal, 0)
	n += 1

sf.write('./Sound_Noise.wav', input_signal, sampl_freq)

#Passing butterworth filter
b, a = signal.butter(order, Wn, 'low')

# Computing H(z)
w = 2 * np.pi * np.arange(n)/n
z = np.exp(-1j * w)
H = np.polyval(b, z)/np.polyval(a, z)

X = np.fft.fft(input_signal)
Y1 = H * X
output_signal_1 = np.fft.ifft(Y1).real
output_signal_2 = signal.filtfilt(b, a, input_signal)
Y2 = np.fft.fft(output_signal_2)
Y1 = np.fft.fftshift(Y1)
Y2 = np.fft.fftshift(Y2)

sf.write('./Sound_With_ReducedNoise_1.wav', output_signal_1, sampl_freq)
sf.write('./Sound_With_ReducedNoise_2.wav', output_signal_2, sampl_freq)

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
plt.title("With own routine")

plt.subplot(2, 1, 2)
plt.plot(t, output_signal_2)
plt.xlabel("t (in sec)")
plt.ylabel("y(t)")
plt.title("With library function")
plt.grid()
plt.savefig('../figs/ee18btech11039_1.eps')

# Frequency domain
f1 = plt.figure(figsize = (10, 10))
plt.subplot(2, 1, 1)
plt.plot(w, abs(Y1))
plt.xlabel("w (in rad)")
plt.ylabel("|Y(w)|")
plt.grid()
plt.title("With own routine")

plt.subplot(2, 1, 2)
plt.plot(w, abs(Y2))
plt.xlabel("w (in rad)")
plt.ylabel("|Y(w)|")
plt.grid()
plt.title("With library function")
plt.savefig('../figs/ee18btech11039_2.eps')
plt.show()

#If using termux
#subprocess.run(shlex.split("termux-open ../figs/ee18btech11039_1.pdf"))
#subprocess.run(shlex.split("termux-open ../figs/ee18btech11039_2.pdf"))
