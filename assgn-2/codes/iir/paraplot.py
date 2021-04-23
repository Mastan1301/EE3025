import numpy as np  
import matplotlib.pyplot as plt
#If using termux
import subprocess
import shlex
#end if

N = 4

for epsilon in np.arange(0.35, 0.55+0.05, 0.05): 
    Omega1 = np.arange(0, 1, 0.01)
    H1 = 1/np.sqrt(1 + epsilon**2 * np.cos(N*np.arccos(Omega1))**2)
    Omega2 = np.arange(1,3+0.01,0.01) 
    H2 = 1/np.sqrt(1 + epsilon**2 * np.cosh(N*np.arccosh(Omega2))**2)
    H = np.concatenate((H1, H2))
    Omega = np.concatenate((Omega1, Omega2))
    plt.plot(Omega, H, label='$\epsilon$ = '+str(round(epsilon,2)))
    
plt.grid()
plt.xlabel('$\Omega$')
plt.ylabel('$|H_{a,LP}(j\Omega)|$')
plt.legend()
plt.title("Analog LowPass filter for varying epsilon")
# plt.savefig("../../figs/iir/analog_lp_epsilons.eps")
# plt.savefig("../../figs/iir/analog_lp_epsilons.pdf")
subprocess.run(shlex.split("termux-open ../figs/iir/analog_lp_epsilons.pdf"))
# plt.show()
