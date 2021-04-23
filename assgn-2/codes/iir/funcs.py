import numpy as np

# This fuctions adds two polynomials defined by vectors x and y 
def add(x, y):
    m, n = len(x), len(y)
    if(m == n):
        return x + y

    z = np.zeros(abs(m - n))
    if m < n:
        res = np.concatenate((z, x), axis = 0) + y
    else:
        res = x + np.concatenate((z, y), axis = 0)
    return res


def polypower(v, N):
    y = np.array([1])
    if N > 0:
        for i in range(int(N)):
            y = np.convolve(y, v)
    return y

# This function transforms the stable bandpass  filter obtained
# from the Chebyschev approximation to the equivalent bandpass
# digital filter using the bilinear transformation
# [dignum,digden,G_bp] = bilin(p,om)
# H_bp(s) = G_bp*num(s)/den(s) is the analog bandpass filter
# obtained through the Chebyschev filter design
# H(z) = G*dignum(z)/digden(z) is digital bandpass filter
# obtained from H_bp(s) by substituting s = (z-1)/(z+1)
# G is obtained from the condition H(om) = 1

def bilin(p, om):
    N = len(p)
    const = np.array([-1, 1])
    v = np.array([1])
    if N > 2:            
        for i in range(1, N):
            v = np.convolve(v, const)
            v = add(v, p[i] * polypower(np.array([1, 1]), i))
        
        digden = v

    elif N == 2:
        M = len(v)
        v[M - 2] += p[N - 1]
        v[M - 1] += p[N - 1]
        digden = v

    else:
        digden = p

    dignum = polypower(np.array([-1, 0, 1]), int((N - 1) / 2))
    G_bp = abs(np.polyval(digden, np.exp(-1j * om)) / np.polyval(dignum, np.exp(-1j * om)))
    
    return (dignum, digden, G_bp)

# THIS FUNCTION GENERATES THE CHEBYSCHEV POLYNOMIAL COEFFICIENTS OF ORDER N
def cheb(N):
    v = np.array([1, 0])
    u = np.array([1])

    if N == 0:
        w = u
    elif N == 1:
        w = v
    else:
        for i in range(N - 1):
            p = np.convolve(np.array([2, 0]), v)
            m, n = len(p), len(u)
            w = p + np.concatenate((np.zeros(m - n), u), axis = 0)
            u = v
            v = w

    return w

# The function
# [C K] = lattice(N,D)
# computes the lattice parameters K and the ladder parameters C for the
# system function H(z) = N(z)/D(z), where both the numerator and denominator
# are of the same order
def lattice(c, v):
    u = np.fliplr(v)
    m = len(v)
    K[m - 2] = v[m - 1]
    C[m - 1] = c[m - 1]
    while m > 2 and K[m - 2] != 1:
        c = c - C[m - 1] * u
        v = (v - K[m - 2] * u) / (1 - K[m - 2]**2)
        m -= 1
        v = v[0 : m]
        c = c[0 : m]
        u = np.fliplr(v)

        if m > 2:
            K[m - 2] = v[m - 1]

        C[m - 1] = c[m - 1]

    return (C, K)

# This function gives the low pass stable filter
# for the Chebyschev approximation based upon
# the design parameters epsilon and N
# H(s) = G/p(s)
# [p,G] = lp_stable_cheb(epsilon,N)
def lp_stable_cheb(epsilon, N):
    # Analytically obtaining the roots of the Chebyschev polynomial
    # in the left half of the complex plane

    beta = ((np.sqrt(1 + epsilon**2)+ 1)/epsilon) ** (1.0 / N)
    r1 = (beta**2 - 1)/(2 * beta)
    r2 = (beta**2 + 1)/(2 * beta)

    # Obtaining the polynomial approximation for the low pass
    # Chebyschev filter to obtain a stable filter
    u = np.array([1])
    for n in range(int((N / 2))):
        phi = np.pi/2 + ((2 * n + 1) * np.pi / (2 * N))
        v = np.array([1, - 2 * r1 * np.cos(phi), r1 ** 2 * np.cos(phi)**2 + r2**2 * np.sin(phi)**2])
        p = np.convolve(v, u)
        u = p
    
    # Evaluating the gain of the stable lowpass filter
    # The gain has to be 1/sqrt(1+epsilon^2) at Omega = 1
    G = abs(np.polyval(p, 1j))/np.sqrt(1 + epsilon**2)

    return (p, G)

# This function transforms the lowpass stable filter obtained
# from the Chebyschev approximation to the bandpass
# equivalent
# [num,den,G] = lpbp(p,Omega0,B,Omega_p2)
# Omega0 and B are the lowpass-bandpass transformation parameters
# and Omega_p2 is the lower limit of the passband, used
# to evaluate the gain G_bp
# H(s) = G/p(s) is the stable low pass Cheybyschev approximation
# Hbp(s) = G_bp*num(s)/den(s) is the corresponding bandpass stable
# filter
def lpbp(p, Omega0, B, Omega_p2):
    N = len(p)
    const = np.array([1, 0, Omega0 ** 2])
    v = const.copy()

    if N > 2:
        for i in range(1, N):
            M = len(v)
            v[M - i - 1] = v[M-i-1] + (p[i])*(B**i)

            if i < N - 1:
                v = np.convolve(const, v)
            
        den = v
    
    elif N == 2:
        M = len(v)
        v[M - 2] += p[N - 1] * B
        den = v
    
    else:
        den = p

    num = np.concatenate((np.array([1]),np.zeros(N-1)))          
    G_bp = abs(np.polyval(den, 1j*Omega_p2)/(np.polyval(num, 1j*Omega_p2)))
    
    return (num, den, G_bp)

    
