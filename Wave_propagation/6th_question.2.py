import numpy as np
import matplotlib.pyplot as plt

rho = 2700  # Bar's density
E = 70e9 # Bar's Young Modulos

L = 0.3 # Bar length
r = 0.0125 # Bar's radio
A = np.pi*(r**2) # Area
rho = rho*A
c_0 = np.sqrt(E/rho) # Bar's propagation velocity

P = 100 # Load at right end
M = rho*A*L # Bar's mass
I = 0.5*M*r**2 # Moment of inertia

N = 1000
t = np.linspace(0, 1, N)
f = np.fft.fftfreq(N, t[1] - t[0])

u = P/(E*A*L)*np.exp(-1j*c_0*f*L)/(1 - np.exp(-1j*c_0*f*L))

y = u*(1 - np.exp(-1j*c_0*f*L))/(1 + np.exp(-1j*c_0*f*L))

FRF = y/u

plt.plot(f, np.abs(FRF))
plt.xlabel('Frequency [Hz]')
plt.yscale('log')
plt.ylabel('FRF')
plt.show()
