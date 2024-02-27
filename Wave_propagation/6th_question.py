import numpy as np
import matplotlib.pyplot as plt

rho = 2700  # Bar's density
E = 70e9 # Bar's Young Modulos
# c_0 = np.sqrt(E/(rho)) # Bar's propagation velocity

# Fist case - 1 element bar

L = 0.3 # Bar length

r = 0.0125 # Bar's radio

A = np.pi*(r**2) # Area

rho = rho*A

c_0 = np.sqrt(E/rho) # Bar's propagation velocity


P = 1 # Load at right end

print('teste',P*L/(A*E))

M = rho*A*L # Bar's mass

I = 0.5*M*r**2 # Moment of inertia

N = 1000

def u_2(w): # Calculates u_2 for given angular frequency

    k = w/c_0 # wave number

    u = 1/((E*A*1j*k/(1-np.exp(-1j*2*k*L)))*(1+np.exp(-1j*2*k*L)))

    return u

def H(w): # Calculates H for given angular frequency

    k = w/c_0 # wave number

    h = 1/((E*A*1j*k/(1-np.exp(-1j*2*k*L)))*(1+np.exp(-1j*2*k*L)))

    # h = (1-np.exp(-1j*2*k*L))*(1+np.exp(-1j*2*k*L)))/(E*A*1j*k)

    return h

H_w = np.zeros(4) # Stores H_2 for given u_w omega
H_w_plus = np.zeros(N)


w = np.zeros(4) # Stores the frequencies
w_plus = np.zeros(N) # Stores the frequencies

def w_n(n): # Calculates the n natural frequency

    # w = (n/(2*L))*np.sqrt(E/rho)*2*np.pi

    w = (2*n - 1)*np.pi*c_0/(2*L)

    return w

# n = [1001,1002,1003,1004]
n = [1,2,3,4]

for i in range(4):
    w[i] = w_n(n[i])

print('Omega:', w)

for i in range(4):
    H_w[i] = np.abs(H(w[i]))

w_plus = np.linspace(w[0],w[-1],N)

for i in range(len(H_w_plus)):
    H_w_plus[i] = np.abs(H(w_plus[i]))

title = str('1st to 4th Natural Frequency')
# title = str('1001th to 1004th Natural Frequency')
plt.figure(title)
plt.plot(w_plus, H_w_plus)
# plt.scatter(w, H_w, color = 'green')
plt.vlines(w, ymin=0, ymax=10E8, color = 'red', linestyles='dotted', label = 'Natural Frequencies')
plt.yscale('log')
plt.xlabel('\u03C9 [rads/s]')
plt.ylabel('H(\u03C9) [m/N]')
plt.legend()
plt.title(title)
plt.show()