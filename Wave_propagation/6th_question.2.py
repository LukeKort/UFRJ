import numpy as np
import matplotlib.pyplot as plt

rho = 2700  # Bar's density
E = 70e9 # Bar's Young Modulos
# rho = 7850  # Bar's density
# E = 210e9 # Bar's Young Modulos


# Fist case - 1 element bar

L = 0.3 # Bar length

r = 0.0125 # Bar's radio

A = np.pi*(r**2) # Area

rho = rho*A

c_0 = np.sqrt(E/rho) # Bar's propagation velocity

M = rho*A*L # Bar's mass

I = 0.5*M*r**2 # Moment of inertia

N = 1000

# Elementary length
L_1 = 0.150
L_2 = 0.03
L_3 = 0.120

#Elementay Young Modulus 
E_1 = E
E_2 = 0.5*E
E_3 = E

#Elementay C_0
c_1 = np.sqrt(E_1/rho)
c_2 = np.sqrt(E_2/rho)
c_3 = np.sqrt(E_3/rho)

def alpha(L, E, c, omega):
  k = omega/c
  aux = E*A/L * ((1j*k*L)/(1-np.exp(-2*1j*k*L)))
  return aux

def beta(L, c, omega):
  k = omega/c
  aux = 1 + np.exp(-2*1j*k*L)
  return aux

def gamma(L, c, omega):
  k = omega/c
  aux = -2*np.exp(-1j*k*L)
  return aux

# Question 6A

# def H(w): # Calculates H for given angular frequency

#     k = w/c_0 # wave number

#     h = 1/((E*A*1j*k/(1-np.exp(-1j*2*k*L)))*(1+np.exp(-1j*2*k*L)))

#     # h = (1-np.exp(-1j*2*k*L))*(1+np.exp(-1j*2*k*L)))/(E*A*1j*k)

#     return h

# Question 6B

def H(w): # Calculates H for given angular frequency



    beta0 = beta(L_1, c_1, w)
    alpha0 = alpha(L_1, E_1, c_1, w)
    gamma0 = gamma(L_1, c_1, w)

    beta1 = beta(L_2, c_2, w)
    alpha1 = alpha(L_2, E_2, c_2, w)
    gamma1 = gamma(L_2, c_2, w)

    beta2 = beta(L_3, c_3, w)
    alpha2 = alpha(L_3, E_3, c_3, w)
    gamma2 = gamma(L_3, c_3, w)

    h = (alpha0*alpha1*beta0*beta1 + alpha0*alpha2*beta0*beta2 + alpha1**2*beta1**2 - alpha1**2*gamma1**2 + alpha1*alpha2*beta1*beta2)/(alpha0*alpha1*alpha2*beta0*beta1*beta2 + alpha0*alpha2**2*beta0*beta2**2 - alpha0*alpha2**2*beta0*gamma2**2 + alpha1**2*alpha2*beta1**2*beta2 - alpha1**2*alpha2*beta2*gamma1**2 + alpha1*alpha2**2*beta1*beta2**2 - alpha1*alpha2**2*beta1*gamma2**2)


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

# w_plus = np.linspace(w[0]*(1-6.87e-4), w[-1]*(1+6.87e-4),N)
w_plus = np.linspace(w[0]*(1-9.7e-1), w[-1]*(1+6.87e-4),N)

for i in range(len(H_w_plus)):
    H_w_plus[i] = np.abs(H(w_plus[i]))

title = str('1st to 4th Natural Frequency')
# title = str('1001th to 1004th Natural Frequency')
plt.figure(title)
plt.plot(w_plus, H_w_plus)
# plt.scatter(w, H_w, color = 'green')
plt.vlines(w, ymin=0, ymax=10E-7, color = 'red', linestyles='dotted', label = 'Natural Frequencies')
plt.yscale('log')
plt.xlabel('\u03C9 [rads/s]')
plt.ylabel('H(\u03C9) [m/N]')
plt.legend()
plt.title(title)
plt.show()