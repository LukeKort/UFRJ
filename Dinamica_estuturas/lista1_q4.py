# Lista 1, questão 4. Enunciado nas anotações da lista.
# Baseado no código do Yuri

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

# Parametros

M = 1 #kg
m = 2 #kg
k = 100 #N/m
c = 0.4 #N.s/m
L = 1 #m
g = 9.81 #m/s²

# Condição inicial

u0 = 0 #m
v0 = 0 #m/s
theta0 = np.pi/4 #rad
omega0 = 0 #rad/s

save_charts = 'y'

S_0 = np.array([u0, v0, theta0, omega0])

def dS_dt(t, S):
    
  u, v, theta, omega = S

  alpha = (np.cos(theta)/(L*((2*(M+m) - m*((np.cos(theta))**2)))))*(- m*L*(omega**2)*np.sin(theta) + k*u + c*v - g*np.tan(theta)*(M+m))

  a = (1/(M + m))*(- m*L*alpha*((np.cos(theta))) + m*L*omega**2*np.sin(theta) - k*u - c*v)

  return [v, a, omega, alpha]

t_max = 100
n_points = 10000
period = t_max/n_points
fps = 1/period

t_span = np.linspace(0, t_max, n_points)

sol = odeint(dS_dt, y0=S_0, t=t_span, tfirst=True) 

u = sol.T[0]
v = sol.T[1]
theta = sol.T[2]
omega = sol.T[3]

fig_size = (8,4)

# Velocidade linear -------------------------------------------------------------------------------------------------------------------

title = 'u x t'
plt.figure(title,figsize=fig_size)
plt.plot(t_span, u, '-', linewidth=1.2, color = 'blue')
plt.title(title)
plt.ylabel('$u$ [m]', fontsize=12)
plt.xlabel('$t$ [s]', fontsize=12)
# plt.grid(True)
# plt.show()

if save_charts == 'y':
    
    import matplotlib.pyplot as plt
    import os

    from time import strftime
    folder_name = 'Charts' + ' ' + strftime("%Y-%m-%d %H-%M")

    if not os.path.exists(folder_name):
        os.makedirs(folder_name)

    plt.savefig(folder_name + '\\' + title + '.png')

# Velocidade linear vs deslocamento----------------------------------------------------------------------------------------------------

title = 'u x v'
plt.figure(title,figsize=fig_size)
plt.plot(u, v, '-', linewidth=1.2, color = 'blue')
plt.title(title)
plt.ylabel('$v$ [m/s]', fontsize=12)
plt.xlabel('$u$ [m]', fontsize=12)
# plt.grid(True)
# plt.show()

if save_charts == 'y':
    plt.savefig(folder_name + '\\' + title + '.png')


# Velocidade angular-------------------------------------------------------------------------------------------------------------------

title = 'theta x t'
plt.figure(title,figsize=fig_size)
plt.plot(t_span, theta, '-', linewidth=1.2, color = 'blue')
plt.title('$\\theta$ x t')
plt.ylabel('$\\theta$ [rad]', fontsize=12)
plt.xlabel('t [s]', fontsize=12)
# plt.grid(True)
# plt.show()

if save_charts == 'y':
    plt.savefig(folder_name + '\\' + title + '.png')

# Velocidade angular vs angulo---------------------------------------------------------------------------------------------------------

title = 'theta x omega'
plt.figure(title,figsize=fig_size)
plt.plot(theta, omega, '-', linewidth=1.2, color = 'blue')
plt.title('$\\theta$ x $\\omega$')
plt.ylabel('$\\omega$ [rad/s]', fontsize=12)
plt.xlabel('$\\theta$ [rad]', fontsize=12)
# plt.grid(True)
# plt.show()

if save_charts == 'y':
    plt.savefig(folder_name + '\\' + title + '.png')

# Energia mecânica----------------------------------------------------------------------------------------------------------------------

Mec_energy = np.zeros(len(t_span))
for i in range(len(t_span)):
  Mec_energy[i] = (1/2)*(M + m)*v[i]**2 + m*L*omega[i]*(v[i]*np.cos(theta[i]) + L*omega[i]) + (1/2)*k*(u[i]**2) + (m*g*L*( 1 - np.cos(theta[i])))
#   Mec_energy[i] = (1/2)*m*(v[i]**2 + 2*v[i]*omega[i]*L*np.cos(theta[i]) + (omega[i]**2)*(L**2)) + (1/2)*M*(v[i]**2) + (1/2)*k*(u[i]**2) + (m*g*L - m*g*L*np.cos(theta[i]))

title = 'Energia Mecânica [J]'
plt.figure(title,figsize=fig_size)
plt.plot(t_span, Mec_energy, '-', linewidth=1.2, color = 'blue')
plt.ylabel(title, fontsize=12)
plt.xlabel('$t$ [s]', fontsize=12)
# plt.grid(True)
# plt.show()

if save_charts == 'y':
    plt.savefig(folder_name + '\\' + title + '.png')

# Trajetórias---------------------------------------------------------------------------------------------------------------------------

m_x_path = np.zeros(len(t_span))
m_y_path = np.zeros(len(t_span))
M_x_path = np.zeros(len(t_span))
M_y_path = np.zeros(len(t_span))

for i in range(len(t_span)):
  m_x_path[i] = L*np.sin(theta[i]) + u[i]
  m_y_path[i] = -L*np.cos(theta[i])
  M_x_path[i] = u[i]

title = 'Trajetórias das massas M e m'
plt.figure(title,figsize=(8,6))

plt.plot(m_x_path, m_y_path, '-', linewidth=1.0, color = 'blue', label='$m$')
plt.plot(M_x_path, M_y_path, '-', linewidth=1.8, color = 'green', label='$M$')
plt.title(title, fontsize=12)
plt.ylim((-1.1,0.1))
plt.ylabel('y [m]', fontsize=12)
plt.xlabel('x [m]', fontsize=12)
# plt.grid(True)
plt.legend(fontsize=12)

if save_charts == 'y':
    plt.savefig(folder_name + '\\' + title + '.png')

plt.show()