import matplotlib.pyplot as plt # Importar bliblioteca para plotagem de gráficos
import numpy as np # Importar bliblioteca para realização de cálculos matemáticos complexos

c_1 = 38.845 # Valor de entrada para a constante do material C1
c_2 = 77.69 # Valor de entrada para a constante do material C2

# # Para cálculo do C1/C2 apenas
# n = 20
# c_1 = n*38.845 # Valor de entrada para a constante do material C1
# c_2 = 77.690 # Valor de entrada para a constante do material C2


T11_MR = np.zeros(30) # Pré-alocagem do vetor T11 - ganho de perfomance
T11_nH = np.zeros(30) # Pré-alocagem do vetor T11 - ganho de perfomance
T33_MR = np.zeros(30) # Pré-alocagem do vetor T11 - ganho de perfomance
lambda_1 = np.ones(30) # Pré-alocagem do vetor T11 - ganho de perfomance

for i in range(1,30): # Função para gerar o lambda em forma de vetor 
    lambda_1[i] = lambda_1[i-1] + 0.1 # Os diferentes lambdas que serâo testados tem incremento de 0,1

for i in range(30): # Função para cálculo da componente (1,1) do tensor T

    lambda_3 = lambda_1[i]**(-1) # Cálculo do lambda_3 em função do lambda 1

    alpha_0 = 2*c_1 # Cálculo do alpha 0 e alpha 1 em função das constantes do material pela formulação de Mooney-Rivlin
    alpha_1 = -2*c_2

    p = alpha_0*lambda_1[i]**(-2)+alpha_1*lambda_1[i]**2 # Cálculo da carga de estiramento

    B = np.array([[lambda_1[i]**2, 0 , 0],[0, 1, 0],[0, 0, lambda_3**2]]) # Cálculo da matriz B
    B_inv = np.linalg.inv(B) # # Cálculo da matriz inversa de B

    T11_MR[i] = -p*1 + alpha_0*B[0][0] + alpha_1*B_inv[0][0] # Cálculo da componente (1,1) do tensor T por Mooney - Rivlin
    T11_nH[i] = alpha_0*(B[0][0] - B_inv[0][0]) # Cálculo da componente (1,1) do tensor T por neo-Hooken
    T33_MR[i] = -p*1 + alpha_0*B_inv[0][0] + alpha_1*B[0][0] # Cálculo da componente (1,1) do tensor T por Mooney - Rivlin
    # O T11 é calculado em forma de vetor, que armazena diversos valores, necessário para plotar o gráfico.

# Plot Gráfico 2D Mooney vs ne-Hookean
plt.figure('Gráfico Tensão vs Estiramento') # Cria uma figura para armazenar o gráfico
plt.title('Tensão vs Estiramento') # Título do gráfico
plt.xlabel('\u03BB') # Rótulo do eixo x
plt.ylabel('$T_{11}$') # Rótulo do eixo x
g_label = 'Mooney: \u03B1\u2080 =' + str(alpha_0) + ' \u03B1\u2081 =' + str(alpha_1) +'\nC\u2081 =' + str(c_1) + ' C\u2082 =' + str(c_2) # Texto da legenda
plt.plot(lambda_1,T11_MR, 'b', lambda_1,T11_nH, 'b--') # Função para plotar o gráfico
# plt.axis([1,4,0,3500])
plt.legend([
    'Mooney: \u03B1\u2080 =' + str(alpha_0) + ' \u03B1\u2081 =' + str(alpha_1) +'\nC\u2081 =' + str(c_1) + ' C\u2082 =' + str(c_2), 
    'neo-Hookean'
    ]) # Função de chamamento da legenda
plt.grid('on')
# plt.show()  # Função de chamamento do gráfico

# Para plotar gráfico C1/C2 vs Lambda apenas

c1vc2 = np.zeros(30) # Pré-alocagem do vetor T11 - ganho de perfomance

for i in range(30): # Função para cálculo da componente (1,1) do tensor T

    c_1_ = (i+1)*c_1 # Valor de entrada para a constante do material C1
    c_2_ = c_2 # Valor de entrada para a constante do material C2

    c1vc2[i] = c_1_/c_2_

    lambda_1_ = 4
    lambda_3_ = 1/4

    alpha_0 = 2*c_1_ # Cálculo do alpha 0 e alpha 1 em função das constantes do material pela formulação de Mooney-Rivlin
    alpha_1 = -2*c_2_

    p = alpha_0*lambda_1_**(-2)+alpha_1*lambda_1_**2 # Cálculo da carga de estiramento

    B = np.array([[lambda_1_**2, 0 , 0],[0, 1, 0],[0, 0, lambda_3_**2]]) # Cálculo da matriz B
    B_inv = np.linalg.inv(B) # # Cálculo da matriz inversa de B

    T11_MR[i] = -p*1 + alpha_0*B[0][0] + alpha_1*B_inv[0][0] # Cálculo da componente (1,1) do tensor T por Mooney - Rivlin
    T11_nH[i] = alpha_0*(B[0][0] - B_inv[0][0]) # Cálculo da componente (1,1) do tensor T por neo-Hooken
    # O T11 é calculado em forma de vetor, que armazena diversos valores, necessário para plotar o gráfico.

# Plot Gráfico razão C1/C2
plt.figure('Gráfico razão C1/C2') # Cria uma figura para armazenar o gráfico
plt.title('Tensão vs $C_{1}$/$C_{2}$') # Título do gráfico
plt.xlabel('$C_{1}$/$C_{2}$ \t \u03BB = 4') # Rótulo do eixo x
plt.ylabel('$T_{11}$') # Rótulo do eixo x
plt.plot(c1vc2,T11_MR, 'b', c1vc2,T11_nH, 'b--') # Função para plotar o gráfico
plt.legend([
    'Mooney', 
    'neo-Hookean'
    ]) # Função de chamamento da legenda
plt.grid('on')
#plt.show()  # Função de chamamento do gráfico

Tdiv = ((T11_MR - T11_nH)/T11_nH)*100

# Plot Gráfico erro Mooney vs neo-Hookean
plt.figure('Gráfico do erro Mooney vs neo-Hookean') # Cria uma figura para armazenar o gráfico
plt.title('Erro Mooney-neo-Hookean vs Tensão') # Título do gráfico
plt.xlabel('$T_{11}$ neo-Hookean \t \u03BB = 4') # Rótulo do eixo x
plt.ylabel('Erro %') # Rótulo do eixo x
plt.plot(T11_nH,Tdiv, 'b') # Função para plotar o gráfico
plt.grid('on')
plt.show()  # Função de chamamento do gráfico