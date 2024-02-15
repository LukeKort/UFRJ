import matplotlib.pyplot as plt # Importar bliblioteca para plotagem de gráficos
import numpy as np # Importar bliblioteca para realização de cálculos matemáticos complexos

c_1 = 0.118 # Valor de entrada para a constante do material C1
c_2 = 0.474 # Valor de entrada para a constante do material C2


T11 = np.zeros(30) # Pré-alocagem do vetor T11 - ganho de perfomance
lambda_1 = np.ones(30) # Pré-alocagem do vetor T11 - ganho de perfomance

for i in range(1,30): # Função para gerar o lambda em forma de vetor 
    lambda_1[i] = lambda_1[i-1] + 0.1 # Os diferentes lambdas que serâo testados tem incremento de 0,1

for i in range(30): # Função para cálculo da componente (1,1) do tensor T

    lambda_3 = lambda_1[i]**(-2) # Cálculo do lambda_3 em função do lambda 1

    alpha_0 = 2*c_1 # Cálculo do alpha 0 e alpha 1 em função das constantes do material pela formulação de Mooney-Rivlin
    alpha_1 = -2*c_2

    p = alpha_1*lambda_1[i]**(-2)+alpha_1*lambda_1[i]**2 # Cálculo da carga de estiramento

    B = np.array([[lambda_1[i]**2, 0 , 0],[0, 4, 0],[0, 0, lambda_3]]) # Cálculo da matriz B
    B_inv = np.linalg.inv(B) # # Cálculo da matriz inversa de B

    T11[i] = p*1 + alpha_0*B[0][0] + i*B_inv[0][0] # Cálculo da componente (1,1) do tensor T

    # O T11 é calculado em forma de vetor, que armazena diversos valores, necessário para plotar o gráfico.

# Plot
plt.figure('Gráfico 2') # Cria uma figura para armazenar o gráfico
plt.title('Tensor vs Taxa de estiramento') # Título do gráfico
plt.xlabel('Taxa de estiramento \u03BB') # Rótulo do eixo x
plt.ylabel('Tensão') # Rótulo do eixo x
g_label = 'Mooney: \u03B1\u2080 =' + str(alpha_0) + ' \u03B1\u2081 =' + str(alpha_1) +'\nC\u2081 =' + str(c_1) + ' C\u2082 =' + str(c_2) # Texto da legenda
plt.plot(lambda_1,T11, label = g_label ) # Função para plotar o gráfico
plt.legend() # Função de chamamento da legenda
plt.show()  # Função de chamamento do gráfico