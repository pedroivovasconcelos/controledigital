import numpy as np
from scipy.signal import cont2discrete
from modelo_planta_ph.setup_pH import par  # Importa variáveis do setup_pH

# Matriz dinâmica
A = np.zeros((7, 7))

# Matriz de entrada
B = np.zeros((7, 2))

# Valores no ponto de operação
x1 = 5.4064
x2 = -3.5872e-05
x3 = 1.1120e-03
ux = 3 # acido forte
uy = 2 # base forte

A[0, 0] = -par[8] / (2 * par[2] * np.sqrt(x1 - par[9]))
A[1, 1] = (-ux - par[0] - uy) / par[3]
A[2, 2] = (-ux - par[0] - uy) / par[3]

B[0, 0] = 1 / par[2]
B[0, 1] = 1 / par[2]
B[1, 0] = (par[10] - x2) / par[3]
B[1, 1] = (par[12] - x2) / par[3]
B[2, 0] = (par[13] - x3) / par[3]
B[2, 1] = (par[15] - x3) / par[3]
B[3, 0] = -1 / par[4]
B[5, 1] = -1 / par[6]

C = np.array([[1, 0, 0, 0, 0, 0, 0]])
D = np.array([[0, 0]])

# Sistema contínuo
sysc = (A, B, C, D)

# Discretização do sistema com tempo de amostragem de 40 segundos
Ts = 40
sysd = cont2discrete(sysc, Ts, method='zoh')
Ad = sysd[0]
Bd = sysd[1]

# Mostra as matrizes discretizadas (opcional)
print("Ad = \n", Ad)
print("Bd = \n", Bd)
