import numpy as np
from scipy.signal import cont2discrete
from scipy.linalg import solve_continuous_are

# Definindo os parâmetros do modelo
par = {
    1: 1.0, 2: 2.0, 3: 3.0, 4: 4.0, 5: 5.0,
    6: 6.0, 7: 7.0, 8: 8.0, 9: 9.0, 10: 10.0,
    11: 11.0, 12: 12.0, 13: 13.0, 14: 14.0, 15: 15.0,
    16: 16.0
}

# Matriz dinâmica
A = np.zeros((7, 7))

# Matriz de entrada
B = np.zeros((7, 2))

# Valores no ponto de operação
x1 = 5.4064
x2 = -3.5872e-05
x3 = 1.1120e-03
ux = 3
uy = 2

A[0, 0] = -par[9] / (2 * par[3] * np.sqrt(x1 - par[10]))
A[1, 1] = (-ux - par[1] - uy) / par[4]
A[2, 2] = (-ux - par[1] - uy) / par[4]

B[0, 0] = 1 / par[3]
B[0, 1] = 1 / par[3]
B[1, 0] = (par[11] - x2) / par[4]
B[1, 1] = (par[13] - x2) / par[4]
B[2, 0] = (par[14] - x3) / par[4]
B[2, 1] = (par[16] - x3) / par[4]
B[3, 0] = -1 / par[5]
B[5, 1] = -1 / par[7]

C = np.array([[1, 0, 0, 0, 0, 0, 0]])
D = np.array([[0, 0]])

# Sistema contínuo
sysc = (A, B, C, D)

# Discretização do sistema com tempo de amostragem de 40 segundos
Ts = 40
sysd = cont2discrete(sysc, Ts, method='zoh')
Ad = sysd[0]
Bd = sysd[1]

# Mostra as matrizes discretizadas
print("Ad = \n", Ad)
print("Bd = \n", Bd)
