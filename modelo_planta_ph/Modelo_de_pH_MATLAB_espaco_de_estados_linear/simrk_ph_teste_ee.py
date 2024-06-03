import numpy as np
import matplotlib.pyplot as plt

# Função que carrega os parâmetros do modelo contínuo (setup_pH)
def setup_pH():
    # Exemplo de parâmetros, ajuste conforme necessário
    global par, Kas, x0, Ad, Bd
    par = {'param1': 1, 'param2': 2}  # Exemplo de parâmetros
    Kas = [1, 2, 3]  # Exemplo de constantes
    x0 = np.zeros(7)  # Exemplo de estado inicial
    Ad = np.eye(7)  # Matriz de estado linearizada
    Bd = np.zeros((7, 2))  # Matriz de controle linearizada

# Função que carrega matrizes em espaço de estado discreto linearizadas (pH_linear)
def pH_linear():
    # Exemplo de matrizes, ajuste conforme necessário
    global Ad, Bd
    Ad = np.eye(7)  # Matriz de estado linearizada
    Bd = np.zeros((7, 2))  # Matriz de controle linearizada

# Função de simulação (substituto para simrk_pH)
def simrk_pH(x, Q1, Q3, h, t, par, Kas, Ts):
    # Exemplo de implementação, ajuste conforme necessário
    xc = x + h * np.random.randn(*x.shape)  # Simulação simples
    pHc = -np.log10(np.random.rand())  # pH simulado
    return xc, pHc

# Configurações iniciais
setup_pH()
pH_linear()

h = 10  # Intervalo de integração em segundos
t0 = h  # Tempo inicial em segundos
tm = 50  # Tempo de simulação em minutos
tf = 60 * tm  # Tempo de simulação em segundos
t = np.arange(t0, tf + h, h)  # Vetor de tempo "contínuo" em segundos

Ts = 40  # Tempo de amostragem
T = t[::Ts // h]  # Vetor de tempo "discreto" em segundos

Q1 = 3 * np.ones(len(T))  # Vazão de entrada ácido forte
Q3 = 2 * np.ones(len(T))  # Vazão de entrada base forte
u = np.zeros(len(T))

xc = np.zeros((len(x0), len(t)))
pHc = np.zeros(len(t))
x = np.zeros((len(x0), len(T)))
pH = np.zeros(len(T))
x[:, 0] = x0

x0l = np.zeros(7)
xl = np.zeros((len(x0), len(T)))
pHl = np.zeros(len(T))
xl[:, 0] = x0l

for k in range(1, 600 // Ts + 1):
    kc = (k - 1) * Ts // h
    x[:, k], pH[k] = simrk_pH(x[:, k - 1], Q1[k], Q3[k], h, t[kc], par, Kas, Ts)

    for i in range(7):
        xl[i, k] = Ad[i, i] * xl[i, k - 1]
    
    p = np.roots([1, (Kas[0] - (xl[1, k] + x0[1])), (Kas[0] * Kas[1] - Kas[0] * (xl[1, k] + x0[1]) - Kas[2] - Kas[0] * (xl[2, k] + x0[2])), -(Kas[0] * Kas[2] + Kas[0] * Kas[1] * (xl[1, k] + x0[1]) + 2 * Kas[0] * Kas[1] * (xl[2, k] + x0[2])), -Kas[0] * Kas[1] * Kas[2]])
    ch = np.max(np.real(p))
    pHl[k] = -np.log10(ch)

ini = k

u1 = Q1
u2 = Q3
du1 = np.zeros(len(Q3))
du2 = np.zeros(len(Q3))

du2[ini:] = 0.25 * np.ones(len(T) - ini)
u2 = Q3 + du2

for k in range(ini, len(T)):
    kc = (k - 1) * Ts // h
    x[:, k], pH[k] = simrk_pH(x[:, k - 1], u1[k], u2[k], h, t[kc], par, Kas, Ts)

    for i in range(7):
        xl[i, k] = Ad[i, i] * xl[i, k - 1] + Bd[i, 0] * du1[k] + Bd[i, 1] * du2[k]
    
    p = np.roots([1, (Kas[0] - (xl[1, k] + x0[1])), (Kas[0] * Kas[1] - Kas[0] * (xl[1, k] + x0[1]) - Kas[2] - Kas[0] * (xl[2, k] + x0[2])), -(Kas[0] * Kas[2] + Kas[0] * Kas[1] * (xl[1, k] + x0[1]) + 2 * Kas[0] * Kas[1] * (xl[2, k] + x0[2])), -Kas[0] * Kas[1] * Kas[2]])
    ch = np.max(np.real(p))
    pHl[k] = -np.log10(ch)

pplot = 1
if pplot == 1:
    plt.figure()
    plt.subplot(211)
    plt.plot((T + Ts) / 60, pH, 'k', label='pH')
    plt.plot((T + Ts) / 60, pHl, 'b*', label='pHl')
    plt.ylabel('pH')
    plt.legend()
    
    plt.subplot(212)
    plt.plot((T + Ts) / 60, u2, 'r', label='u2')
    plt.xlabel('t (min)')
    plt.ylabel('vazão de base (mL/seg)')
    plt.legend()
    plt.show()
