import numpy as np

# Definição de parâmetros para o modelo de pH
# Deve ser rodado uma vez no início de cada simulação

ka1 = 10**-6.35  # ka1=[HCO3-][H+]/[H2CO3-]
ka2 = 10**-10.33 # ka2=[CO3^2-][H+]/[HCO3-]
kw  = 10**-14    # kw1=[H+][OH-]
Kas = [ka1, ka2, kw]

# Condições Iniciais de H, Wa, Wb, Hta, Htt, Htb e Htc
x0 = np.array([12.3, 10**-14, 0, 45.45, 36.37, 45.45, 0])

#--------------------------------------------------------------------------
# Constantes do modelo Q1, Q2, Q3, Q4, A, c, h0, wa1, wa2, wa3, wb1, wb2, wb3
#--------------------------------------------------------------------------

Q2 = 0.1  # Vazão de tampão em mL/seg
Q4 = 0    # Vazão de saída quando o nível é constante

Ar = 72    # Área da base do reator
Vr = 870
Ata = 1320  # Área da base do tanque de reagente ácido (cm^2)
Att = 1320  # Área da base do tanque de reagente tampão (cm^2)
Atb = 1320  # Área da base do tanque de reagente base (cm^2)
Atc = 1963.5  # Área da base do tanque de coletor (rejeito) (cm^2)
c = 8      # Constante da válvula de saída (Q4)
h0 = 5     # Altura para diferencial de pressão saída do tanque
wa1 = 0.006    # Concentração ácida em Q1 -> [HNO3]
wa2 = -0.06    # Concentração ácida em Q2 -> -[NaHCO3]
wa3 = -0.0061  # Concentração ácida em Q3 -> -[NaHCO3]-[NaOH]
wb1 = 0        # Concentração base em Q1
wb2 = 0.06     # Concentração base em Q2 -> [NaHCO3]
wb3 = 0.0001   # Concentração base em Q3 -> [NaHCO3]
par = [Q2, Q4, Ar, Vr, Ata, Att, Atb, Atc, c, h0, wa1, wa2, wa3, wb1, wb2, wb3]
