# -*- coding: utf-8 -*-

# Pacotes utilizados
import numpy as np
import matplotlib.pyplot as plt

#**************************************************************************
#* Simulacao de planta de neutralizacao de pH  - com entrada PRBS -       *
#* definida em termos de equacoes diferenciais e controle em termos de    *
#* equacao de diferencas. Portanto a malha eh hibrida, como em um caso    *
#* real. 
#
# Código original em Matlab: LAA 04/08/22
# Traduzido para Python: PASB 29/08/22

#%% Definição dos parâmetros para o modelo de pH 
#Código orignal: setup_pH.m


ka1 = 10**-6.35  # ka1=[HCO3-][H+]/[H2CO3-]
ka2 = 10**-10.33 # ka2=[CO3^2-][H+]/[HCO3-]
kw  = 10**-14    # kw1=[H+][OH-]
Kas = np.array([ka1, ka2, kw])

# Condicoes Iniciais de H, Wa, Wb, Hta, Htt, Htb e Htc
x0 = np.array([12.3, 10**-14, 0, 45.45, 36.37, 45.45, 0])

#--------------------------------------------------------------------------
# Constantes do modelo Q1, Q2, Q3, Q4, A, c, h0, wa1, wa2, wa3, wb1, wb2, wb3
#--------------------------------------------------------------------------

Q2 = 0.1 #Vazao de tampao em mL/seg
# Q4 = Q1+Q2+Q3; # Vazao de saida quando o nivel eh constante
Q4 = 0

Ar = 72 # Area da base do reator
Vr = 870
Ata = 1320 # Area da base do tanque de reagente acido(cm^2) 1320
Att = 1320 # Area da base do tanque de reagente tampao(cm^2)
Atb = 1320 # Area da base do tanque de reagente base(cm^2)
Atc = 1963.5 # Area da base do tanque de coletor (rejeito) (cm^2)1963.5
c = 8 # Constante da valvula de saida(Q4)
h0 = 5 # Altura para diferencial de pressao saida do tanque
wa1 = 0.006 # concentracao acida em Q1 -> [HNO3] 
wa2 = -0.06 # concentracao acida em Q2 -> -[NaHCO3]
wa3 = -0.0061 # concentracao acida em Q3 -> -[NaHCO3]-[NaOH]
wb1 = 0 # concentracao base em Q1
wb2 = 0.06 # concentracao base em Q2 -> [NaHCO3]
wb3 = 0.0001 # concentracao base em Q3 -> [NaHCO3]
par = np.array([Q2, Q4, Ar, Vr, Ata, Att, Atb, Atc, c, h0, wa1, wa2, wa3, wb1, wb2, wb3])

#%%

#Função rkpH
#Código orignal: rkpH.m
def rkpH(x0,ux,uy,h,t,par):
    #function [x,xc] = rkpH(x0,ux,uy,h,t,par)
    # function x = rkpH(x0,ux,uy,h,t,par%
    # x0 eh a condicao inicial para a iteracao em questao
    # ux vazao de entrada acido forte (mL/seg)
    # xy vazao de entrada base forte (mL/seg)
    # h eh o intervalo de integracao (nao deve ser alterado pelo usuario)
    # t eh o instante de tempo atual
    # par eh o vetor de parametros do modelo 
    #
    # x eh o vetor de estado em tempo discreto
    # xc eh o vetor de estado em tempo continuo
    
    # 1st evaluation at 
    xd = dvpH(x0,ux,uy,t,par)
    savex0 = x0
    phi = xd
    x0 = savex0+0.5*h*xd
    
    # 2nd evaluation 
    xd = dvpH(x0,ux,uy,t+0.5*h,par)
    phi = phi+2*xd
    x0 = savex0+0.5*h*xd
    
    # 3rd evaluation
    xd = dvpH(x0,ux,uy,t+0.5*h,par)
    phi = phi+2*xd
    x0 = savex0+h*xd
    
    # 4th evaluation
    xd = dvpH(x0,ux,uy,t+h,par)
    x = savex0+(phi+xd)*h/6
    
    return x

#Função dvpH
#Código orignal: dvpH.m
def dvpH(x,ux,uy,t,par):
    # function xdot = dvpH(x,ux,uy,t,par)
    # xdot eh o vetor de derivadas (campo vetorial) usado pelo algoritmo de
    # integracao 
    # ux e uy sao os valores das entradas mantidas constantes durante o
    # intervalor de amostragem
    # t eh o instante de tempo atual. Usado apenas como registro do instante a
    # que se refere a avaliacao do campo vetorial. Normalmente nao eh usado.
    # par eh o vetor de parametros do modelo
    
    #-------------------------------------------------------------------------%
    # Equacoes diferenciais                                                   %
    #-------------------------------------------------------------------------%
    # x(1) Nivel do Reator
    # x(2) Concentracao de acido
    # x(3) Concentracao de Base
    # x(4) Nivel Tanque de acido
    # x(5) Nivel Tanque de Tampao
    # x(6) Nivel Tanque de Base
    # x(7) Nivel Tanque do Coletor
    # ux vazao de entrada acido forte (mL/seg)
    # xy vazao de entrada base forte (mL/seg)
    
    xd=np.zeros(7)
    
    # Equacao diferencial do nivel do reator
    xd[0] = (ux+par[0]+uy-(par[8]*(np.sqrt(x[0]-par[9]))))/par[2]; 
    
    # Equacao diferencial de Wa
    xd[1] = ((ux*(par[10]-x[1]))+(par[0]*(par[11]-x[1]))+(uy*(par[12]-x[1])))/(par[3])
    
    # Equacao diferencial de Wb
    xd[2] = ((ux*(par[13]-x[2]))+(par[0]*(par[14]-x[2]))+(uy*(par[15]-x[2])))/(par[3])
    
    # Equacao diferencial altura tanque acido
    xd[3] = -ux/par[4] 
    
    # Equacao diferencial altura tanque buffer
    xd[4] = -par[0]/par[5] 
    
    # Equacao diferencial altura tanque base 
    xd[5] = -uy/par[6]
    
    # Equacao diferencial altura tanque coletor (rejeito) - esta equacao eh usada quando houver controle de nivel
    xd[6] = par[1]/par[7] 
    
    
    xdot=xd
    
    return xdot


#%%

#Função simrk_pH
#Código orignal: simrk_pH.m

def simrk_pH(x0,u1,u2,h,t,par,Kas,Ts):

    # function [x,pH] = simrk_pH(x0,ux,uy,h,t,par,Kas,Ts)
    
    # x0 na entrada eh a condicao inicial para a iteracao em questao
    # u1 vazao de entrada acido forte (mL/seg)
    # u2 vazao de entrada base forte (mL/seg)
    # h eh o intervalo de integracao (nao deve ser alterado pelo usuario)
    # t eh o instante de tempo atual
    # par eh o vetor de parametros do modelo 
    # Kas parametros para calculo do pH
    # Ts tempo de amostragem
    # Ts/h deve ser um numero inteiro
    #
    # x na saida eh o vetor de variaveis de estado um tempo de amostragem
    # depois
    # pH eh o valor de pH um tempo de amostragem depois
    
    # simula a planta "continua" por um intervalo de amostragem Ts
    # durante o qual assume-se que as entradas permanecem contantes
    xc = np.zeros((len(x0),int(Ts/h)))

    pHc = np.zeros(int(Ts/h))
    xc[:,0] = x0
        
    coeff = [1, Kas[0] - xc[1,0], Kas[0]*Kas[1]-Kas[0]*xc[1,0]-Kas[2]-Kas[0]*xc[2,0], -(Kas[0]*Kas[2]+Kas[0]*Kas[1]*xc[1,0]+2*Kas[0]*Kas[1]*xc[2,0]), -Kas[0]*Kas[1]*Kas[2]]
    p=np.roots(coeff)
    
    ch = np.max(np.real(p))
    
    pHc[0] = -np.log10(ch)
    
    for i in range(1,int(Ts/h)):
        xc[:,i] = rkpH(xc[:,i-1],u1,u2,h,t,par)
        t = t+h
        
        #==========================================================================
        # Calculo do pH - Esse calculo nao faz parte da integracao
        #==========================================================================
        # Polinomio com a concentracao hidrogenio
        # [H+]^4+(ka1-x(i,2))*[H+]^3+(ka1*ka2-ka1*x(i,2)-kw-ka1*x(i,3))*....
        # [H+]^2-(ka1*kw+ka1*ka2*x(i,2)+2*ka1*ka2*x(i,3))*[H+]-ka1*ka2*kw=0;
        #==========================================================================
        
        coeff = [1, Kas[0] - xc[1,i], Kas[0]*Kas[1]-Kas[0]*xc[1,i]-Kas[2]-Kas[0]*xc[2,i], -(Kas[0]*Kas[2]+Kas[0]*Kas[1]*xc[1,i]+2*Kas[0]*Kas[1]*xc[2,i]), -Kas[0]*Kas[1]*Kas[2]]
        p=np.roots(coeff)
        
        ch = np.max(np.real(p))
        pHc[i] = -np.log10(ch)
        
    # O ultimo estado continuo eh "amostrado"
    x = xc[:,i]
    
    # O ultimo valore de pH eh "amostrado"
    pH = pHc[i]
    
    
    # Lembretes pessoais:
    #
    # indice do vetor de tempo continuo (t): kc
    # indice do vetor de tempo discreto (T): k
    # tempo discreto (nao usado no codigo): n
    # o mesmo instante de tempo eh alcancado fazendo
    # T(k)-h corresponde ao tempo discreto n ou seja t= nTs
    # t((k-1)*Ts/h + 1)-h corresponde ao mesmo instante da linha anterior,
    # ou seja, kc = (k-1)*Ts/h + 1.
    
        
    simrk_pH_dict = dict()
    simrk_pH_dict['x']=x
    simrk_pH_dict['pH']=pH
    simrk_pH_dict['xc']=xc
    simrk_pH_dict['pHc']=pHc

    return simrk_pH_dict
    

#%%
# a definicao de h deve ser "transparente ao usuaria", pois nao eh
# algo que ocorra em uma situacao experimental.
h = 10 # Intervalo de integracao em segundos

# Configuracoes de simulacao
t0 = h # Tempo inicial em segundos
tm = 50 # Tempo de simulacao em minutos 
tf = 60*tm # Tempo de simulacao em segundos
t=np.arange(t0,tf+h,h) # Vetor de tempo "continuo" em segundos


# tempo de amostragem. Este eh o valor que o usuario deve escolher para
# operar o sistema de controle em tempo discreto. NOTAR que Ts/h
# deve ser um numero inteiro
Ts = 20

# indexação python : <Elemento inicial>:<Elemento final>:<De quantos em quantos>
T = t[::int(Ts/h)] # Vetor de tempo "discreto" em segundos


#==========================================================================
# Definicao do ponto de operacao, em termos das entradas. Variacoes,
# sejam em malha aberta ou fechada ocorrem "em torno" de Q1 e Q3.
#==========================================================================
# 
# vazao de entrada acido forte (ml/seg)
Q1 = 3*np.ones(len(T))

# vazao de entrada base forte (ml/seg)
Q3 = 2*np.ones(len(T))
u = np.zeros(len(T))

##
#==========================================================================
# Integracao Numerica inicial, somente para atingir o poonto de operacao
# qualquer teste em malha aberta ou fechada deve ser realizado apos esse
# periodo que foi escolhido igual a 10 minutos = 600 segundos
#==========================================================================
# Vetor de estado e pH continuo
xc = np.zeros((len(x0),len(t)))
pHc = np.zeros(len(t))
# Vetor de estado e pH discreto
x = np.zeros((len(x0),len(T)))
y1 = np.zeros(len(T))
y2 = np.copy(y1)
y3 = np.copy(y1)
y4 = np.copy(y1)
y5 = np.copy(y1)
pH = np.zeros(len(T))

x[:,0] = x0

# regime transitorio
for k in range(1,int(600/Ts+1)):
    
    kc = int(k*Ts/h)
    
    simrk_pH_dict=simrk_pH(x[:,k-1],Q1[k],Q3[k],h,t[kc],par,Kas,Ts)
    
    x[:,k]=simrk_pH_dict['x']
    pH[k]=simrk_pH_dict['pH'] 
    xc[:,kc:kc+int(Ts/h)]=simrk_pH_dict['xc']
    pHc[kc:kc+int(Ts/h)]=simrk_pH_dict['pHc']
    


for k in range(5,int(600/Ts+1)):
    y1[k] = 0.7770*y1[k-1] + 0.6177*u[k-1] 
    y2[k] = 0.8521*y2[k-1] + 0.6234*u[k-1] -0.4327*u[k-3] + 0.2190*u[k-4]
    y3[k] = 1.4042*y3[k-1] -0.4983*y3[k-2] +0.3569*u[k-1] -0.5219*u[k-4] +0.4256*u[k-5]
    y4[k] = 0.8888*y4[k-1] +0.5221*u[k-1] -0.4247*u[k-4] +0.2106*u[k-5]
    y5[k] = 1.3316*y5[k-1] -0.4182*y5[k-2] +0.2399*u[k-1]



# valor de kc inicial para o periodo de testes
ini = k
##########################################################################
# Testes em malha aberta ou fechada a partir daqui

# ponto de operacao
u1 = Q1
u2 = Q3
# delta u2, variacoes de u2 em torno do ponto de operacao
du2 = np.zeros(len(Q3))

# mudanca em degrau na vazao de base forte a partir de 10min
du2[ini:] = 0.5*np.ones(len(T)-ini)
u2 = Q3 + du2

# regime de teste
for k in range(ini,len(T)):
    
    kc = int(k*Ts/h)
    
    simrk_pH_dict = simrk_pH(x[:,k-1],u1[k],u2[k],h,t[kc],par,Kas,Ts)
    
    x[:,k]=simrk_pH_dict['x']
    pH[k]=simrk_pH_dict['pH']
    xc[:,kc:kc+int(Ts/h)]=simrk_pH_dict['xc']
    pHc[kc:kc+int(Ts/h)]=simrk_pH_dict['pHc']
  

u = du2
for k in range(ini,len(T)):
    y1[k] = 0.7770*y1[k-1] + 0.6177*u[k-1] 
    y2[k] = 0.8521*y2[k-1] + 0.6234*u[k-1] -0.4327*u[k-3] + 0.2190*u[k-4]
    y3[k] = 1.4042*y3[k-1] -0.4983*y3[k-2] +0.3569*u[k-1] -0.5219*u[k-4] +0.4256*u[k-5]
    y4[k] = 0.8888*y4[k-1] +0.5221*u[k-1] -0.4247*u[k-4] +0.2106*u[k-5]
    y5[k] = 1.3316*y5[k-1] -0.4182*y5[k-2] +0.2399*u[k-1]

#%% 


pplot = 1
if pplot == 1:
    
    #gera_graficos 
    plt.figure(1, figsize=(10, 8.107))
    plt.plot(t/60,pHc,color='black',linewidth=1.5)
    plt.scatter((T+Ts)/60,pH,color='blue',marker='*',s=100)      
    plt.step((T+Ts)/60,u2,color='red',linewidth=1.5)
    plt.xlim(0, 50)
    plt.ylim(0, 7)
    plt.xlabel('t(min)', fontsize=17.5)
    plt.ylabel('pH', fontsize=17.5, labelpad=15)
    plt.xticks(fontsize=17.5) # Tamanho da fonte dos números do eixo x
    plt.yticks(fontsize=17.5)
    #plt.savefig('fig_1.eps', format='eps')
    plt.show()

  
plt.figure(1, figsize=(10, 8.107))
plt.plot(t/60,pHc,color='blue',linewidth=1.5)
plt.scatter(T/60,y1+5,color='black',marker='+',s=100)
plt.scatter(T/60,y2+5,color='red',marker='+',s=100)  
plt.scatter(T/60,y3+5,color='green',marker='+',s=100)  
plt.scatter(T/60,y4+5,color='blue',marker='+',s=100)  
plt.scatter(T/60,y5+5,color='cyan',marker='+',s=100)  
plt.xlim(9, T[-1]/60)
plt.ylim(4.9, 6.5)
plt.xlabel('t(min)', fontsize=17.5)
plt.ylabel('pH', fontsize=17.5, labelpad=15)
plt.xticks(fontsize=17.5)
plt.yticks(fontsize=17.5)
#plt.savefig('fig_2.eps', format='eps')
plt.show()

