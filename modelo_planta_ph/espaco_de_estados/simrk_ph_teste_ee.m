%**************************************************************************
%* Simulacao de planta de neutralizacao de pH  - em espaco de estados     *
%* definida em termos de equacoes diferenciais e controle em termos de    *
%* equacao de diferencas. Portanto a malha eh hibrida, como em um caso    *
%* real. 
%
%* LAA 05/12/22
%**************************************************************************
clearvars; close all; clc

% carrega os parametros do modelo continuo
setup_pH;

% carrega matrizes em espaco de estado discreto linearizadas em torno
% do ponto de operacao u1 = 3 e u2 = 2;
pH_linear;

% a definicao de h deve ser "transparente ao usuaria", pois nao eh
% algo que ocorra em uma situacao experimental.
h = 10; % Intervalo de integracao em segundos

% Configuracoes de simulacao
t0 = h; % Tempo inicial em segundos
tm = 50; % Tempo de simulacao em minutos 
tf = 60*tm; % Tempo de simulacao em segundos
t = t0:h:tf; % Vetor de tempo "continuo" em segundos

% tempo de amostragem. Este eh o valor que o usuario deve escolher para
% operar o sistema de controle em tempo discreto. NOTAR que Ts/h
% deve ser um numero inteiro
Ts = 40;
T = t(1:Ts/h:end); % Vetor de tempo "discreto" em segundos


%==========================================================================
% Definicao do ponto de operacao, em termos das entradas. Variacoes,
% sejam em malha aberta ou fechada ocorrem "em torno" de Q1 e Q3.
%==========================================================================
% 
% vazao de entrada acido forte (ml/seg)
Q1 = 3*ones(1,length(T)); 

% vazao de entrada base forte (ml/seg)
Q3 = 2*ones(1,length(T)); 
u = zeros(1,length(T)); 

%%
%==========================================================================
% Integracao Numerica inicial, somente para atingir o poonto de operacao
% qualquer teste em malha aberta ou fechada deve ser realizado apos esse
% periodo que foi escolhido igual a 10 minutos = 600 segundos
%==========================================================================
% Vetor de estado e pH continuo
xc = zeros(length(x0),length(t)); 
pHc = zeros(1,length(t));

% Vetor de estado e pH discreto
x = zeros(length(x0),length(T));
pH = zeros(1,length(T));
x(:,1) = x0;

% Vetor de estado e pH discreto LINEAR eh igual 
x0l = [5.4064e+00
  -3.5872e-05
   1.1120e-03
   4.4495e+01
   3.6338e+01
   4.4814e+01
            0];
x0l = zeros(7,1);
xl = zeros(length(x0),length(T));
pHl = zeros(1,length(T));
xl(:,1) = x0l;

% regime transitorio
for k = 2:600/Ts+1
    kc = (k-1)*Ts/h + 1;
    [x(:,k),pH(k),xc(:,kc:kc+Ts/h-1),pHc(kc:kc+Ts/h-1)] = simrk_pH(x(:,k-1),Q1(k),Q3(k),h,t(kc),par,Kas,Ts);
   
    % modelo linear  
    xl(1,k) = Ad(1,1)*xl(1,k-1);
    xl(2,k) = Ad(2,2)*xl(2,k-1);
    xl(3,k) = Ad(3,3)*xl(3,k-1);  
    xl(4,k) = Ad(4,4)*xl(4,k-1);
    xl(5,k) = Ad(5,5)*xl(5,k-1);
    xl(6,k) = Ad(6,6)*xl(6,k-1);
    xl(7,k) = Ad(7,7)*xl(7,k-1);
    
    p = roots([1 (Kas(1) - (xl(2,k)+x2)) (Kas(1)*Kas(2)-Kas(1)*(xl(2,k)+x2)-Kas(3)-Kas(1)*(xl(3,k)+x3)) -(Kas(1)*Kas(3)+Kas(1)*Kas(2)*(xl(2,k)+x2)+2*Kas(1)*Kas(2)*(xl(3,k)+x3)) -Kas(1)*Kas(2)*Kas(3)]);
    ch = max(real(p));
    pHl(k) = -log10(ch);
end



% valor de kc inicial para o periodo de testes
ini = k;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Testes em malha aberta ou fechada a partir daqui

% ponto de operacao
u1 = Q1;
u2 = Q3;
% delta u1, variacoes de u1 em torno do ponto de operacao
du1 = zeros(size(Q3));
% delta u2, variacoes de u2 em torno do ponto de operacao
du2 = zeros(size(Q3));

% mudanca em degrau na vazao de base forte a partir de 10min
du2(ini:end) = 0.25*ones(1,length(T)-ini+1); 
u2 = Q3 + du2;

% regime de teste
for k=ini:length(T)
    kc = (k-1)*Ts/h + 1;
    [x(:,k),pH(k),xc(:,kc:kc+Ts/h-1),pHc(kc:kc+Ts/h-1)] = simrk_pH(x(:,k-1),u1(k),u2(k),h,t(kc),par,Kas,Ts);

  
    % modelo linear  
    xl(1,k) = Ad(1,1)*xl(1,k-1) +Bd(1,1)*du1(k) +Bd(1,2)*du2(k);
    xl(2,k) = Ad(2,2)*xl(2,k-1) +Bd(2,1)*du1(k) +Bd(2,2)*du2(k);
    xl(3,k) = Ad(3,3)*xl(3,k-1) +Bd(3,1)*du1(k) +Bd(3,2)*du2(k);  
    xl(4,k) = Ad(4,4)*xl(4,k-1) +Bd(4,1)*du1(k);
    xl(5,k) = Ad(5,5)*xl(5,k-1);
    xl(6,k) = Ad(6,6)*xl(6,k-1) +Bd(6,2)*du2(k);
    xl(7,k) = Ad(7,7)*xl(7,k-1);
         
    p = roots([1 (Kas(1) - (xl(2,k)+x2)) (Kas(1)*Kas(2)-Kas(1)*(xl(2,k)+x2)-Kas(3)-Kas(1)*(xl(3,k)+x3)) -(Kas(1)*Kas(3)+Kas(1)*Kas(2)*(xl(2,k)+x2)+2*Kas(1)*Kas(2)*(xl(3,k)+x3)) -Kas(1)*Kas(2)*Kas(3)]); 
    ch = max(real(p));
    pHl(k) = -log10(ch);
end



pplot = 1;
if pplot == 1

   %gera_graficos 
   figure(1)
   yyaxis left
   plot((T+Ts)/60,pH,'k',(T+Ts)/60,pHl,'b*');
   set(gca,'FontSize',16)
   ylabel('pH')
   yyaxis right
   plot((T+Ts)/60,u2,'r');
   set(gca,'FontSize',16)
   xlabel('t (min)')
   ylabel('vazao de base (mL/seg)')
   axis([0 (t(end)+h)/60 0 12])
end




