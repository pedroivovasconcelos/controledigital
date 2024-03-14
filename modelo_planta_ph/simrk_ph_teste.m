%**************************************************************************
%* Simulacao de planta de neutralizacao de pH  - com entrada PRBS -       *
%* definida em termos de equacoes diferenciais e controle em termos de    *
%* equacao de diferencas. Portanto a malha eh hibrida, como em um caso    *
%* real. 
%
%* LAA 04/08/22
%**************************************************************************
clearvars; close all; clc

% carrega os parametros do modelo]
setup_pH;

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
y1 = zeros(1,length(T));
y2 = y1;
y3 = y1;
y4 = y1;
y5 = y1;
pH = zeros(1,length(T));

x(:,1) = x0;


% regime transitorio
for k = 2:600/Ts+1
    kc = (k-1)*Ts/h + 1;
    [x(:,k),pH(k),xc(:,kc:kc+Ts/h-1),pHc(kc:kc+Ts/h-1)] = simrk_pH(x(:,k-1),Q1(k),Q3(k),h,t(kc),par,Kas,Ts);

end


for k = 6:600/Ts+1
    y1(k) = 0.7770*y1(k-1) + 0.6177*u(k-1); 
    y2(k) = 0.8521*y2(k-1) + 0.6234*u(k-1) -0.4327*u(k-3) + 0.2190*u(k-4);
    y3(k) = 1.4042*y3(k-1) -0.4983*y3(k-2) +0.3569*u(k-1) -0.5219*u(k-4) +0.4256*u(k-5);
    y4(k) = 0.8888*y4(k-1) +0.5221*u(k-1) -0.4247*u(k-4) +0.2106*u(k-5);
    y5(k) = 1.3316*y5(k-1) -0.4182*y5(k-2) +0.2399*u(k-1);
end

% valor de kc inicial para o periodo de testes
ini = k;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Testes em malha aberta ou fechada a partir daqui

% ponto de operacao
u1 = Q1;
u2 = Q3;
% delta u2, variacoes de u2 em torno do ponto de operacao
du2 = zeros(size(Q3));

% mudanca em degrau na vazao de base forte a partir de 10min
du2(ini:end) = 0.5*ones(1,length(T)-ini+1); 
u2 = Q3 + du2;

% regime de teste
for k=ini:length(T)
    kc = (k-1)*Ts/h + 1;
    [x(:,k),pH(k),xc(:,kc:kc+Ts/h-1),pHc(kc:kc+Ts/h-1)] = simrk_pH(x(:,k-1),u1(k),u2(k),h,t(kc),par,Kas,Ts);
end

u = du2;
for k = ini:length(T)
   y1(k) = 0.7770*y1(k-1) + 0.6177*u(k-1); 
    y2(k) = 0.8521*y2(k-1) + 0.6234*u(k-1) -0.4327*u(k-3) + 0.2190*u(k-4);
    y3(k) = 1.4042*y3(k-1) -0.4983*y3(k-2) +0.3569*u(k-1) -0.5219*u(k-4) +0.4256*u(k-5);
    y4(k) = 0.8888*y4(k-1) +0.5221*u(k-1) -0.4247*u(k-4) +0.2106*u(k-5);
    y5(k) = 1.3316*y5(k-1) -0.4182*y5(k-2) +0.2399*u(k-1);
end


pplot = 1;
if pplot == 1

   %gera_graficos 
   figure(1)
   yyaxis left
   plot(t/60,pHc,'k',(T+Ts)/60,pH,'b*');
   set(gca,'FontSize',16)
   ylabel('pH')
   yyaxis right
   plot((T+Ts)/60,u2,'r');
   set(gca,'FontSize',16)
   xlabel('t (min)')
   ylabel('vazao de base (mL/seg)')
   axis([0 (t(end)+h)/60 0 12])
end

%%
figure(2)
plot((T+Ts)/60,pH,'b',T/60,y1+5,'k+',T/60,y2+5,'r+',T/60,y3+5,'g+',T/60,y4+5,'b+',T/60,y5+5,'c+');
axis([9 T(end)/60 4.9 6.5])
set(gca,'FontSize',16)
ylabel('pH')
xlabel('t (min)')


