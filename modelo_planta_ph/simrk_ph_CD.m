%**************************************************************************
%* Simulacao de malha de controle digital: planta de neutralizacao de pH  *
%* definida em termos de equacoes diferenciais e controle em termos de    *
%* equacao de diferencas. Portanto a malha eh hibrida, como em um caso    *
%* real. O usuario idealmente deve mexer somente na parte correspondente  *
%* ao controle digital.
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
tr = 60*15; % instante em que ocorre a mudanca na referencia
td = 60*10; % quanto tempo depois ocorre o disturbio. Tempo do disturbio eh tr+td
t = t0:h:tf; % Vetor de tempo "continuo" em segundos

% tempo de amostragem. Este eh o valor que o usuario deve escolher para
% operar o sistema de controle em tempo discreto. NOTAR que Ts/h
% deve ser um numero inteiro
Ts = 40;
T = t(1:Ts/h:end); % Vetor de tempo "discreto" em segundos

% constantes do controlador
K = 2;

%
%==========================================================================
% Definicao do ponto de operacao, em termos das entradas. Variacoes,
% sejam em malha aberta ou fechada ocorrem "em torno" de Q1 e Q3.
%==========================================================================
% 
% vazao de entrada acido forte (ml/seg)
Q1 = 3*ones(1,length(T)); 

% disturbio em Q1
d1 = zeros(1,length(T)); 

% vazao de entrada base forte (ml/seg)
Q3 = 2*ones(1,length(T)); 

% referencia (setpoint)
r = 5*ones(1,length(T));

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
% erro
e = zeros(1,length(T));
% acao de controle
m = zeros(1,length(T));

x(:,1) = x0;


% regime transitorio
for k = 2:floor(tr/Ts)+1
    kc = (k-1)*Ts/h + 1;
    [x(:,k),pH(k),xc(:,kc:kc+Ts/h-1),pHc(kc:kc+Ts/h-1)] = simrk_pH(x(:,k-1),Q1(k),Q3(k),h,t(kc),par,Kas,Ts);

end



% valor de kc inicial para o periodo de testes
ini = k;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Testes em malha aberta ou fechada a partir daqui

% ponto de operacao
u1 = Q1;
u2 = Q3;
% m sao variacoes de u2 em torno do ponto de operacao que serao
% determinadas pelo controle
m = zeros(size(Q3));

% referencia (setpoint) 
r(ini:end) = 5.5*ones(1,length(T)-ini+1); 

% disturbio em degrau 
d1(ini+floor(td/Ts):end) = 0.3*ones(1,length(T)-ini-floor(td/Ts)+1); 

% regime de teste
for k=ini:length(T)
    % erro
    e(k-1) = r(k-1)-pH(k-1);
      
    % controlador proporcional
    m(k-1) = K*e(k-1);
    u2(k-1) = Q3(k-1)+m(k-1);
    
    % disturbio
    u1(k-1) = Q1(k-1) - d1(k-1);
    
    % planta
    kc = (k-1)*Ts/h + 1;
    [x(:,k),pH(k),xc(:,kc:kc+Ts/h-1),pHc(kc:kc+Ts/h-1)] = simrk_pH(x(:,k-1),u1(k-1),u2(k-1),h,t(kc),par,Kas,Ts);
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
   stairs((T+Ts)/60,u2,'r');
   set(gca,'FontSize',16)
   xlabel('t (min)')
   ylabel('vazao de base (mL/seg)')
   axis([0 (t(end)+h)/60 0 12])
end


