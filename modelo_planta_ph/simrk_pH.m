function [x,pH,xc,pHc] = simrk_pH(x0,u1,u2,h,t,par,Kas,Ts)
% function [x,pH] = simrk_pH(x0,ux,uy,h,t,par,Kas,Ts)

% x0 na entrada eh a condicao inicial para a iteracao em questao
% u1 vazao de entrada acido forte (mL/seg)
% u2 vazao de entrada base forte (mL/seg)
% h eh o intervalo de integracao (nao deve ser alterado pelo usuario)
% t eh o instante de tempo atual
% par eh o vetor de parametros do modelo 
% Kas parametros para calculo do pH
% Ts tempo de amostragem
% Ts/h deve ser um numero inteiro
%
% x na saida eh o vetor de variaveis de estado um tempo de amostragem
% depois
% pH eh o valor de pH um tempo de amostragem depois

% simula a planta "continua" por um intervalo de amostragem Ts
% durante o qual assume-se que as entradas permanecem contantes
xc = zeros(length(x0),Ts/h);
pHc = zeros(1,Ts/h);
xc(:,1) = x0;
p = roots([1 (Kas(1) - xc(2,1)) (Kas(1)*Kas(2)-Kas(1)*xc(2,1)-Kas(3)-Kas(1)*xc(3,1)) -(Kas(1)*Kas(3)+Kas(1)*Kas(2)*xc(2,1)+2*Kas(1)*Kas(2)*xc(3,1)) -Kas(1)*Kas(2)*Kas(3)]);
ch = max(real(p));
pHc(1) = -log10(ch);

for i=2:Ts/h
    xc(:,i) = rkpH(xc(:,i-1),u1,u2,h,t,par);
    t = t+h;
    
    %==========================================================================
    % Calculo do pH - Esse calculo nao faz parte da integracao
    %==========================================================================
    % Polinomio com a concentracao hidrogenio
    % [H+]^4+(ka1-x(i,2))*[H+]^3+(ka1*ka2-ka1*x(i,2)-kw-ka1*x(i,3))*....
    % [H+]^2-(ka1*kw+ka1*ka2*x(i,2)+2*ka1*ka2*x(i,3))*[H+]-ka1*ka2*kw=0;
    %==========================================================================

    p = roots([1 (Kas(1) - xc(2,i)) (Kas(1)*Kas(2)-Kas(1)*xc(2,i)-Kas(3)-Kas(1)*xc(3,i)) -(Kas(1)*Kas(3)+Kas(1)*Kas(2)*xc(2,i)+2*Kas(1)*Kas(2)*xc(3,i)) -Kas(1)*Kas(2)*Kas(3)]);
    ch = max(real(p));
    pHc(i) = -log10(ch); 
end

% O ultimo estado continuo eh "amostrado"
x = xc(:,i);

% O ultimo valore de pH eh "amostrado"
pH = pHc(i);




% Lembretes pessoais:
%
% indice do vetor de tempo continuo (t): kc
% indice do vetor de tempo discreto (T): k
% tempo discreto (nao usado no codigo): n
% o mesmo instante de tempo eh alcancado fazendo
% T(k)-h corresponde ao tempo discreto n ou seja t= nTs
% t((k-1)*Ts/h + 1)-h corresponde ao mesmo instante da linha anterior,
% ou seja, kc = (k-1)*Ts/h + 1.


