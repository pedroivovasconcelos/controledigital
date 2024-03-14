function xdot = dvpH(x,ux,uy,t,par)
% function xdot = dvpH(x,ux,uy,t,par)
% xdot eh o vetor de derivadas (campo vetorial) usado pelo algoritmo de
% integracao 
% ux e uy sao os valores das entradas mantidas constantes durante o
% intervalor de amostragem
% t eh o instante de tempo atual. Usado apenas como registro do instante a
% que se refere a avaliacao do campo vetorial. Normalmente nao eh usado.
% par eh o vetor de parametros do modelo

%-------------------------------------------------------------------------%
% Equacoes diferenciais                                                   %
%-------------------------------------------------------------------------%
% x(1) Nivel do Reator
% x(2) Concentracao de acido
% x(3) Concentracao de Base
% x(4) Nivel Tanque de acido
% x(5) Nivel Tanque de Tampao
% x(6) Nivel Tanque de Base
% x(7) Nivel Tanque do Coletor
% ux vazao de entrada acido forte (mL/seg)
% xy vazao de entrada base forte (mL/seg)

% Equacao diferencial do nivel do reator
xd(1) = (ux+par(1)+uy-(par(9)*(sqrt(x(1)-par(10)))))/par(3); 

% Equacao diferencial de Wa
xd(2) = ((ux*(par(11)-x(2)))+(par(1)*(par(12)-x(2)))+(uy*(par(13)-x(2))))/(par(4)); 

% Equacao diferencial de Wb
xd(3) = ((ux*(par(14)-x(3)))+(par(1)*(par(15)-x(3)))+(uy*(par(16)-x(3))))/(par(4));

%Equacao diferencial altura tanque acido
xd(4) = -ux/par(5); 

%Equacao diferencial altura tanque buffer
xd(5) = -par(1)/par(6); 

%Equacao diferencial altura tanque base 
xd(6) = -uy/par(7); 

% Equacao diferencial altura tanque coletor (rejeito) - esta equacao eh usada quando houver controle de nivel
xd(7) = par(2)/par(8); 

xdot=xd';
end
