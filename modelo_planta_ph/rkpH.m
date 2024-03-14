function [x,xc] = rkpH(x0,ux,uy,h,t,par)
% function x = rkpH(x0,ux,uy,h,t,par%
% x0 eh a condicao inicial para a iteracao em questao
% ux vazao de entrada acido forte (mL/seg)
% xy vazao de entrada base forte (mL/seg)
% h eh o intervalo de integracao (nao deve ser alterado pelo usuario)
% t eh o instante de tempo atual
% par eh o vetor de parametros do modelo 
%
% x eh o vetor de estado em tempo discreto
% xc eh o vetor de estado em tempo continuo

% 1st evaluation at 
xd = dvpH(x0,ux,uy,t,par);
savex0 = x0;
phi = xd;
x0 = savex0+0.5*h*xd;

% 2nd evaluation 
xd = dvpH(x0,ux,uy,t+0.5*h,par);
phi = phi+2*xd;
x0 = savex0+0.5*h*xd;

% 3rd evaluation
xd = dvpH(x0,ux,uy,t+0.5*h,par);
phi = phi+2*xd;
x0 = savex0+h*xd;

% 4th evaluation
xd = dvpH(x0,ux,uy,t+h,par);
x = savex0+(phi+xd)*h/6;
end