% linearizacao do modelo de pH

% matriz dinamica
A = zeros(7,7);
% matriz de entrada
B = zeros(7,2);

% valoes no ponto de operacao
x1 = 5.4064;
x2 = -3.5872e-05;
x3 = 1.1120e-03;
ux = 3;
uy = 2;

A(1,1) = -par(9)/(2*par(3)*sqrt(x1-par(10)));
A(2,2) = (-ux-par(1)-uy)/par(4);
A(3,3) = (-ux-par(1)-uy)/par(4);

B(1,1) = 1/par(3); 
B(1,2) = 1/par(3);
B(2,1) = (par(11)-x2)/par(4);
B(2,2) = (par(13)-x2)/par(4);
B(3,1) = (par(14)-x3)/par(4);
B(3,2) = (par(16)-x3)/par(4);
B(4,1) = -1/par(5);
B(6,2) = -1/par(7); 

C = [1 0 0 0 0 0 0];
D = [0 0];

sysc = ss(A,B,C,D); 
%[Ad,Bd,C,D] = c2d(A,B,C,D,40);
sysd = c2d(sysc,40,'zoh');
Ad = sysd.a;
Bd = sysd.b;

