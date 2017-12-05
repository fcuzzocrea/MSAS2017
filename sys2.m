function [A,B,C,D,K,x0] = sys2(par,t)
% State space representation of the model

% Parameters
C = par(1);
Tem = par(2);
R = par(3);
m = par(4);
b = par(5);
k = par(6);
V0 = par(7);

% Model
A = [-1/(C*R), Tem/(C*R) 0;...
    0, 0, 1;...
    Tem/(C*m), (-k/m)-Tem^2/(C*m), -b/m];

B = [1/R; 0; 0];

C = [0, 1, 0];

D = 0;

K = zeros(3,1);

x0 = [0; 0; V0];