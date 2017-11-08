%% function to define parametric system for identification
function [A,B,C,D,K,X0] = par_system(par,t)

m = par(1);
mp = par(2);
k = par(3);
kp = par(4);
b = par(5); 
bp = par(6);
c = par(7);
Tem = par(8);
R = par(9);
V0 = par(10);

A = [-8*bp/mp, 4*bp/mp, -8*kp/mp, 4*kp/mp, 0;...
    2*bp/(m+mp/4), -(2*bp+b)/(m+mp/4), 2*kp/(m+mp/4), -(2*kp+k+Tem^2/c)/(m+mp/4), Tem/c/(m+mp/4);...
    1, 0, 0, 0, 0;...
    0, 1, 0, 0, 0;...
    0, 0, 0, Tem/(R*c), -1/(R*c)];
 B = [0; 0; 0; 0; 1/R];
 C = [0, 0, 0, 1, 0];
 D = 0;
 K = zeros(5,1);
 X0 = [0; V0/2; V0; 0; 0];