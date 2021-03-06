%% function to define parametric system for identification
function [A,B,C,D] = par_system(par,T)

m = par(1); c = par(2); Tem = par(3); kp = par(4);
bp = par(5); mp = par(6); k = par(7); b = par(8); R = par(9);
A = [-8*bp/mp, 4*bp/mp, -8*kp/mp, 4*kp/mp, 0;...
    2*bp/(m+mp/4), -(2*bp+b)/(m+mp/4), 2*kp/(m+mp/4), -(2*kp+k+Tem^2/c)/(m+mp/4), Tem/c/(m+mp/4);...
    1, 0, 0, 0, 0;...
    0, 1, 0, 0, 0;...
    0, 0, 0, Tem/(R*c), -1/(R*c)];
 B = [0; 0; 0; 0; 1/R];
 C = [0, 0, 0, 1, 0];
 D = 0;