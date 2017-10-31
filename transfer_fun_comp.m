clearvars
close all

syms s Tem c kp k b bp R m mp
A = [-8*bp/mp, 4*bp/mp, -8*kp/mp, 4*kp/mp, 0;...
    2*bp/(m+mp/4), -(2*bp+b)/(m+mp/4), 2*kp/(m+mp/4), -(2*kp+k+Tem^2/c)/(m+mp/4), Tem/c/(m+mp/4);...
    1, 0, 0, 0, 0;...
    0, 1, 0, 0, 0;...
    0, 0, 0, Tem/(R*c), -1/(R*c)];
B = [0; 0; 0; 0; 1/R];
C = [0, 0, 0, 1, 0];
F = C*inv((s*eye(5)-A))*B;