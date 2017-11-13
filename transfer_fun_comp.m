clearvars
close all

syms s Tem c k b R m
A = [-1/(c*R), Tem/(c*R) 0;...
    0, 0, 1;...
    Tem/(c*m), (-k/m)-Tem^2/(c*m), -b/m];
B = [1/R; 0; 0];
C = [0, 1, 0];
F = C*inv((s*eye(3)-A))*B;