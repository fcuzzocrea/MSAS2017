function [ dx ] = odefun(t,x,p)

% parameters
C = p(1);               % Capacitance of the piezo stack
Tem = p(2);             % To take in account the piezo effect
R = p(3);               % Resistance of the electrical circuit
m = p(4);               % Mass of the Relase Tip
b = p(5);               % Damping coefficient 
k = p(6);               % Stiffness coefficient

% LTI state space model definition
A = [-1/(C*R),     Tem/(C*R)     , 0;...
            0,         0         , 1;...
    Tem/(C*m), (-k/m)-Tem^2/(C*m), -b/m];

B = [1/R;              0;             0];

% Time Shift (lag time)
T = 5e-6;

% Input modeled as a step
u = 60*(sign(t-T)+1);     % This returns -1 if t < 5e-6 otherwise +1

% State space model
dx = A*x+B*u;               

end