function dx = odefun(t,x,p)

% parameters
c = p(1);
Tem = p(2);
R = p(3);
m = p(4);
b = p(5);
k = p(6);

% LTI model definition
A = [-1/(c*R), Tem/(c*R) 0;...
    0, 0, 1;...
    Tem/(c*m), (-k/m)-Tem^2/(c*m), -b/m];
B = [1/R; 0; 0];

u = 60*(sign(t-5e-6)+1);     % input
dx = A*x+B*u;                % ODE