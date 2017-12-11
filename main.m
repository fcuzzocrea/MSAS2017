%{
   MODELING AND SIMULATION OF AEROSPACE SYSTEMS, Fall 2017
   (C) Collogrosso Alfonso, Cuzzocrea Francescodario, Mastrantuono Andrea
   Politecnico di Milano
   WEB : https://github.com/fcuzzocrea/MSAS2017.git
%}

clearvars
close all
clear 
clc

% Preliminary Datas
M = 1.9;                % TM mass in [kg]
I = 6.9e-4;             % Inertia [kg*m^2] w.r.t. any of the 3 axes

%% IDENTIFICATION

% Define system to identify
par =  [4.8e-7, 1.5, 400, 1e-3, 150, 1.5e7, 0]; % par = [Ca Tem R m b k V0]
aux = {};                                       % Estimation optional arguments
T = 0;                                          % Parameter to specify model as continous-time
sys = idgrey('sys2',par,'c',aux,T);

% Set constraints on parameters
sys.structure.Parameters.Minimum(1) = 0;        % min value for C
sys.structure.Parameters.Minimum(2) = 0;        % min value for Tem
sys.structure.Parameters.Minimum(3) = 100;      % min value for R
sys.structure.Parameters.Minimum(4) = 1e-4;     % min value for m
sys.structure.Parameters.Minimum(5) = 100;      % min value for b
sys.structure.Parameters.Minimum(6) = 0;        % min value for k
sys.structure.Parameters.Minimum(7) = -0.05;    % min value for V0
sys.structure.Parameters.Maximum(1) = 1e-6;     % max value for C
sys.structure.Parameters.Maximum(3) = 800;      % max value for R
sys.structure.Parameters.Maximum(4) = 0.014;    % max value for m
sys.structure.Parameters.Maximum(7) = 0.05;     % max value for V0


% Load and resample datas to fit in identification
load('responce_ref2.mat');                              % Displacements in microm, time in micros
ti = t_exp(1);
tf = t_exp(end);
ts = 1e-6;                                              % Sampling time (f=1MHz)
t_sim = ti:ts:tf;
y_exp = interp1(t_exp,d_exp,t_sim,'linear');
n = 10;                                                       % Time lag (in time samples)
input_voltage = 120;                                          % Volt
input = [zeros(1,n),input_voltage*ones(1,length(t_sim)-n)];   % Input signal
data = iddata(y_exp',input',ts);                              % Experimental response

% Identification
opt = greyestOptions('InitialState','model','Display','on','SearchMethod','lm','DisturbanceModel','none');
id_sys = pem(data,sys,opt);
p = getpvec(id_sys);                                    % Parameters of the identified system

% Compare results to validate
[A,B,C,~,~,x0] = sys2(pvec,0);                          % Evaluate the system with estimated datas
figure(1)
opt = compareOptions('InitialCondition',x0);
compare(data,id_sys,Inf,opt)
title('Experimental response vs. Identified model response')
xlabel('Time')
ylabel('Elongation [m]')
grid on


%% Compute adhesion forces fitting
% Load experimental datas

load('adhesion')            %load adhesion experimental data:
                            %ydata = adhesion force
                            %xdata = elongation
                            
% Data Fitting
X_adh     = zeros(3,3);     
sigma_adh = zeros(3,3);     

X     = zeros(5,3);         
sigma = zeros(5,1);         

for n = 1:3                     
    for m = 1:5                 
        [X(m,:), sigma(m)] = lsqcurvefit(@(x,xdata) x(1)*xdata...
            .*exp(-x(2)*(xdata.^x(3))),ones(3,1),max(set(n).exp(m).el,zeros(size(set(n)...
            .exp(m).el))),max(set(n).exp(m).F,zeros(size(set(n).exp(m).F))));
    end
    X_adh(n,:) = sum(X)/5;
    sigma_adh(n,:) = sqrt(sum((X-X_adh(n,:)).^2)/4);
end


%% Simulation with parameters covariance matrix estimation

% Residuals computation
[t_i,x_i] = ode113(@(t,x) odefun(t,x,p),t_sim,x0);
y_i = C*x_i';                       % Modeled tip responce
res = y_i - y_exp;                  % Model error w.r.t. experiment
s_y = res*res'/length(res);         % Variance of the residual 
y_rms = sqrt(s_y);                  % R.M.S. error
figure(2)
plot(t_i,res)
grid on
hold on
title('Residuals : $$ \varepsilon = y - \hat{y}$$','interpreter','latex')
xlabel('Time [s]')
ylabel('Residual value [m]')
xL = xlim;
line(xL, [y_rms y_rms]);  %x-axis
str = '$$ \varepsilon_{RMS} = 1.3967e-07 $$';
text(3*1e-4,y_rms+0.2*y_rms,str,'Interpreter','latex')


% Stability Analysis
figure(3)
l = eig(A);
plot(real(l),imag(l),'*')
stiff = max(abs(real(l)))/min(abs(real(l)));
axis square
grid on
xlabel('Re')
ylabel('Im')

% Plot Runge-Kutta stability regions
figure(4)
load('chebfun');
z = exp(1i*t);
w = z-1;
for i = 1:3
  w = w-(1+w+.5*w.^2-z.^2)./(1+w);
end
plot(w)                                % order 2
hold on
for i = 1:4
  w = w-(1+w+.5*w.^2+w.^3/6-z.^3)./(1+w+w.^2/2);
end
plot(w)                                % order 3
t_step = 1.0000e-06;
plot(real(l)*(t_step),imag(l)*(t_step),'*')
axis([-5 2 -3.5 3.5])
grid on

% Jacobian estimation:
% This estimation is performed by simulating the responce varing one
% single parameter each time and comparing it with the nominal parameters
% system responce

J = zeros(length(y_i),length(p));
dp = 0.01;                                % ? Perchè ? 

for i = 1:length(p)
    p1 = p;
    p1(i) = p(i)*(1+dp);
    [~,~,C,~,~,x0] = sys2(p1,0);
    [T,X] = ode113(@(t,x) odefun(t,x,p1),t_sim,x0);
    y = C*X';
    J(:,i) = (y-y_i)./(dp*p(i));
end
 
r = diag(1./((res.^2)));                    % Residual's matrix, standard deviation not available

stdev_p = 0.3*pvec;                         % A priori initial standard deviation  % Assumed to be 1/3
var_p = stdev_p.^2;  
s_p_0 = diag(1./var_p);                     % A priori covariance matrix
s_p = ( s_p_0 + J'*r*J)\eye(7);             % Covariance matrix
st_p = sqrt(diag(s_p));                     % Standard deviation on parameters


%% Monte Carlo Simulation

n = 1000;                          % Number of simulations
X_l = st_p.*randn(7,n) + pvec;     % Random parameters generation for left hand tip
X_r = st_p.*randn(7,n) + pvec;     % Random parameters generation for right hand tip
Tlag_l = 3*rand(1,n)*1e-5;         % Random lag time (30 microseconds max)
Tlag_r = 0;
mis_l = 100*randn(1,n)*1e-6;       % Misalignment of left tip w.r.t. TM baricenter  % um
mis_r = 100*randn(1,n)*1e-6;       % Misalignment of right tip w.r.t. TM baricenter % um

X_el_l = X_adh(1,:) + randn(n,3).*sigma_adh(1,:);
X_el_r = X_adh(1,:) + randn(n,3).*sigma_adh(1,:);

% Integrators testing
[A, ~,~,~,~,~] = sys2(p, 0);    
xt0 = A\[-120/p(3);0;0.3];
x0 = [xt0;xt0;0;0;0];

% ODE23
tic
[T,X1] = ode23(@(t,x) rel(t,x,p,p,M,I,0,0,0,0,xt0,xt0,X_el_l(1,:),X_el_r(1,:),0),[0 0.0015],x0);
cputime1 = toc;

% ODE23s
tic
[T,X2] = ode23s(@(t,x) rel(t,x,p,p,M,I,0,0,0,0,xt0,xt0,X_el_l(1,:),X_el_r(1,:),0),[0 0.0015],x0);
cputime2 = toc;

% ODE23t
tic
[T,X3] = ode23t(@(t,x) rel(t,x,p,p,M,I,0,0,0,0,xt0,xt0,X_el_l(1,:),X_el_r(1,:),0),[0 0.0015],x0);
cputime3 = toc;

% ODE23tb
tic
[T,X4] = ode23tb(@(t,x) rel(t,x,p,p,M,I,0,0,0,0,xt0,xt0,X_el_l(1,:),X_el_r(1,:),0),[0 0.0015],x0);
cputime4 = toc;

fprintf('ode23: \n     No. of time steps = %d, \t  CPUtime = %f \n',length(X1),cputime1)
fprintf('ode23s: \n    No. of time steps = %d, \t  CPUtime = %f \n',length(X2),cputime2)
fprintf('ode23t: \n    No. of time steps = %d, \t  CPUtime = %f \n',length(X3),cputime3)
fprintf('ode23tb: \n   No. of time steps = %d, \t  CPUtime = %f \n',length(X3),cputime3)

% Initialize and perform simulation

V_res = zeros(1,n);                % Residual velocity
W_res = zeros(1,n);
h = waitbar(0,'Monte Carlo Simulation');

t_end = 0.00015;                    % Simulation end time

for i = 1:n
    waitbar(i/n)
    
    % Determine initial conditions
    [Al, Bl,~,~,~,~] = sys2(X_l(:,i), 0);
    [Ar, Br,~,~,~,~] = sys2(X_r(:,i), 0);
    xl0 = Al\[-120/X_l(3,i);0;0.3];
    xr0 = Ar\[-120/X_r(3,i);0;0.3];
    x0 = [xl0;xr0;0;0;0];
    
    % Perform simulation
    [t_s,x_s] = ode23t(@(t,x) rel(t,x,X_l(:,i),X_r(:,i),M,I,Tlag_l(i),Tlag_r,mis_l(i),mis_r(i),xl0,xr0,X_el_l(i,:),X_el_r(i,:),0),[0,t_end],x0);
    
    % Evaluate results
    V_res(i) = x_s(end,8);  %residual linear velocity at each simulation
    W_res(i) = x_s(end,9);  %residual angular velocity at each simulation
end

close(h)

% Statistical analysis
figure(5)
h = histogram((V_res),'Normalization','Probability','NumBins',50);
title('Residual TM linear velocity PDF')
xlabel('Residual velocity [m/s]')

figure(6)
k = histogram(abs(W_res),'Normalization','Probability','NumBins',50);

% Compute 3 sigma residual velocity
H = 0;
i = 1;

while H<0.997 && i<n
    H = H + h.Values(i);
    i = i + 1;
end

V_res_3s = h.BinEdges(i-1);         % 3 sigma residual linear velocity

K = 0;
i = 1;

while K<0.997 && i<n
    K = K + k.Values(i);
    i = i + 1;
end

W_res_3s = h.BinEdges(i-1);         % 3 sigma residual angular velocity

title('Residual TM angular velocity PDF')
xlabel('Residual angular velocity [rad/s]')

%% Simulation with noise evaluation
