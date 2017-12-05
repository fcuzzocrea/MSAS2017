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
pvec = getpvec(id_sys);                                 % Parameters of the identified system

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

load('adhesion_dataset_1')      % Load adhesion force experimental data:
                                %
                                % ydata = Force
                                % xdata = Elongation
                            
% Data Fitting
[X_el,sigma] = lsqcurvefit(@(x,xdata) x(1)*xdata.*exp(-x(2)*(xdata.^x(3))),ones(3,1),xdata,ydata);


%% Simulation with parameters covariance matrix estimation

% Residuals computation
[t_i,x_i] = ode113(@(t,x)odefun(t,x,pvec),t_sim,x0);
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
l = eig(A);
figure(3)
plot(real(l),imag(l),'x')
stiff = max(abs(real(l)))/min(abs(real(l)));
title('Eigenvalues location')
xlabel('Re')
ylabel('Im')
grid on

% Jacobian estimation:
% This estimation is performed by simulating the responce varing one
% single parameter each time and comparing it with the nominal parameters
% system responce

J = zeros(length(y_i),length(pvec));
dp = 0.01;                                  % ? Perchè ? 

for i = 1:length(pvec)
    p1 = pvec;
    p1(i) = pvec(i)*(1+dp);
    [~,~,C,~,~,x0] = sys2(p1,0);
    [T,X] = ode113(@(t,x) odefun(t,x,p1),t_sim,x0);
    y = C*X';
    J(:,i) = (y-y_i)./(dp*pvec(i));
end
 
r = diag(1./((res.^2)));                    % Residual's matrix   

stdev_p = 0.3*pvec;                         % A priori initial standard deviation  % ? Perchè ?
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
mis_l = 100*randn(1,n)*1e-6;       % Misalignment of left tip w.r.t. TM baricenter
mis_r = 100*randn(1,n)*1e-6;       % Misalignment of right tip w.r.t. TM baricenter
V_res = zeros(1,n);
W_res = zeros(1,n);
h = waitbar(0,'Monte Carlo Simulation');

for i = 1:n
    waitbar(i/n)
    % determine initial conditions
    [Al, Bl,~,~,~,~] = sys2(X_l(:,i), 0);
    [Ar, Br,~,~,~,~] = sys2(X_r(:,i), 0);
    xl0 = Al\[-120/X_l(3,i);0;0.3];
    xr0 = Ar\[-120/X_r(3,i);0;0.3];
    x0 = [xl0;xr0;0;0;0];
    
    % perform simulation
    [t_s,x_s] = ode23(@(t,x) rel(t,x,X_l(:,i),X_r(:,i),M,I,Tlag_l(i),Tlag_r,mis_l(i),mis_r(i),xl0,xr0,X_el),[0,0.00015],x0);
    % evaluate results
    V_res(i) = x_s(end,8);  %residual linear velocity at each simulation
    W_res(i) = x_s(end,9);  %residual angular velocity at each simulation
end
close(h)

% statistical analysis
figure
h = histogram((V_res),'Normalization','Probability','NumBins',50);
title('residual TM linear velocity PDF')
xlabel('residual velocity [m/s]')

figure
k = histogram(abs(W_res),'Normalization','Probability','NumBins',50);

% compute 3 sigma residual velocity
H = 0;
i = 1;
while H<0.997 && i<n
    H = H + h.Values(i);
    i = i + 1;
end
V_res_3s = h.BinEdges(i-1);         %3 sigma residual linear velocity

K = 0;
i = 1;
while K<0.997 && i<n
    K = K + k.Values(i);
    i = i + 1;
end
W_res_3s = h.BinEdges(i-1);         %3 sigma residual angular velocity
title('residual TM angular velocity PDF')
xlabel('residual angular velocity [rad/s]')

%% Simulation with noise evaluation
