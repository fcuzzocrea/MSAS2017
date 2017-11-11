%{
   Modeling and Simlation of Aerospace Systems, Fall 2017
   Politecnico di Milano
%}

% Script files required to run main: sys2.m, odefun.m, rel.m
% Workspace files required: adhesion_dataset_1.mat, responce_ref2.mat

clearvars
close all
clear 
clc


%% identification
% define system to identify
par = [4.8e-7, 1.5, 400, 1e-3, 150, 1.5e7, 0];
aux = {};           %estimation optional arguments
T = 0;              %parameter to specify model as continous-time
sys = idgrey('sys2',par,'c',aux,T);

% set constraints on parameters
sys.structure.Parameters.Minimum(1) = 0;    %min value for C
sys.structure.Parameters.Minimum(2) = 0;    %min value for Tem
sys.structure.Parameters.Minimum(3) = 100;  %min value for R
sys.structure.Parameters.Minimum(4) = 1e-4; %  "   "    "  m
sys.structure.Parameters.Minimum(5) = 100;  %  "   "    "  b
sys.structure.Parameters.Minimum(6) = 0;    %  "   "    "  k
sys.structure.Parameters.Minimum(7) = -0.05; %  "   "    "  V0

sys.structure.Parameters.Maximum(1) = 1e-6; %max value for C
sys.structure.Parameters.Maximum(3) = 800;  % "    "    "  R
sys.structure.Parameters.Maximum(4) = 0.014;% "    "    "  m
sys.structure.Parameters.Maximum(7) = 0.05;  % "    "    "  V0


% load and resample data to fit in identification
load('responce_ref2.mat')
tf = t_exp(end);
ts = 1e-6;          %sampling time (f=1MHz)
t_sim = t_exp(1):ts:tf;
y_exp = interp1(t_exp,d_exp,t_sim,'linear');
n = 4;
input = [zeros(1,n),120*ones(1,length(t_sim)-n)];
data = iddata(y_exp',input',ts);

% identification
opt = greyestOptions('SearchMethod','lm','InitialState','model','DisturbanceModel','none');
id_sys = pem(data,sys,opt);

p = getpvec(id_sys);

% compare results to validate
[A,B,C,~,~,x0] = sys2(p,0);
figure
subplot(2,1,1)
opt = compareOptions('InitialCondition',x0);
compare(data,id_sys,Inf,opt)
title('experimental vs. model responce')
ylabel('elongation [m]')
grid on


%% Compute adhesion forces fitting
% retrieve experimental data

load('adhesion_dataset_1')  %load adhesion experimental data:
                            %ydata = adhesion force
                            %xdata = elongation
                            
% data fitting
[X_el,sigma] = lsqcurvefit(@(x,xdata) x(1)*xdata.*exp(-x(2)*(xdata.^x(3))),ones(3,1),xdata,ydata);


%% simulation
% parameters covariance matrix estimation
% residuals computation
[t_i,x_i] = ode23(@(t,x) odefun(t,x,p),t_sim,x0);
y_i = C*x_i';           %modeled tip responce
res = y_i - y_exp;               %model error w.r.t. experiment
s_y = res*res'/length(res);
y_rms = sqrt(s_y);              %r.m.s. error
subplot(2,1,2)
plot(t_i,res)
grid on
title('residual') 
xlabel('time [s]')
ylabel('res')
r = diag(1./((res.^2)));


% Jacobian estimation:
% this estimation is performed by simulating the responce varing one
% single parameter each time and comparing it with the nominal parameters
% system responce

J = zeros(length(y_i),length(p));
dp = 0.01;

for i = 1:length(p)
    p1 = p;
    p1(i) = p(i)*(1+dp);
    [~,~,C,~,~,x0] = sys2(p1,0);
    [T,X] = ode23(@(t,x) odefun(t,x,p1),t_sim,x0);
    y = C*X';
    J(:,i) = (y-y_i)./(dp*p(i));
end

stdev_p = 0.3*p;
var_p = stdev_p.^2;
s_p = (diag(1./var_p)+J'*r*J)\eye(7);       %covariance matrix
st_p = sqrt(diag(s_p));                     %standard deviation on parameters

%% Monte Carlo Simulation
n = 10000;                     %number of simulations
X_l = st_p.*randn(7,n) + p;    %random parameters generation for left hand tip
X_r = st_p.*randn(7,n) + p;    %random parameters generation for right hand tip
Tlag_l = 3*rand(1,n)*1e-5;     %random lag time (30 microseconds max)
Tlag_r = 0;
M = 1.9;
V_res = zeros(1,n);
h = waitbar(0,'Monte Carlo Simulation');

for i = 1:n
    waitbar(i/n)
    % determine initial conditions
    [Al, Bl,~,~,~,~] = sys2(X_l(:,i), 0);
    [Ar, Br,~,~,~,~] = sys2(X_r(:,i), 0);
    xl0 = Al\[-120/X_l(3,i);0;0.3];
    xr0 = Ar\[-120/X_r(3,i);0;0.3];
    x0 = [xl0;xr0;0;0];
    
    % perform simulation
    [t_s,x_s] = ode23(@(t,x) rel(t,x,X_l(:,i),X_r(:,i),M,Tlag_l(i),Tlag_r,xl0,xr0,X_el),[0,0.00015],x0);
    % evaluate results
    V_res(i) = x_s(end,8);  %residual velocity at each simulation
end
close(h)

% statistical analysis
figure
h = histogram(V_res,'Normalization','Probability','NumBins',500);

% compute 3 sigma residual velocity
H = 0;
i = 1;
while H<0.997 && i<n
    H = H + h.Values(i);
    i = i + 1;
end
V_res_3s = h.BinEdges(i-1);         %3 sigma residual velocity

% plot histogram
h = histogram(V_res,'Normalization','Probability','NumBins',100);
title('residual TM velocity PDF')
xlabel('residual velocity [m/s]')