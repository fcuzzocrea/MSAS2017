%{
   Modeling and Simlation of Aerospace Systems, Fall 2017
   Politecnico di Milano
%}

% Script files required to run main: sys2.m, odefun.m, rel.m
% Workspace files required: adhesion.mat, responce_ref2.mat

clearvars
close all
clear 
clc

% preliminary data
M = 1.9;            %TM mass in [kg]
I = 6.9e-4;         %inertia [kg*m^2] w.r.t. any of the 3 axes
rng('default')
%% Compute adhesion forces fitting
% retrieve experimental data

load('adhesion')            %load adhesion experimental data:
                            %ydata = adhesion force
                            %xdata = elongation
                            
% data fitting
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
%fitting plots
%set 1
l = linspace(0,3,100);
figure
plot(l*1e-6,X_adh(1,1)*l.*exp(-X_adh(1,2)*l.^(X_adh(1,3))),'r',set(1).exp(1).el*1e-6,set(1).exp(1).F,'.r',set(1).exp(2).el*1e-6,set(1).exp(2).F,'.r'...
    ,set(1).exp(3).el*1e-6,set(1).exp(3).F,'.r',set(1).exp(4).el*1e-6,set(1).exp(4).F,'.r',set(1).exp(5).el*1e-6,set(1).exp(5).F,'.r')
grid on
hold on
%set 2
plot(l*1e-6,X_adh(2,1)*l.*exp(-X_adh(2,2)*l.^(X_adh(2,3))),'g',set(2).exp(1).el*1e-6,set(2).exp(1).F,'.g',set(2).exp(2).el*1e-6,set(2).exp(2).F,'.g'...
    ,set(2).exp(3).el*1e-6,set(2).exp(3).F,'.g',set(2).exp(4).el*1e-6,set(2).exp(4).F,'.g',set(2).exp(5).el*1e-6,set(2).exp(5).F,'.g')
%set 3
plot(l*1e-6,X_adh(3,1)*l.*exp(-X_adh(3,2)*l.^(X_adh(3,3))),'b',set(3).exp(1).el*1e-6,set(3).exp(1).F,'.b',set(3).exp(2).el*1e-6,set(3).exp(2).F,'.b'...
    ,set(3).exp(3).el*1e-6,set(3).exp(3).F,'.b',set(3).exp(4).el*1e-6,set(3).exp(4).F,'.b',set(3).exp(5).el*1e-6,set(3).exp(5).F,'.b')
grid on
title('Experimental and fitted data')
xlabel('Elongation [m]')
ylabel('Force [N]')
axis([0 3e-6 0 0.16]);

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
sys.structure.Parameters.Minimum(7) = -0.05;%  "   "    "  V0

sys.structure.Parameters.Maximum(1) = 1e-6; %max value for C
sys.structure.Parameters.Maximum(3) = 800;  % "    "    "  R
sys.structure.Parameters.Maximum(4) = 0.014;% "    "    "  m
sys.structure.Parameters.Maximum(7) = 0.05; % "    "    "  V0


% load and resample data to fit in identification
load('responce_ref2.mat')
tf = t_exp(end);
ts = 1e-6;          %sampling time (f=1MHz)
t_sim = t_exp(1):ts:tf;
y_exp = interp1(t_exp,d_exp,t_sim,'linear');
n = 4;              %lag time (in time samples)
input = [zeros(1,n),120*ones(1,length(t_sim)-n)];
data = iddata(y_exp',input',ts);

% identification
opt = greyestOptions('SearchMethod','lm','InitialState','model','DisturbanceModel','none');
id_sys = pem(data,sys,opt);

p = getpvec(id_sys);

% compare results to validate
[A,~,C,~,~,x0] = sys2(p,0);
figure
subplot(2,1,1)
opt = compareOptions('InitialCondition',x0);
PEM = id_sys;
compare(data,PEM,Inf,opt)
legend('location','NorthWest')
title('experimental vs. model responce')
ylabel('elongation [m]')
grid on



%% Simulation with parameters covariance matrix estimation
 
% residuals computation
[t_i,x_i] = ode113(@(t,x) odefun(t,x,p),t_sim,x0);
y_i = C*x_i';                    %modeled tip responce
res = y_i - y_exp;               %model error w.r.t. experiment
s_y = res*res'/length(res);      %variance of the residual 
y_rms = sqrt(s_y);               %r.m.s. error
subplot(2,1,2)
plot(t_i,res)
grid on
title('residuals')
xlabel('Time (seconds)')
ylabel('residual value [m]')
r = diag(1./((res.^2)));

% stability analysis
l = eig(A);
figure
plot(real(l),imag(l),'x')
stiff = max(abs(real(l)))/min(abs(real(l)));
title('eigenvalues location')
xlabel('Re')
ylabel('Im')
grid on


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
    [T,X] = ode113(@(t,x) odefun(t,x,p1),t_sim,x0);
    y = C*X';
    J(:,i) = (y-y_i)./(dp*p(i));
end

stdev_p = 0.3*p;
var_p = stdev_p.^2;
s_p = (diag(1./var_p) + J'*r*J)\eye(7);     %covariance matrix
st_p = sqrt(diag(s_p));                     %standard deviation on parameters

%% Monte Carlo Simulation
n = 5000;                      %number of simulations
X_l = st_p.*randn(7,n) + p;    %random parameters generation for left hand tip
X_r = st_p.*randn(7,n) + p;    %random parameters generation for right hand tip
max_lag_time = 3e-5;           %random lag time (30 microseconds max)
Tlag_l = rand(1,n)*max_lag_time;
Tlag_r = 0;
mis_l = 100*randn(1,n)*1e-6;    %misalignment of left tip w.r.t. TM baricenter
mis_r = 100*randn(1,n)*1e-6;    %misalignment of right tip w.r.t. TM baricenter

%adhesion coefficients 1st set
X_el_l1 = X_adh(1,:) + randn(n,3).*sigma_adh(1,:);
X_el_r1 = X_adh(1,:) + randn(n,3).*sigma_adh(1,:);

%adhesion coefficients 2nd set
X_el_l2 = X_adh(2,:) + randn(n,3).*sigma_adh(2,:);
X_el_r2 = X_adh(2,:) + randn(n,3).*sigma_adh(2,:);

%adhesion coefficients 3rd set
X_el_l3 = X_adh(3,:) + randn(n,3).*sigma_adh(3,:);
X_el_r3 = X_adh(3,:) + randn(n,3).*sigma_adh(3,:);

% integrators testing
cputime1 = zeros(1,50);
cputime2 = zeros(1,50);
cputime3 = zeros(1,50);
cputime4 = zeros(1,50);
step1    = zeros(1,50);
step2    = zeros(1,50);
step3    = zeros(1,50);
step4    = zeros(1,50);
h = waitbar(0,'integrators testing');
for i = 1:50
    waitbar(i/50)
    [Al, ~,~,~,~,~] = sys2(X_l(:,i), 0);
    [Ar, ~,~,~,~,~] = sys2(X_r(:,i), 0);
    xt0_l = Al\[-120/X_l(3,i);0;0.3];
    xt0_r = Ar\[-120/X_r(3,i);0;0.3];
    x0 = [xt0_l;xt0_r;0;0;0];

% ode23
    tic
    [T,X1] = ode23(@(t,x) rel(t,x,X_l,X_r,M,I,Tlag_l(i),0,mis_l(i),mis_r(i),xt0_l,xt0_r,X_el_l1(i,:),X_el_r1(i,:),0),[0 0.0015],x0);
    cputime1(i) = toc;
    step1(i) = length(T);
% ode23s
    tic
    [T,X2] = ode23s(@(t,x) rel(t,x,X_l,X_r,M,I,Tlag_l(i),0,mis_l(i),mis_r(i),xt0_l,xt0_r,X_el_l1(i,:),X_el_r1(i,:),0),[0 0.0015],x0);
    cputime2(i) = toc;
    step2(i) = length(T);
%ode23t
    tic
    [T,X3] = ode23t(@(t,x) rel(t,x,X_l,X_r,M,I,Tlag_l(i),0,mis_l(i),mis_r(i),xt0_l,xt0_r,X_el_l1(i,:),X_el_r1(i,:),0),[0 0.0015],x0);
    cputime3(i) = toc;
    step3(i) = length(T);
%ode23tb
    tic
    [T,X4] = ode23tb(@(t,x) rel(t,x,X_l,X_r,M,I,Tlag_l(i),0,mis_l(i),mis_r(i),xt0_l,xt0_r,X_el_l1(i,:),X_el_r1(i,:),0),[0 0.0015],x0);
    cputime4(i) = toc;
    step4(i) = length(T);
end
close(h)

% print results
cputime23     = sum(cputime1)/length(cputime1);
cputime23s    = sum(cputime2)/length(cputime2);
cputime23t    = sum(cputime3)/length(cputime3);
cputime23tb   = sum(cputime4)/length(cputime4);
timesteps23   = round(sum(step1)/length(step1));
timesteps23s  = round(sum(step2)/length(step2));
timesteps23t  = round(sum(step3)/length(step3));
timesteps23tb = round(sum(step4)/length(step4));

fprintf(  'ode23: \n     No. of time steps (mean) = %d, \t  CPUtime = %f \n',timesteps23,cputime23)
fprintf( 'ode23s: \n     No. of time steps (mean) = %d, \t  CPUtime = %f \n',timesteps23s,cputime23s)
fprintf( 'ode23t: \n     No. of time steps (mean) = %d, \t  CPUtime = %f \n',timesteps23t,cputime23t)
fprintf('ode23tb: \n     No. of time steps (mean) = %d, \t  CPUtime = %f \n',timesteps23tb,cputime23tb)

% initialize and perform simulation
V_res = zeros(3,n);
W_res = zeros(3,n);
h = waitbar(0,'Monte Carlo Simulation (with parameters variance)');
t_end = 0.00015;                %simulation end time
for i = 1:n
    waitbar(i/n)
    % determine initial conditions
    [Al, Bl,~,~,~,~] = sys2(X_l(:,i), 0);
    [Ar, Br,~,~,~,~] = sys2(X_r(:,i), 0);
    xl0 = Al\[-120/X_l(3,i);0;0.3];
    xr0 = Ar\[-120/X_r(3,i);0;0.3];
    x0 = [xl0;xr0;0;0;0];
    
    % perform simulation
    % adhesion set 1
    [t_s,x_s] = ode23tb(@(t,x) rel(t,x,X_l(:,i),X_r(:,i),M,I,Tlag_l(i),Tlag_r,mis_l(i),mis_r(i),xl0,xr0,X_el_l1(i,:),X_el_r1(i,:),0),[0,t_end],x0);
    % evaluate results
    V_res(1,i) = x_s(end,8);  %residual linear velocity at each simulation
    W_res(1,i) = x_s(end,9);  %residual angular velocity at each simulation
    % adhesion set 2
    [t_s,x_s] = ode23tb(@(t,x) rel(t,x,X_l(:,i),X_r(:,i),M,I,Tlag_l(i),Tlag_r,mis_l(i),mis_r(i),xl0,xr0,X_el_l2(i,:),X_el_r2(i,:),0),[0,t_end],x0);
    V_res(2,i) = x_s(end,8);  %residual linear velocity at each simulation
    W_res(2,i) = x_s(end,9);  %residual angular velocity at each simulation
    %adhesion set 3
    [t_s,x_s] = ode23tb(@(t,x) rel(t,x,X_l(:,i),X_r(:,i),M,I,Tlag_l(i),Tlag_r,mis_l(i),mis_r(i),xl0,xr0,X_el_l3(i,:),X_el_r3(i,:),0),[0,t_end],x0);
    V_res(3,i) = x_s(end,8);  %residual linear velocity at each simulation
    W_res(3,i) = x_s(end,9);  %residual angular velocity at each simulation
end
close(h)

% statistical analysis
figure
hold on
h1 = histogram(V_res(1,:),'Normalization','Probability','NumBins',50,'BinLimits',[-2e-6,7e-6]);
h2 = histogram(V_res(2,:),'Normalization','Probability','NumBins',50,'BinLimits',[-2e-6,7e-6]);
h3 = histogram(V_res(3,:),'Normalization','Probability','NumBins',50,'BinLimits',[-2e-6,7e-6]);
title('residual TM linear velocity PDF')
xlabel('residual velocity [m/s]')
legend('set 1','set 2','set 3')

figure
hold on
k1 = histogram(abs(W_res(1,:)),'Normalization','Probability','NumBins',50,'BinLimits',[0, 3e-6]);
k2 = histogram(abs(W_res(2,:)),'Normalization','Probability','NumBins',50,'BinLimits',[0, 3e-6]);
k3 = histogram(abs(W_res(3,:)),'Normalization','Probability','NumBins',50,'BinLimits',[0, 3e-6]);
title('residual TM angular velocity PDF')
xlabel('residual angular velocity [rad/s]')
legend('set 1','set 2','set 3')
% compute 3 sigma residual velocity
V_res_3s = zeros(3,1);
W_res_3s = zeros(3,1);
% adhesion set 1
H = 0;
i = 1;
while H<0.997 && i<length(h1.BinEdges)
    H = H + h1.Values(i);
    i = i + 1;
end
V_res_3s(1) = h1.BinEdges(i-1);         %3 sigma residual linear velocity

K = 0;
i = 1;
while K<0.997 && i<n
    K = K + k1.Values(i);
    i = i + 1;
end
W_res_3s(1) = k1.BinEdges(i-1);         %3 sigma residual angular velocity

% adhesion set 2
H = 0;
i = 1;
while H<0.997 && i<n
    H = H + h2.Values(i);
    i = i + 1;
end
V_res_3s(2) = h2.BinEdges(i-1);         %3 sigma residual linear velocity

K = 0;
i = 1;
while K<0.997 && i<n
    K = K + k2.Values(i);
    i = i + 1;
end
W_res_3s(2) = k2.BinEdges(i-1);         %3 sigma residual angular velocity

% adhesion set 3
H = 0;
i = 1;
while H<0.997 && i<n
    H = H + h3.Values(i);
    i = i + 1;
end
V_res_3s(3) = h3.BinEdges(i-1);         %3 sigma residual linear velocity

K = 0;
i = 1;
while K<0.997 && i<n
    K = K + k3.Values(i);
    i = i + 1;
end
W_res_3s(3) = k3.BinEdges(i-1);         %3 sigma residual angular velocity

%% Simulation with noise evaluation
% weighting filter characterization
C_rr= conv(res, res);
PSD_2sided = fft(C_rr);
d = ceil(length(PSD_2sided)/2);
PSD = sqrt(2*PSD_2sided(1:d));
figure
fr = linspace(1,5e5,length(PSD));
loglog(fr,real(PSD))
title('post-fit residuals PSD')
xlabel('Freq. [Hz]')
ylabel('PSD [m/sqrt(Hz)]')
grid on

% System FRF
c   = p(1);
Tem = p(2);
R   = p(3);
m   = p(4);
b   = p(5);
k   = p(6);

H =@(w) (R*c*1i*w + 1)./(k + b*1i*w - m*w.^2 + R*Tem^2*1i*w - R*b*c*w.^2 - R*c*m*1i*w.^3 + R*c*k*1i.*w);

% introduce white noise
N = fft(wgn(length(PSD_2sided),1,0,'dbW'));
N = N(1:d)';

% residual and force
x = N.*PSD;
F_n = x./H(fr);
figure
subplot(3,1,1)
loglog(fr,abs(H(fr)))
xlabel('frequency [Hz]')
ylabel('|H|')
grid on
title('FRF')
subplot(3,1,2)
loglog(fr,abs(x))
grid on
title('simulated residual')
ylabel('residual [m]')
xlabel('frequency [Hz]')
subplot(3,1,3)
loglog(fr,abs(F_n))
grid on
title('simulated force in Fourier domain')
ylabel('Force [N]')
xlabel('frequency [Hz]')
%residual in time domain
X = ifft([x,x(end:-1:1)],'symmetric');
figure
plot(linspace(0,4e-4,length(X)),abs(X))
title('simulated residual in time domain')
xlabel('t [s]')
ylabel('res [m]')

% noise force in time domain
f_n = ifft([F_n,F_n(end:-1:1)],'symmetric');        %force as function of time
figure
plot(linspace(0,4e-4,length(f_n)),f_n)
title('simulated random force in time domain')
xlabel('t [s]')
ylabel('force [N]')
grid on

%% Montecarlo simulation: noise force approach
Noise = wgn(length(PSD_2sided),n,0,'dbW');    %white noise matrix

% initial conditions
[A, B,~,~,~,~] = sys2(p, 0);
xt0 = A\[-120/p(3);0;0.3];
x0  = [xt0;xt0;0;0;0];

V_res_noise = zeros(3,n);
W_res_noise = zeros(3,n);

% simulation
h = waitbar(0,'Monte-Carlo Simulation (with noise filtering)');
for i = 1:n
    waitbar(i/n)
    N_2sided = fft(Noise(:,i));
    N = N_2sided(1:d)';
    X = N.*PSD;
    F_n = X./H(linspace(0,5e5,length(X)));
    f_n = ifft(F_n,'symmetric');
    
    % adhesion set1
    [t_s,x_s] = ode23tb(@(t,x) rel(t,x,p,p,M,I,Tlag_l(i),0,mis_l(i),mis_r(i),xt0,xt0,X_el_l1(i,:),X_el_r1(i,:),f_n),[0,t_end],x0);
    V_res_noise(1,i) = x_s(end,8);
    W_res_noise(1,i) = x_s(end,9);
    
    % adhesion set1
    [t_s,x_s] = ode23tb(@(t,x) rel(t,x,p,p,M,I,Tlag_l(i),0,mis_l(i),mis_r(i),xt0,xt0,X_el_l2(i,:),X_el_r2(i,:),f_n),[0,t_end],x0);
    V_res_noise(2,i) = x_s(end,8);
    W_res_noise(2,i) = x_s(end,9);
    
    % adhesion set1
    [t_s,x_s] = ode23tb(@(t,x) rel(t,x,p,p,M,I,Tlag_l(i),0,mis_l(i),mis_r(i),xt0,xt0,X_el_l3(i,:),X_el_r3(i,:),f_n),[0,t_end],x0);
    V_res_noise(3,i) = x_s(end,8);
    W_res_noise(3,i) = x_s(end,9);
end
close(h)
figure 
plot(t_s,x_s(:,8))

% statistical analysis: plot PDF and compute 3 sigma residual velocity
figure
hold on
hn_1 = histogram(V_res_noise(1,:),'Normalization','Probability','NumBins',50,'BinLimits',[-1e-6,7e-6]);
hn_2 = histogram(V_res_noise(2,:),'Normalization','Probability','NumBins',50,'BinLimits',[-1e-6,7e-6]);
hn_3 = histogram(V_res_noise(3,:),'Normalization','Probability','NumBins',50,'BinLimits',[-1e-6,7e-6]);
title('residual TM velocity PDF (through noise filtering)')
xlabel('residual velocity [m/s]')
legend('set 1','set 3','set 3')

% 3 sigma quantities
% linear velocity
V_res_noise_3s = zeros(1,3);
% set 1
H = 0;
i = 1;
while H<0.997 && i<length(hn_1.BinEdges)
    H = H + hn_1.Values(i);
    i = i + 1;
end
V_res_noise_3s(1) = hn_1.BinEdges(i-1);         %3 sigma residual linear velocity

% set 2
H = 0;
i = 1;
while H<0.997 && i<n
    H = H + hn_1.Values(i);
    i = i + 1;
end
V_res_noise_3s(2) = hn_2.BinEdges(i-1);         %3 sigma residual linear velocity

% set 3
H = 0;
i = 1;
while H<0.997 && i<n
    H = H + hn_3.Values(i);
    i = i + 1;
end
V_res_noise_3s(3) = hn_3.BinEdges(i-1);         %3 sigma residual linear velocity

% angular velocity
figure
hold on
kn_1 = histogram(W_res_noise(1,:),'Normalization','Probability','NumBins',50,'BinLimits',[-3e-6,3e-6]);
kn_2 = histogram(W_res_noise(2,:),'Normalization','Probability','NumBins',50,'BinLimits',[-3e-6,3e-6]);
kn_3 = histogram(W_res_noise(3,:),'Normalization','Probability','NumBins',50,'BinLimits',[-3e-6,3e-6]);
title('residual TM angular velocity PDF (through noise filtering)')
xlabel('residual angular velocity [rad/s]')
legend('set 1','set 3','set 3')

% set 1
K = 0;
i = 1;
while K<0.997 && i<length(kn_1.Values)
    K = K + kn_1.Values(i);
    i = i + 1;
end
W_res_noise_3s(1) = kn_1.BinEdges(i-1);         %3 sigma residual angular velocity

% set 2
K = 0;
i = 1;
while K<0.997 && i<length(kn_2.Values)
    K = K + kn_2.Values(i);
    i = i + 1;
end
W_res_noise_3s(2) = kn_2.BinEdges(i-1);         %3 sigma residual angular velocity

% set 3
K = 0;
i = 1;
while K<0.997 && i<length(kn_3.Values)
    K = K + kn_3.Values(i);
    i = i + 1;
end
W_res_noise_3s(3) = kn_3.BinEdges(i-1);         %3 sigma residual angular velocity

%% sensitivity tests
V_res_test = zeros(3,n);
Tlag_l1 = rand(1,n)*2*max_lag_time;
Tlag_l2 = rand(1,n)*max_lag_time/2;
h = waitbar(0, 'sensitivity Monte-Carlo');
for i = 1:n
    waitbar(i/n)
    % determine initial conditions
    [Al, Bl,~,~,~,~] = sys2(X_l(:,i), 0);
    [Ar, Br,~,~,~,~] = sys2(X_r(:,i), 0);
    xl0 = Al\[-120/X_l(3,i);0;0.3];
    xr0 = Ar\[-120/X_r(3,i);0;0.3];
    x0 = [xl0;xr0;0;0;0];
    
    % perform simulation
    [t_s,x_s] = ode23tb(@(t,x) rel(t,x,X_l(:,i),X_r(:,i),M,I,Tlag_l(i),Tlag_r,mis_l(i),mis_r(i),xl0,xr0,X_el_l1(i,:),X_el_r1(i,:),0),[0,t_end],x0);
    % evaluate results
    V_res_test(1,i) = x_s(end,8);  %residual linear velocity at each simulation
    
    % perform simulation
    [t_s,x_s] = ode23tb(@(t,x) rel(t,x,X_l(:,i),X_r(:,i),M,I,Tlag_l1(i),Tlag_r,mis_l(i),mis_r(i),xl0,xr0,X_el_l1(i,:),X_el_r1(i,:),0),[0,t_end],x0);
    % evaluate results
    V_res_test(2,i) = x_s(end,8);  %residual linear velocity at each simulation
    
    % perform simulation
    [t_s,x_s] = ode23tb(@(t,x) rel(t,x,X_l(:,i),X_r(:,i),M,I,Tlag_l2(i),Tlag_r,mis_l(i),mis_r(i),xl0,xr0,X_el_l1(i,:),X_el_r1(i,:),0),[0,t_end],x0);
    % evaluate results
    V_res_test(3,i) = x_s(end,8);  %residual linear velocity at each simulation
end
close(h)
figure
hold on
l2 = histogram(V_res_test(3,:),'Normalization','Probability','NumBins',50,'BinLimits',[-1e-6,10e-6]);
l  = histogram(V_res_test(1,:),'Normalization','Probability','NumBins',50,'BinLimits',[-1e-6,10e-6]);
l1 = histogram(V_res_test(2,:),'Normalization','Probability','NumBins',50,'BinLimits',[-1e-6,10e-6]);
title('residual TM linear velocity PDF with different max time lag')
xlabel('residual linear velocity [m/s]')
legend('15e-6 [s]','30e-6 [s]','60e-6 [s]')