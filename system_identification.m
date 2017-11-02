clearvars
close all
%% linear system identification

%model data structure
% the parameter guess are in the order: m[kg], mp[kg], k[N/micron], kp[N/micron]
%, b[N*ms/micron], bp[N*ms/micron], c[F], Tem[C/micron], R[ohm]

par = [1e-3, 1e-3, 1e-3, 1e1, 1e-2, 0.1, 4.8*1e-7 , 1e-6, 400];   %set of parameters initial guess
aux = {};           %estimation optional arguments
T = 0;              %parameter to specify model as continous-time
m = idgrey('par_system',par,'c',aux,T);     
%constraints on parameters
m.Structure.Parameters.Minimum = [1e-5; 1e-5; 1e-4; 1e-1; 1e-4; 1e-2; 1e-9; 1e-8; 350];
m.Structure.Parameters.Maximum = [1e-2; 1e-2; 1e-2; 10; 1; 1; 1e-4; 1e-4; 450];
experimental_responce;      %load data

% resample output
T_fin = time_experimental(end); %final time in [ms]
Ts = 1e-3;      %sampling time in [ms], consistent with sample frequency = 1MHz
T_sim = 0:Ts:T_fin; %resampled experiment time vector [ms]
y_exp = interp1(time_experimental,x_experimental,T_sim,'spline');
n=6;
input_voltage = [120*ones(1,n), zeros(1,length(T_sim)-n)];
% input_voltage = zeros(size(T_sim));

%experimental data structure
data = iddata(y_exp', input_voltage', Ts);

%estimation
opt = greyestOptions('SearchMethod','lm','InitialState','model','DisturbanceModel','none');
par_est = greyest(data,m,opt);

p = getpvec(par_est);
s = diag(getcov(par_est));
e = sqrt(s)./p;

[~,~,~,~,~,x0] = par_system(p,0);
opt = compareOptions('InitialCondition',x0);
figure(2)
compare(data,par_est,Inf,opt)