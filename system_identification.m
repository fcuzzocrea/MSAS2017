clearvars
close all
%% linear system identification

% the parameters guess are in the order: m[kg], c[F], Tem[C/micron], kp[N/micron]
%bp[N*ms/micron], mp[kg], b[N*ms/micron], R[ohm]
par = [1e-3, 1e-7, 1e-6, 10, 1e-1, 1e-3, 1e-3, 1e-2, 400];   %set of parameters initial guess
aux = {};           %estimation optional arguments
T = 0;              %parameter to specify model as continous-time
m = idgrey('par_system',par,'c',aux,T);

experimental_responce;      %load data

% resample output
T_fin = time_experimental(end); %final time in [ms]
Ts = 1e-3;      %sampling time in [ms]
T_sim = 0:Ts*10:T_fin; %resampled experiment time vector [ms]
y_exp = interp1(time_experimental,displacement_experimental,T_sim,'spline');
n=6;
input_voltage = [120*ones(1,n), zeros(1,length(T_sim)-n)];
% input_voltage = zeros(size(T_sim));

%experimental data structure
data = iddata(y_exp', input_voltage', Ts);

%estimation
opt = greyestOptions('SearchMethod','lm','InitialState','estimate');
par_est = greyest(data,m,opt);
opt = compareOptions('InitialCondition','Estimate');
figure(2)
compare(data,par_est,Inf,opt)