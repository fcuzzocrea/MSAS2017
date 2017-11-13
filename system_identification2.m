clearvars
close all
%% linear system identification

experimental_responce;      %load data

% resample data
T_fin = time_experimental(end); %final time [ms]
Ts = 1e-2;      %sampling time [ms]
T_sim = 0:Ts:T_fin; %resampled experiment time vector [ms]
y_exp = interp1(time_experimental,displacement_experimental,T_sim,'spline');
n=6;
input_voltage = [120*ones(1,n), zeros(1,length(T_sim)-n)];

%experimental data structure
data = iddata(y_exp', input_voltage', Ts);

%identification
est = tfest(data,5,2);
figure(2)
opt = compareOptions('InitialCondition','Estimate');
compare(data,est,Inf)
v = getpvec(est);
w = [-v(3:-1:1);v(8:-1:4)];
prm = fsolve(@(x) param_fun(x,w),ones(8,1));