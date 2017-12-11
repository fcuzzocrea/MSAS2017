clearvars
close all
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

sys = init(sys,0.01*ones(size(par)).*par,par,'b');

% load and define data to fit
responce_ref2
tf = t_exp(end);
ts = 1e-6;          %sampling time (f=1MHz)
t_sim = t_exp(1):ts:tf;
y_exp = interp1(t_exp,d_exp,t_sim,'linear');
n = 4;
input = [zeros(1,n),120*ones(1,length(t_sim)-n)];
data = iddata(y_exp',input',ts);

%system identification
opt = greyestOptions('SearchMethod','gn','InitialState','model','DisturbanceModel','none');
id_sys = pem(data,sys,opt);

p = getpvec(id_sys);
s = diag(getcov(id_sys));
e = sqrt(s)./p;

%compare results
[A,B,C,~,~,x0] = sys2(p,0);
opt = compareOptions('InitialCondition',x0);
figure(1)
compare(data,id_sys,Inf,opt)

%integration
% [t_i,x_i] = ode45(@(t,x) odefun(t,x,A,B,input),t_sim,x0);
