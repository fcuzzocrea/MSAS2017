clearvars
close all
%% linear system identification

%load and prepare data to fit
experimental_responce;
T_fin = time_experimental(end);     %final time in [s]
Ts = 1e-6;                          %sampling time in [s], consistent with sample frequency = 1MHz
T_sim = 0:Ts:T_fin;                 %resampled experiment time vector [s]
y_exp = interp1(time_experimental,x_experimental,T_sim,'linear');
n=70;
input_voltage = [120*ones(1,n), zeros(1,length(T_sim)-n)];
for i = 1:length(y_exp)
    if isnan(y_exp(i))
        y_exp(i) = 0;
    end
end

data = iddata(y_exp', input_voltage', Ts);

%random guess parameters identification
par = [5e-4, 1e-3, 1000, 1e7, 10, 100, 4.8e-7, 1.229, 400, 0];   %set of nominal parameters
% x = rand(100,length(par_nom));
% for i = 1:100
% par = par_nom.*(2*x(i));    %random initial guess parameters vector
aux = {};           %model optional arguments
T = 0;              %parameter to specify model as continous-time

%model data structure:
% the parameter guess are in the order: m[kg], mp[kg], k[N/m], kp[N/m]
%, b[N*s/m], bp[N*s/m], c[F], Tem[C/m], R[ohm]

sys = idgrey('par_system',par,'c',aux,T);     

%constraints on parameters
sys.Structure.Parameters.Free(10) = false;

sys.Structure.Parameters.Minimum(1) = 0;   %min value for m
sys.Structure.Parameters.Minimum(2) = 0;   %min value for mp
sys.Structure.Parameters.Minimum(3) = 0;      %min value for k
sys.Structure.Parameters.Minimum(4) = 0;      % "     "   "  kp
sys.Structure.Parameters.Minimum(5) = 0;        % "     "   "  b
sys.Structure.Parameters.Minimum(6) = 0;        % "     "   "  bp
sys.Structure.Parameters.Minimum(7) = 0;        % "     "   "  c
sys.Structure.Parameters.Minimum(8) = 0;        % "     "   "  Tem
sys.Structure.Parameters.Minimum(9) = 0;        % "     "   "  R
sys.Structure.Parameters.Minimum(10)= -1e-6;     % "     "   "  V0

sys.Structure.Parameters.Maximum(1) = 0.014;    %max value for m
sys.Structure.Parameters.Maximum(2) = 0.014;    % "     "   "  mp
sys.Structure.Parameters.Maximum(7) = 1e-6;     % "     "   "  c
sys.Structure.Parameters.Maximum(9) = 800;      % "     "   "  R
sys.Structure.Parameters.Maximum(10)= 1e-6;      % "     "   "  V0


%estimation
opt = greyestOptions('SearchMethod','lm','InitialState','model','DisturbanceModel','none');
id_par = greyest(data,sys,opt);

p = getpvec(id_par);
s = diag(getcov(id_par));
e = sqrt(s)./p;

[~,~,~,~,~,x0] = par_system(p,0);
opt = compareOptions('InitialCondition',x0);
figure(2)
compare(data,id_par,Inf,opt)