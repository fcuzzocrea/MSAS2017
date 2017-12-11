system_id3
[t_i,x_i] = ode23(@(t,x) odefun(t,x,p),t_sim,x0);
y_i = C*x_i';
figure(2)
subplot(2,1,1)
plot(t_i,y_i,t_sim,y_exp,'--')
title('experimental vs. model responce')
legend('modeled','experimental','location','northwest')
ylabel('elongation [m]')
grid on
res = y_i-y_exp;
s_y = res*res'/length(res);
y_rms = sqrt(s_y);
subplot(2,1,2)
plot(t_i,res)
title('residual') 
xlabel('time [s]')
ylabel('res')
r = diag(1./((res.^2)));

% parameters covariance estimation

%Jacobian
J = zeros(length(y_i),length(p));
dp = 0.01;
% s_yp = zeros(size(p));

for i = 1:length(p)
    p1 = p;
    p1(i) = p(i)*(1+dp);
    [~,~,C,~,~,x0] = sys2(p1,0);
    [T,X] = ode45(@(t,x) odefun(t,x,p1),t_sim,x0);
    y = C*X';
%     s_yp(i) = ((y-y_i)*(y-y_i)')/length(res);
    J(:,i) = (y-y_i)./(dp*p(i));
end

% s_p = s_y./s_yp;
% s_p = 1./diag(J'*r*J);

stdev_p = 0.3*p;
var_p = stdev_p.^2;
s_p = (diag(1./var_p)+J'*r*J)\eye(7);
st_p = sqrt(diag(s_p));

% C_res = conv([res,res,res],res);
% C_res = C_res(floor(length(C_res)/2)-ceil(length(res)/2)+1:floor(length(C_res)/2)+ceil(length(res)/2));
% PSD_res = fftshift(fft(C_res));
% w = -0.5/ts:1/t_sim(end):0.5/ts;
% figure
% loglog(w,abs(PSD_res))

%% Compute adhesion forces fitting
% retrieve experimental data

adhesion_dataset_1
xdata = deltal;
ydata = adhesion_force;

% data fitting
[X_el,sigma] = lsqcurvefit(@(x,xdata) x(1)*xdata.*exp(-x(2)*(xdata.^x(3))),ones(3,1),xdata,ydata);

%% Monte Carlo Simulation
n = 10000;                      %number of simulations
X_l = st_p.*randn(7,n) + p;    %random parameters generation for left hand tip
X_r = st_p.*randn(7,n) + p;    %random parameters generation for right hand tip
Tlag_l = 3*rand(1,n)*1e-5;
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
    
    %perform simulation
    [t_s,x_s] = ode23(@(t,x) rel(t,x,X_l(:,i),X_r(:,i),M,Tlag_l(i),Tlag_r,xl0,xr0,X_el),[0,0.00015],x0);
    %statistical analysis
    V_res(i) = x_s(end,8);
end
close(h)
figure
h = histogram(V_res,'Normalization','Probability','NumBins',500);
H = 0;
while H<0.997 && i<n
    H = H + h.Values(i);
    i = i + 1;
end
V_res_3s = h.BinEdges(i-1);
h = histogram(V_res,'Normalization','Probability','NumBins',100);