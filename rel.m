% ODE function for tip release

function dx = rel(t,x,par_l,par_r,M,I,Tlag_l,Tlag_r,mis_l,mis_r,xl0,xr0,X_el,f_n)

% state definition
xl = x(1:3);
xr = x(4:6);
x_tm = x(7:8);
x_ang = x(9);
noise = 0;
%% system of equations
% adhesion force model
e_l = 1e6*(xl0(2) - xl(2) + x_tm(1));         %left tip elongation
e_r = 1e6*(xr0(2) - xr(2) - x_tm(1));         %right tip elongation
% input voltage definition
ul = 120 - 60*(1+sign(t-Tlag_l));
ur = 120 - 60*(1+sign(t-Tlag_r));

% force input definition
%left hand tip
if t-Tlag_l<1e-4
    if e_l<3 && e_l>1e-7
        F_l = -X_el(1)*e_l*exp(-X_el(2)*e_l^X_el(3));
    elseif e_l<1e-7
        F_l = min([-3e6*e_l+0.3,0.3]);
    else
        F_l = 0;
    end
else
    F_l = 0;
end
%add noise force
if t > Tlag_l
    noise = ceil(length(f_n)*t/0.00015);
end

% right hand tip
if t - Tlag_r < 1e-4
    if e_r < 3 && e_r > 1e-7
        F_r = -X_el(1)*e_r*exp(-X_el(2)*e_r^X_el(3));
    elseif e_r < 1e-7
        F_r = min([-3e6*e_l+0.3,0.3]);
    else
        F_r = 0;
    end
else
     F_r = 0;
end

% tips equations
[Al, Bl,~,~,~,~] = sys2(par_l, 0);
[Ar, Br,~,~,~,~] = sys2(par_r, 0);
dxl = Al*xl + Bl*ul - [0;0;-1]*(F_l*noise);
dxr = Ar*xr + Br*ur + [0;0;-1]*F_r;

% test mass equations
dx_tm = [0, 1; 0, 0]*x_tm + [0; 1/M]*(F_l - F_r);
dx_ang = (mis_r*F_r - mis_l*F_l)/I;
dx = [dxl;dxr;dx_tm;dx_ang];