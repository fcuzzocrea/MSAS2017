function dx = rel(t,x,par_l,par_r,M,I,Tlag_l,Tlag_r,mis_l,mis_r,xl0,xr0,X_el)
% ODE function for tip release

% State definition
xl = x(1:3);                % Displacement of the left tip
xr = x(4:6);                % Displacement of the right tip 
x_tm = x(7:8);              % Displacement of the Test Mass
x_ang = x(9);               % Angular displacemente of the Test Mass

% Adhesion force model
e_l = 1e6*(xl0(2) - xl(2) + x_tm(1));         % Left tip elongation
e_r = 1e6*(xr0(2) - xr(2) - x_tm(1));         % Right tip elongation

% Input voltage definition
ul = 120 - 60*(1+sign(t-Tlag_l));             % Volt 
ur = 120 - 60*(1+sign(t-Tlag_r));             % Volt

% Force input definition

% Left hand tip
if t-Tlag_l < 1e-4                                      % ??????????    
    if e_l < 3 && e_l > 1e-7                            % ?????????? Ci vuole un 10^-6 ?
        F_l = -X_el(1)*e_l*exp(-X_el(2)*e_l^X_el(3));
    elseif e_l < 1e-7                                   % ??????????
        F_l = min([-3e6*e_l+0.3,0.3]);
    else
        F_l = 0;
    end
else
    F_l = 0;
end

% Right hand tip
if t-Tlag_r < 1e-4
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

% Tips equations
[Al, Bl,~,~,~,~] = sys2(par_l, 0);
[Ar, Br,~,~,~,~] = sys2(par_r, 0);
dxl = Al*xl + Bl*ul - [0;0;-1]*F_l;
dxr = Ar*xr + Br*ur + [0;0;-1]*F_r;

% Test mass equations     
dx_tm = [0, 1; 0, 0]*x_tm + [0; 1/M]*(F_l - F_r);
dx_ang = (mis_r*F_r - mis_l*F_l)/I;
dx = [dxl;dxr;dx_tm;dx_ang];

end