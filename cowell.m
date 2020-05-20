%% Two Body Problem
function [dstatedt] = cowell(t, state, date, date_1, jd, A, m)
global mue solar drag sunbody moonbody jupiterbody J26
%% Position [km]
% Component
x = state(1);
y = state(2);
z = state(3);
% Vector
r = [x y z]';
% Magnitude
R = norm(r);

%% Velocity [km/s]
% Component
dx = state(4);
dy = state(5);
dz = state(6);
% Vector
v = [dx dy dz]';
% Magnitude
V = norm(v);

%% Acceleration [km/s^2]
% Component
ddx = -mue*x/R^3;
ddy = -mue*y/R^3;
ddz = -mue*z/R^3;
% Vector
a = [ddx ddy ddz]';

% Perturbing forces in same frame
if drag == 1
    [a_drag] = drag_exp(v, r, A, m);
else
    a_drag = [0 0 0]';
end
if J26 == 1
    [a_J2_6] = J2_6(r);
else
    a_J2_6 = [0 0 0]';
end
if sunbody == 1
[a_sun_grav] = sun_N_body(t, jd, r);
else
    a_sun_grav = [0 0 0]';
end
if moonbody == 1
[a_lunar_grav] = lunar_N_body(t, jd, r);
else
    a_lunar_grav = [0 0 0]';
end
if jupiterbody == 1
[a_jupiter_grav] = jupiter_N_body(t, jd, r);
else
    a_jupiter_grav = [0 0 0]';
end
if solar == 1
    [a_srp] = srp(t, date_1, jd, r, A, m);
else
    a_srp = [0 0 0]';
end
a_pert = a_sun_grav + a_lunar_grav + a_jupiter_grav + a_drag + a_J2_6 + a_srp;

a = a + a_pert;
dstatedt = [v;a]; % Return state vector for next step
end

