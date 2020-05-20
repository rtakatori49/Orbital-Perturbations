%% Variation of Parameter
function [dstatedt] = VoP(t, state, date, date_1, jd, A, m)
global mue solar drag sunbody moonbody jupiterbody J26
% Get COE
h = state(1);
ecc = state(2);
theta = state(3);
inc = state(4);
RAAN = state(5);
omega = state(6);

% Getting RV vectors
[r_peri,v_peri,reci,veci] = coe2rv_rad(mue,ecc,h,inc,RAAN,omega,theta);

% Compute perturbing force
r = norm(reci);
N_hat = cross(reci,veci)/norm(cross(reci,veci));
R_hat = reci/r;
T_hat = cross(N_hat,reci)/norm(cross(N_hat,reci));

% Perturbing forces in same frame
if drag == 1
    [a_drag] = drag_exp(veci, reci, A, m);
else
    a_drag = [0 0 0]';
end
if J26 == 1
    [a_J2_6] = J2_6(reci);
else
    a_J2_6 = [0 0 0]';
end
if sunbody == 1
[a_sun_grav] = sun_N_body(t, jd, reci);
else
    a_sun_grav = [0 0 0]';
end
if moonbody == 1
[a_lunar_grav] = lunar_N_body(t, jd, reci);
else
    a_lunar_grav = [0 0 0]';
end
if jupiterbody == 1
[a_jupiter_grav] = jupiter_N_body(t, jd, reci);
else
    a_jupiter_grav = [0 0 0]';
end
if solar == 1
    [a_srp] = srp(t, date_1, jd, reci, A, m);
else
    a_srp = [0 0 0]';
end
a_pert = a_sun_grav + a_lunar_grav + a_jupiter_grav + a_drag + a_J2_6 + a_srp;
% a_pert = a_sun_grav + a_drag + a_srp;
R = dot(a_pert,R_hat);
T = dot(a_pert,T_hat);
N = dot(a_pert,N_hat);

% d(element)/dt calculation
% Angular momemtum
dhdt = r*T;
% Eccentricity
dedt = (h/mue)*sin(theta)*R+(1/(mue*h))*((h^2+(mue*r))*cos(theta)+mue*ecc*r)*T;
% True anomaly
dthetadt_two_body = (h/r^2);
dthetadt_pert = (1/(ecc*h))*((((h^2)*R)/mue)*cos(theta)-(((h^2)/mue)+r)*T*sin(theta));
dthetadt_all = dthetadt_two_body+dthetadt_pert;
% Inclination
u = omega+theta;
dincdt = (r/h)*N*cos(u);
% Right ascention of ascending node
dRAANdt = ((r*sin(u))/(h*sin(inc)))*N;
% Argument of periapse
domegadt = ((-r*sin(u))/(h*tan(inc)))*N-dthetadt_pert;

% Return states
dstatedt = [dhdt;dedt;dthetadt_all;dincdt;dRAANdt;domegadt];
end

