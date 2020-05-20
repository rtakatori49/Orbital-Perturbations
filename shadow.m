function [nu] = shadow(date, t, r_sc)
re = 6378; % Earth radius [km]
[r_sun, rtasc, rdecl] = sun(date + t/(24*60*60)); % Vallado sun calculation
% [(sun vector) [km],(RAAN) [rad],(decl) [rad]]
R_earth = norm(r_sc); % Spacecraft position magnitude [km]
R_sun = norm(r_sun)*1.496e+8; % Sun position magnitude [km]
% Angles from center of Earth required to calculate shadow
theta = acosd(dot(r_sun*1.496e+8,r_sc)/(R_sun*R_earth));
theta_1 = acosd(re/R_earth);
theta_2 = acosd(re/R_sun);
% Checking shadow
if theta_1+theta_2 <= theta
    nu = 0;
else
    nu = 1;
end
end