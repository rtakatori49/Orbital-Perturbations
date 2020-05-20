%% N-body (Sun)
function [a_sun_grav] = sun_N_body(t, date, r_earth)
global mus
% Sun position
[r_sun,rtasc,decl] = sun(date + t/(24*60*60));
r_sun = r_sun'*1.496e+8;
R_sun = norm(r_sun);

% Earth to Sun
r_earth_sun = r_sun - r_earth;
R_earth_sun = norm(r_earth_sun);

% Binomial expansion
q = dot(r_earth, ((2*r_sun) - r_earth))/R_sun^2;
F = ((q^2 - 3*q + 3)/(1 + (1 - q)^(3/2)))*q;

% Perturbation
a_sun_grav = (mus/R_earth_sun^3)*((F*r_sun) - r_earth);
end

