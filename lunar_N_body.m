%% N-body (Moon)
function [a_lunar_grav] = lunar_N_body(t, date, r_earth)
global mum
% Moon position
r_moon = lunar_position(date + t/(24*60*60))';
R_moon = norm(r_moon);

% Earth to Moon
r_earth_moon = r_moon - r_earth;
R_earth_moon = norm(r_earth_moon);

% Binomial expansion
q = dot(r_earth, ((2*r_moon) - r_earth))/R_moon^2;
F = ((q^2 - 3*q + 3)/(1 + (1 - q)^(3/2)))*q;

% Perturbation
a_lunar_grav = (mum/R_earth_moon^3)*((F*r_moon) - r_earth);
end
