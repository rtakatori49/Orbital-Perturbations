%% N-body (Jupiter)
function [a_jupiter_grav] = jupiter_N_body(t, date, r_earth)
global muj

% Jupiter position
r_jupiter = planetEphemeris(date + t/(24*60*60), 'Earth', 'Jupiter')';
R_jupiter = norm(r_jupiter);

% Earth to Jupiter
r_earth_jupiter = r_jupiter - r_earth;
R_earth_jupiter = norm(r_earth_jupiter);

% Binomial expansion
q = dot(r_earth, ((2*r_jupiter) - r_earth))/R_jupiter^2;
F = ((q^2 - 3*q + 3)/(1 + (1 - q)^(3/2)))*q;

% Perturbation
a_jupiter_grav = (muj/R_earth_jupiter^3)*((F*r_jupiter) - r_earth);
end
