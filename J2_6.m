%% J2~6
function [a_J2_6] = J2_6(reci)
global re mue
% Zonal harmonics
j2 = 0.00108263;
j3 = -2.33936e-3*j2;
j4 = -1.49601e-3*j2;
j5 = -0.20995e-3*j2;
j6 = 0.49941e-3*j2;

r = norm(reci); % Position magnitude [km]

r_I = reci(1);
r_J = reci(2);
r_K = reci(3);

% Perturbation calculation
% J2
factor_J2 = -(3*j2*mue*re^2)/(2*r^5);
a_I_J2 = factor_J2*r_I*(1 - ((5*(r_K^2))/r^2));
a_J_J2 = factor_J2*r_J*(1 - ((5*(r_K^2))/r^2));
a_K_J2 = factor_J2*r_K*(3 - ((5*(r_K^2))/r^2));

% J3
factor_J3 = -(5*j3*mue*re^3)/(2*r^7);
a_I_J3 = factor_J3*r_I*((3*r_K) - ((7*(r_K^3))/r^2));
a_J_J3 = factor_J3*r_J*((3*r_K) - ((7*(r_K^3))/r^2));
a_K_J3 = factor_J3*(6*r_K^2 - ((7*(r_K^4))/r^2) - ((3/5)*r^2));

% J4
factor_J4 = (15*j4*mue*re^4)/(8*r^7);
a_I_J4 = factor_J4*r_I*(1 - ((14*(r_K^2))/r^2) + ((21*(r_K^4))/r^4));
a_J_J4 = factor_J4*r_J*(1 - ((14*(r_K^2))/r^2) + ((21*(r_K^4))/r^4));
a_K_J4 = factor_J4*r_K*(5 - ((70*(r_K^2))/(3*r^2)) + ((21*(r_K^4))/r^4));

% J5
factor_J5 = (3*j5*mue*re^5*r_K)/(8*r^9);
a_I_J5 = factor_J5*r_I*(35 - ((210*(r_K^2))/r^2) + ((231*(r_K^4))/r^4));
a_J_J5 = factor_J5*r_J*(35 - ((210*(r_K^2))/r^2) + ((231*(r_K^4))/r^4));
a_K_J5 = factor_J5*r_K*(105 - ((315*(r_K^2))/r^2) + ((231*(r_K^4))/r^4))...
    - (15*j5*mue*re^5)/(8*r^7);

% J6
factor_J6 = -(j6*mue*re^6)/(16*r^9);
a_I_J6 = factor_J6*r_I*(35 - ((945*(r_K^2))/r^2) + ((3465*(r_K^4))/r^4)...
    - ((3003*(r_K^6))/r^6));
a_J_J6 = factor_J6*r_J*(35 - ((945*(r_K^2))/r^2) + ((3465*(r_K^4))/r^4)...
    - ((3003*(r_K^6))/r^6));
a_K_J6 = factor_J6*r_K*(245 - ((2205*(r_K^2))/(3*r^2))...
    + ((4851*(r_K^4))/r^4) - ((3003*(r_K^6))/r^6));

% Perturbation
a_I = a_I_J2 + a_I_J3 + a_I_J4 + a_I_J5 + a_I_J6;
a_J = a_J_J2 + a_J_J3 + a_J_J4 + a_J_J5 + a_J_J6;
a_K = a_K_J2 + a_K_J3 + a_K_J4 + a_K_J5 + a_K_J6;

% ECEF to ECI
a_J2_6 = [a_I;a_J;a_K];

end

