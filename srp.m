function [a_srp] = srp(t, date_1, jd, r, A, m)
global C_R j4_doy
% Time adjustment
now = date_1 + seconds(t); % Current time
doy = day(now, 'dayofyear'); % Current time in date of year
dsa = doy - j4_doy;% Days since aphelion (July 4)
if dsa < 0
    dsa = dsa + 365;
end
D_a = 2*pi*dsa;
S = 1358/(1.004 + 0.0534*cos(D_a)); % Solar flux calculation [W/m^2]
c = 2.998e8; % Speed of light [m/s]
P_SR = S/c;
[nu] = shadow(jd, t, r);
a_srp = -((P_SR*C_R*A)/m)*(r/norm(r))*nu; % Perturbation
end

