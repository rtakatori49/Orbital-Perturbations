%% COE to RV
% Degree
% Ryo Takatori
function [r_peri,v_peri,r_eci,v_eci] = coe2rv_deg(mu,ecc,h,inc,RAAN,omega,theta)
    % State vectors in perifocal
    r_peri = ((h^2)/mu)*(1/(1+ecc*cosd(theta)))*[cosd(theta) sind(theta) 0];
    v_peri = (mu/h)*[-sind(theta) ecc+cosd(theta) 0];
    
    % Perifocal to ECI
    [r_eci,v_eci ] = peri2eci_deg(r_peri,v_peri,omega,RAAN,inc,unit);
end