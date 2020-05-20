%% COE to RV
% Radian
% Ryo Takatori

function [r_peri,v_peri,r_eci,v_eci] = coe2rv_rad(mu,ecc,h,inc,RAAN,omega,theta)
    % Radians
    % State vectors in perifocal
    r_peri = ((h^2)/mu)*(1/(1+ecc*cos(theta)))*[cos(theta) sin(theta) 0];
    v_peri = (mu/h)*[-sin(theta) ecc+cos(theta) 0];
    
    % Perifocal to ECI
    [r_eci,v_eci ] = peri2eci_rad(r_peri,v_peri,omega,RAAN,inc );

end
