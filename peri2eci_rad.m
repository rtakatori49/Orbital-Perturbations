function [reci,veci] = peri2eci_rad(rperi,vperi,omega,OMEGA,inc)
% Rotate to ECI
R3omega = [cos(omega) sin(omega) 0;
    -sin(omega) cos(omega) 0;
    0 0 1];

R1inc = [1 0 0;
    0 cos(inc) sin(inc);
    0 -sin(inc) cos(inc)];

R3OMEGA = [cos(OMEGA) sin(OMEGA) 0;
    -sin(OMEGA) cos(OMEGA) 0;
    0 0 1];

Q = (R3omega*R1inc*R3OMEGA)';

reci = Q*rperi';
veci = Q*vperi';
end

