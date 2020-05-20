%% Drag
% Exponential model using US standard
function [a_drag] = drag_exp(v, r, A, m)
global ang_vel re C_D
% Relative velocity
v_rel = v - cross(ang_vel, r);
V_rel = norm(v_rel);
% Altitude
R = norm(r);
alt = R - re;
% Density
[rho] = atmosphere(alt);
% Drag acceleration
a_drag = -(0.5*C_D*A*rho*V_rel*v_rel*1000)/m;
end