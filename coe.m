%% Aero 215 HW3
% Ryo Takatori
% 10/27/17
% Introduction to Aerospace Design

function [a, ecc, h, inc, RAAN, omega, theta ] = coe(r, v)
global mue
%% Given inertial and velocity vectors
R = norm(r); % Inertial postion vector [km]
V = norm(v); % Velocity vector [km/s]

%% Gravitational Constant and Semi-major axis (a)
epsilon = ((V.^2)/2)-(mue/R); % Equation for specific mechanical energy [MJ/kg]
a = -mue/(2*epsilon); % Equation for semi-major axis

%% Eccentricity (e)
E = (1/mue)*((((V^2)-(mue/R))*r)-(dot(r,v)*v)); % Equation for eccentricity
ecc = norm(E); % Magnitude of the eccentricity

%% Inclination (inc)
H = cross(r,v); % Cross product of r and v
h = norm(H); % Magnitude of h
inc = acosd(H(3)/h); % Equation for inclination angle [degree]

%% RAAN
k = [0,0,1]; % K hat vector
I = [1,0,0]; % I hat vector
n = cross(k,H); % Cross product of k hat and h
N = norm(n); % Magnitude of n
RAAN = acosd(dot(I,n)/N); % Equation for RAAN [degree]
if n(2) < 0 % if statement for checking quadrant ambiguity
    RAAN = 360-RAAN;
end

%% Argument of perigee (omega)
omega = acosd(dot(n,E)/(N*ecc)); % Equation for the argument of perigee [degree]
if E(3) < 0 % if statement for checking quadrant ambiguity
    omega = 360-omega;
end
    
%% True anomally (theta)
theta = acosd(dot(E,r)/(ecc*R)); % Equation for true anomally [degree]
if dot(r,v) < 0 % if statement for checking quadrant ambiguity
    theta = 360-theta;
end