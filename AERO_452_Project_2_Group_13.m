%% Project 2
% AERO 452
% Ryo Takatori
% Group 13
% 11/21/2019

clc, clear all, close all
global mue mum mus muj me mm re ang_vel C_D C_R j4_doy solar drag ...
    sunbody moonbody jupiterbody J26 duration
% Commonly used constants
mue = 398600; % Earth gravitational constant [km^3/s^2]
mum = 4902; % Moon gravitational constant [km^3/s^2]
mus = 132.712e9; % Sun gravitational constant [km^3/s^2]
muj = 126686000; % Jupiter gravitational constant [km^3/s^2]
me = 5.974*10^24; % Mass of Earth [kg]
mm = 73.48*10^21; % Mass of Moon [kg]
re = 6378; % Radius of Earth [km]
ang_vel  = [0;0;72.9211e-6]; % Angular velocity of Earth [rad/s]
C_D = 2.2; % Coefficient of drag
C_R = 1.2; % Coefficient of
j4_doy = day(datetime([2019 7 4 0 0 0]), 'dayofyear'); % July 4 in date of year
options = odeset('RelTol',1e-8,'AbsTol',1e-8); % Ode settings
% Short Term
t_short = 30*24*60*60;
tspan_short = [0 t_short]; % Time span [s]
% Long Term
t_long = 365*24*60*60;
tspan_long = [0 t_long]; % Time span [s]

%% Object 1 (LEMUR-2 JOEL)
disp('Object 1 (LEMUR-2 JOEL)')
solar = 1;
drag = 1;
sunbody = 1;
moonbody = 0;
jupiterbody = 0;
J26 = 0;
m_1 = 4; % Mass [kg]
A_1 = 0.1*0.1*3; % Area [m^2] (3U CubeSat 10 [cm] x 10 [cm] x 30 [cm])
[Me0_1, n0_1, ecc0_1, inc0_1, RAAN0_1, omega0_1, epoch0_1, tle0_1, a0_1,...
    E0_1, theta0_1, h0_1, T0_1, reci0_1, veci0_1] = TLE_Reader('LEMUR-2 JOEL.txt');
rp0_1 = a0_1 * (1 - ecc0_1); % Perigee [km]
ra0_1 = (2 * a0_1) - rp0_1; % Apogee [km]
zp0_1 = rp0_1 - re; % Perigee altitude [km]
za0_1 = ra0_1 - re; % Apogee altitude [km]
T_1 = 2*pi*sqrt((a0_1^3)/mue);
coe0_1 = [zp0_1, za0_1, a0_1, ecc0_1, h0_1, inc0_1, RAAN0_1, omega0_1, theta0_1];
date0_1 = datevec(datenum(datevec(seconds(epoch0_1)) + [2019 0 0 0 0 0]));

% In radians
RAAN0_1_rad = deg2rad(RAAN0_1); % RAAN [rad]
inc0_1_rad = deg2rad(inc0_1); % Inclination [rad]
omega0_1_rad = deg2rad(omega0_1); % Arguement of perigee [rad]
theta0_1_rad = deg2rad(theta0_1); % True Anomaly [rad]
state_1 = [h0_1,ecc0_1,theta0_1_rad,inc0_1_rad,RAAN0_1_rad,...
omega0_1_rad]'; % State vector

% Checking Orbit without perturbation short
dt_1 = T_1/2;
tspan_np_s_1 = [0 t_short + dt_1];
state_np_1 = [reci0_1; veci0_1]; % State vector

[tnew_np_s_1, statenew_s_np_1] = ode45(@two_body, tspan_np_s_1, state_np_1, options); % ode45
rnew_np_s_1 = statenew_s_np_1(end,1:3)';
vnew_np_s_1 = statenew_s_np_1(end,4:6)';
tspan_np_L_1 = [0 t_long + dt_1];
[tnew_np_L_1, statenew_L_np_1] = ode45(@two_body, tspan_np_L_1, state_np_1, options); % ode45
rnew_np_L_1 = statenew_L_np_1(end,1:3)';
vnew_np_L_1 = statenew_L_np_1(end,4:6)';

% Plot of orbit
figure
plot3(statenew_s_np_1(:,1), statenew_s_np_1(:,2), statenew_s_np_1(:,3),'Color','r')
hold on
plot3(reci0_1(1), reci0_1(2), reci0_1(3), 'o', 'MarkerSize', 5, 'Color', 'b')
[x,y,z] = sphere;
globe = surf(x*re,y*re,z*re);
cdata = imread('earth_base_high_1.jpg');
sph1 = findobj('Type', 'surface');
set(globe,'FaceColor','texturemap','CData',cdata,'EdgeColor','none');
title('No Perturbation LEMUR-2 JOEL Orbit')
xlabel('x [km]')
ylabel('y [km]')
zlabel('z [km]')
legend('Orbit', 'Starting Point')
grid on
axis equal

% Short term
disp('Short Term')
profile on
[tnew_s_1, rnew_s_1, vnew_s_1, statenew_s_1, diff_COE_s_1, timer_s_1] = proj2_solver_VoP(@VoP, tspan_short,...
    state_1, options, date0_1, A_1, m_1, coe0_1, 'Short Term LEMUR-2 JOEL');
% [tnew_s_1,statenew_s_1,diff_COE_s_1,timer_s_1] = proj2_solver_ode(@cowell ,tspan_short,...
%     state_np_1,options,date0_1,A_1,m_1,coe0_1,'Short Term LEMUR-2 JOEL');
% rnew_s_1 = statenew_s_1(:,1:3);
% vnew_s_1 = statenew_s_1(:,4:6);
fprintf('Time it took to propagate Variation of Parameter: %f [s]\n\n',...
timer_s_1)
profile viewer
[v1_s_1, v2_s_1] = lambert(rnew_s_1(end,:)', rnew_np_s_1, dt_1);
delta_v_s_1_1 = vnew_s_1(end,:)' - v1_s_1;
delta_v_s_2_1 = v2_s_1 - vnew_np_s_1;
delta_v_s_1 = norm(delta_v_s_1_1) +  norm(delta_v_s_2_1);
fprintf('Delta-V to Correct %f [km/s]\n\n',delta_v_s_1)

% Long term
disp('Long Term')
profile on
[tnew_L_1, rnew_L_1, vnew_L_1, statenew_L_1, diff_COE_L_1, timer_L_1] = proj2_solver_VoP(@VoP, tspan_long,...
    state_1, options, date0_1, A_1, m_1, coe0_1, 'Long Term LEMUR-2 JOEL');
% [tnew_L_1,statenew_L_1,diff_COE_L_1,timer_L_1] = proj2_solver_ode(@cowell ,tspan_long,...
%     state_np_1,options,date0_1,A_1,m_1,coe0_1,'Long Term LEMUR-2 JOEL');
% rnew_L_1 = statenew_L_1(:,1:3);
% vnew_L_1 = statenew_L_1(:,4:6);
fprintf('Time it took to propagate Variation of Parameter: %f [s]\n\n',...
timer_L_1)
profile viewer
[v1_L_1, v2_L_1] = lambert(rnew_L_1(end,:)', rnew_np_L_1, dt_1);
delta_v_L_1_1 = vnew_L_1(end,:)' - v1_L_1;
delta_v_L_2_1 = v2_L_1 - vnew_np_L_1;
delta_v_L_1 = norm(delta_v_L_1_1) +  norm(delta_v_L_2_1);
fprintf('Delta-V to Correct %f [km/s]\n\n',delta_v_L_1)

%% Object 2 (BSAT-3C)
disp('Object 2 (BSAT-3C)')
solar = 0;
drag = 0;
sunbody = 0;
moonbody = 0;
jupiterbody = 0;
J26 = 0;
m_2 = 2906; % Mass [kg]
A_2 = 5.3*2.0; % Area [m^2] (dimensions: 5.3 [m] x 2.0 [m] x 1.9 [m])
[Me0_2, n0_2, ecc0_2, inc0_2, RAAN0_2, omega0_2, epoch0_2, tle0_2, a0_2,...
    E0_2, theta0_2, h0_2, T0_2, reci0_2, veci0_2] = TLE_Reader('BSAT-3C.txt');
rp0_2 = a0_2 * (1 - ecc0_2); % Perigee [km]
ra0_2 = (2 * a0_2) - rp0_2; % Apogee [km]
zp0_2 = rp0_2 - re; % Perigee altitude [km]
za0_2 = ra0_2 - re; % Apogee altitude [km]
T_2 = 2*pi*sqrt((a0_2^3)/mue);
coe0_2 = [zp0_2, za0_2, a0_2, ecc0_2, h0_2, inc0_2, RAAN0_2, omega0_2, theta0_2];
date0_2 = datevec(datenum(datevec(seconds(epoch0_2)) + [2019 0 0 0 0 0]));

% In radians
RAAN0_2_rad = deg2rad(RAAN0_2); % RAAN [rad]
inc0_2_rad = deg2rad(inc0_2); % Inclination [rad]
omega0_2_rad = deg2rad(omega0_2); % Arguement of perigee [rad]
theta0_2_rad = deg2rad(theta0_2); % True Anomaly [rad]
state_2 = [h0_2,ecc0_2,theta0_2_rad,inc0_2_rad,RAAN0_2_rad,...
omega0_2_rad]'; % State vector

% Checking Orbit without perturbation
% Short
dt_2 = T_2/2;
tspan_np_s_2 = [0 t_short + dt_2];
state_np_2 = [reci0_2; veci0_2]; % State vector
[tnew_np_s_2, statenew_np_s_2] = ode45(@two_body, tspan_short, state_np_2, options); % ode45
rnew_np_s_2 = statenew_np_s_2(end,1:3)';
vnew_np_s_2 = statenew_np_s_2(end,4:6)';
% Long
tspan_np_L_2 = [0 t_long + dt_2];
[tnew_np_L_2, statenew_np_L_2] = ode45(@two_body, tspan_short, state_np_2, options); % ode45
rnew_np_L_2 = statenew_np_L_2(end,1:3)';
vnew_np_L_2 = statenew_np_L_2(end,4:6)';

duration = 0;
% [tnew_c_s, rnew_c_s, vnew_c_s] = correction(state_np_2, state_2, options, date0_2, A_2, m_2, coe0_2, dt_2);

% Plot of orbit
figure
plot3(statenew_np_s_2(:,1), statenew_np_s_2(:,2), statenew_np_s_2(:,3),'Color','r')
hold on
plot3(reci0_2(1), reci0_2(2), reci0_2(3), 'o', 'MarkerSize', 5, 'Color', 'b')
[x,y,z] = sphere;
globe = surf(x*re,y*re,z*re);
cdata = imread('earth_base_high_1.jpg');
sph1 = findobj('Type', 'surface');
set(globe,'FaceColor','texturemap','CData',cdata,'EdgeColor','none');
title('No Perturbation BSAT-3C Orbit')
xlabel('x [km]')
ylabel('y [km]')
zlabel('z [km]')
legend('Orbit', 'Starting Point')
grid on
axis equal

% Short term
disp('Short Term')
profile on
[tnew_s_2, rnew_s_2, vnew_s_2, statenew_s_2, diff_COE_s_2, timer_s_2] = proj2_solver_VoP(@VoP, tspan_short,...
    state_2, options, date0_2, A_2, m_2, coe0_2, 'Short Term BSAT-3C');
[tnew_s_2,statenew_s_2,diff_COE_s_2,timer_s_2] = proj2_solver_ode(@cowell ,tspan_short,...
    state_np_2,options,date0_2,A_2,m_2,coe0_2,'Short Term BSAT-3C');
rnew_s_2 = statenew_s_2(:,1:3);
vnew_s_2 = statenew_s_2(:,4:6);
fprintf('Time it took to propagate Variation of Parameter: %f [s]\n\n',...
timer_s_2)
profile viewer
[v1_s_2, v2_s_2] = lambert(rnew_s_2(end,:)', rnew_np_s_2, dt_2);
delta_v_s_1_2 = vnew_s_2(end,:)' - v1_s_2;
delta_v_s_2_2 = v2_s_2 - vnew_np_s_2;
delta_v_s_2 = norm(delta_v_s_1_2) +  norm(delta_v_s_2_2);
fprintf('Delta-V to Correct %f [km/s]\n\n',delta_v_s_2)

% Long term
disp('Long Term')
profile on
[tnew_L_2, rnew_L_2, vnew_L_2, statenew_L_2, diff_COE_L_2, timer_L_2] = proj2_solver_VoP(@VoP, tspan_long,...
    state_2, options, date0_2, A_2, m_2, coe0_2, 'Long Term BSAT-3C');
[tnew_L_2,statenew_L_2,diff_COE_L_2,timer_L_2] = proj2_solver_ode(@cowell ,tspan_long,...
    state_np_2,options,date0_2,A_2,m_2,coe0_2,'Long Term BSAT-3C');
rnew_L_2 = statenew_L_2(:,1:3);
vnew_L_2 = statenew_L_2(:,4:6);
fprintf('Time it took to propagate Variation of Parameter: %f [s]\n\n',...
timer_L_2)
profile viewer
[v1_L_2, v2_L_2] = lambert(rnew_L_2(end,:)', rnew_np_L_2, dt_2);
delta_v_L_1_2 = vnew_L_2(end,:)' - v1_L_2;
delta_v_L_2_2 = v2_L_2 - vnew_np_L_2;
delta_v_L_2 = norm(delta_v_L_1_2) +  norm(delta_v_L_2_2);
fprintf('Delta-V to Correct %f [km/s]\n\n',delta_v_L_2)

% %% Object 2 (BSAT-3C) with Jupiter
% disp('Object 2 (BSAT-3C)')
% solar = 1;
% drag = 0;
% sunbody = 1;
% moonbody = 1;
% jupiterbody = 1;
% J26 = 0;
% m_2 = 2906; % Mass [kg]
% A_2 = 5.3*2.0; % Area [m^2] (dimensions: 5.3 [m] x 2.0 [m] x 1.9 [m])
% [Me0_2, n0_2, ecc0_2, inc0_2, RAAN0_2, omega0_2, epoch0_2, tle0_2, a0_2,...
%     E0_2, theta0_2, h0_2, T0_2, reci0_2, veci0_2] = TLE_Reader('BSAT-3C.txt');
% rp0_2 = a0_2 * (1 - ecc0_2); % Perigee [km]
% ra0_2 = (2 * a0_2) - rp0_2; % Apogee [km]
% zp0_2 = rp0_2 - re; % Perigee altitude [km]
% za0_2 = ra0_2 - re; % Apogee altitude [km]
% T_2 = 2*pi*sqrt((a0_2^3)/mue);
% coe0_2 = [zp0_2, za0_2, a0_2, ecc0_2, h0_2, inc0_2, RAAN0_2, omega0_2, theta0_2];
% date0_2 = datevec(datenum(datevec(seconds(epoch0_2)) + [2019 0 0 0 0 0]));
% 
% % In radians
% RAAN0_2_rad = deg2rad(RAAN0_2); % RAAN [rad]
% inc0_2_rad = deg2rad(inc0_2); % Inclination [rad]
% omega0_2_rad = deg2rad(omega0_2); % Arguement of perigee [rad]
% theta0_2_rad = deg2rad(theta0_2); % True Anomaly [rad]
% state_2 = [h0_2,ecc0_2,theta0_2_rad,inc0_2_rad,RAAN0_2_rad,...
% omega0_2_rad]'; % State vector
% 
% % Checking Orbit without perturbation
% % Short
% dt_2 = T_2/2;
% tspan_np_s_2 = [0 t_short + dt_2];
% state_np_2 = [reci0_2; veci0_2]; % State vector
% [tnew_np_s_2, statenew_np_s_2] = ode45(@two_body, tspan_short, state_np_2, options); % ode45
% rnew_np_s_2 = statenew_np_s_2(end,1:3)';
% vnew_np_s_2 = statenew_np_s_2(end,4:6)';
% % Long
% tspan_np_L_2 = [0 t_long + dt_2];
% [tnew_np_L_2, statenew_np_L_2] = ode45(@two_body, tspan_short, state_np_2, options); % ode45
% rnew_np_L_2 = statenew_np_L_2(end,1:3)';
% vnew_np_L_2 = statenew_np_L_2(end,4:6)';

% 
% % Plot of orbit
% figure
% plot3(statenew_np_s_2(:,1), statenew_np_s_2(:,2), statenew_np_s_2(:,3),'Color','r')
% hold on
% plot3(reci0_2(1), reci0_2(2), reci0_2(3), 'o', 'MarkerSize', 5, 'Color', 'b')
% [x,y,z] = sphere;
% globe = surf(x*re,y*re,z*re);
% cdata = imread('earth_base_high_1.jpg');
% sph1 = findobj('Type', 'surface');
% set(globe,'FaceColor','texturemap','CData',cdata,'EdgeColor','none');
% title('No Perturbation BSAT-3C Orbit')
% xlabel('x [km]')
% ylabel('y [km]')
% zlabel('z [km]')
% legend('Orbit', 'Starting Point')
% grid on
% axis equal
% 
% % Short term
% disp('Short Term')
% profile on
% [tnew_s_2, rnew_s_2, vnew_s_2, statenew_s_2, diff_COE_s_2, timer_s_2] = proj2_solver_VoP(@VoP, tspan_short,...
%     state_2, options, date0_2, A_2, m_2, coe0_2, 'Short Term BSAT-3C with Jupiter');
% fprintf('Time it took to propagate Variation of Parameter: %f [s]\n\n',...
% timer_s_2)
% profile viewer
% [v1_s_2, v2_s_2] = lambert(rnew_s_2(end,:)', rnew_np_s_2, dt_2);
% delta_v_s_1_2 = vnew_s_2(end,:)' - v1_s_2;
% delta_v_s_2_2 = v2_s_2 - vnew_np_s_2;
% delta_v_s_2 = norm(delta_v_s_1_2) +  norm(delta_v_s_2_2);
% fprintf('Delta-V to Correct %f [km/s]\n\n',delta_v_s_2)
% 
% % Long term
% disp('Long Term')
% profile on
% [tnew_L_2, rnew_L_2, vnew_L_2, statenew_L_2, diff_COE_L_2, timer_L_2] = proj2_solver_VoP(@VoP, tspan_long,...
%     state_2, options, date0_2, A_2, m_2, coe0_2, 'Long Term BSAT-3C with Jupiter');
% fprintf('Time it took to propagate Variation of Parameter: %f [s]\n\n',...
% timer_L_2)
% profile viewer
% [v1_L_2, v2_L_2] = lambert(rnew_L_2(end,:)', rnew_np_L_2, dt_2);
% delta_v_L_1_2 = vnew_L_2(end,:)' - v1_L_2;
% delta_v_L_2_2 = v2_L_2 - vnew_np_L_2;
% delta_v_L_2 = norm(delta_v_L_1_2) +  norm(delta_v_L_2_2);
% fprintf('Delta-V to Correct %f [km/s]\n\n',delta_v_L_2)

%% Object 3 (MOLNIYA 3-50)
disp('Object 3 (MOLNIYA 3-50)')
solar = 1;
drag = 1;
sunbody = 1;
moonbody = 1;
jupiterbody = 0;
J26 = 1;
m_3 = 1600; % Mass [kg]
A_3 = 4.4*1.4*6; % Area [m^2] (dimensions: 5.3 [m] x 2.0 [m] x 1.9 [m])
[Me0_3, n0_3, ecc0_3, inc0_3, RAAN0_3, omega0_3, epoch0_3, tle0_3, a0_3,...
    E0_3, theta0_3, h0_3, T0_3, reci0_3, veci0_3] = TLE_Reader('MOLNIYA 3-50.txt');
rp0_3 = a0_3 * (1 - ecc0_3); % Perigee [km]
ra0_3 = (2 * a0_3) - rp0_3; % Apogee [km]
zp0_3 = rp0_3 - re; % Perigee altitude [km]
za0_3 = ra0_3 - re; % Apogee altitude [km]
T_3 = 2*pi*sqrt((a0_3^3)/mue);
coe0_3 = [zp0_3, za0_3, a0_3, ecc0_3, h0_3, inc0_3, RAAN0_3, omega0_3, theta0_3];
date0_3 = datevec(datenum(datevec(seconds(epoch0_3)) + [2019 0 0 0 0 0]));

% In radians
RAAN0_3_rad = deg2rad(RAAN0_3); % RAAN [rad]
inc0_3_rad = deg2rad(inc0_3); % Inclination [rad]
omega0_3_rad = deg2rad(omega0_3); % Arguement of perigee [rad]
theta0_3_rad = deg2rad(theta0_3); % True Anomaly [rad]
state_3 = [h0_3,ecc0_3,theta0_3_rad,inc0_3_rad,RAAN0_3_rad,...
omega0_3_rad]'; % State vector

% Checking Orbit without perturbation
% Short
dt_3 = T_3/2;
tspan_np_s_3 = [0 t_short + dt_3];
% tspan_np_s_3 = linspace(0, t_short + dt_3, 40000);
state_np_3 = [reci0_3; veci0_3]; % State vector
[tnew_np_s_3, statenew_np_s_3] = ode45(@two_body, tspan_short, state_np_3, options); % ode45
rnew_np_s_3 = statenew_np_s_3(end,1:3)';
vnew_np_s_3 = statenew_np_s_3(end,4:6)';

% Long
tspan_np_L_3 = [0 t_long + dt_3];
% tspan_np_L_3 = linspace(0, t_long + dt_3, 40000);
[tnew_np_L_3, statenew_np_L_3] = ode45(@two_body, tspan_short, state_np_3, options); % ode45
rnew_np_L_3 = statenew_np_L_3(end,1:3)';
vnew_np_L_3 = statenew_np_L_3(end,4:6)';

% Plot of orbit
figure
plot3(statenew_np_s_3(:,1), statenew_np_s_3(:,2), statenew_np_s_3(:,3),'Color','r')
hold on
plot3(reci0_3(1), reci0_3(2), reci0_3(3), 'o', 'MarkerSize', 5, 'Color', 'b')
[x,y,z] = sphere;
globe = surf(x*re,y*re,z*re);
cdata = imread('earth_base_high_1.jpg');
sph1 = findobj('Type', 'surface');
set(globe,'FaceColor','texturemap','CData',cdata,'EdgeColor','none');
title('No Perturbation MOLNIYA 3-50 Orbit')
xlabel('x [km]')
ylabel('y [km]')
zlabel('z [km]')
legend('Orbit', 'Starting Point')
grid on
axis equal

% Short term
disp('Short Term')
profile on
[tnew_s_3, rnew_s_3, vnew_s_3, statenew_s_3, diff_COE_s_3, timer_s_3] = proj2_solver_VoP(@VoP, tspan_short,...
    state_3, options, date0_3, A_3, m_3, coe0_3, 'Short Term MOLNIYA 3-50');
fprintf('Time it took to propagate Variation of Parameter: %f [s]\n\n',...
timer_s_3)
profile viewer
[v1_s_3, v2_s_3] = lambert(rnew_s_3(end,:)', rnew_np_s_3, dt_3);
delta_v_s_1_3 = vnew_s_3(end,:)' - v1_s_3;
delta_v_s_2_3 = v2_s_3 - vnew_np_s_3;
delta_v_s_3 = norm(delta_v_s_1_3) +  norm(delta_v_s_2_3);
fprintf('Delta-V to Correct %f [km/s]\n\n',delta_v_s_3)

% Long term
disp('Long Term')
profile on
[tnew_L_3, rnew_L_3, vnew_L_3, statenew_L_3, diff_COE_L_3, timer_L_3] = proj2_solver_VoP(@VoP, tspan_long,...
    state_3, options, date0_3, A_3, m_3, coe0_3, 'Long Term MOLNIYA 3-50');
fprintf('Time it took to propagate Variation of Parameter: %f [s]\n\n',...
timer_L_3)
profile viewer
[v1_L_3, v2_L_3] = lambert(rnew_L_3(end,:)', rnew_np_L_3, dt_3);
delta_v_L_1_3 = vnew_L_3(end,:)' - v1_L_3;
delta_v_L_2_3 = v2_L_3 - vnew_np_s_3;
delta_v_L_3 = norm(delta_v_L_1_3) +  norm(delta_v_L_2_3);
fprintf('Delta-V to Correct %f [km/s]\n\n',delta_v_L_3)


% %% Object 3 (RADARSAT 2)
% disp('Object 3 (RADARSAT 2)')
% solar = 1;
% drag = 1;
% sunbody = 1;
% moonbody = 1;
% jupiterbody = 0;
% J26 = 1;
% m_3 = 2200; % Mass [kg]
% A_3 = 4.4*1.4*6; % Area [m^2] (dimensions: 5.3 [m] x 2.0 [m] x 1.9 [m])
% [Me0_3, n0_3, ecc0_3, inc0_3, RAAN0_3, omega0_3, epoch0_3, tle0_3, a0_3,...
%     E0_3, theta0_3, h0_3, T0_3, reci0_3, veci0_3] = TLE_Reader('RADARSAT 2.txt');
% rp0_3 = a0_3 * (1 - ecc0_3); % Perigee [km]
% ra0_3 = (2 * a0_3) - rp0_3; % Apogee [km]
% zp0_3 = rp0_3 - re; % Perigee altitude [km]
% za0_3 = ra0_3 - re; % Apogee altitude [km]
% T_3 = 2*pi*sqrt((a0_3^3)/mue);
% coe0_3 = [zp0_3, za0_3, a0_3, ecc0_3, h0_3, inc0_3, RAAN0_3, omega0_3, theta0_3];
% date0_3 = datevec(datenum(datevec(seconds(epoch0_3)) + [2019 0 0 0 0 0]));
% 
% % In radians
% RAAN0_3_rad = deg2rad(RAAN0_3); % RAAN [rad]
% inc0_3_rad = deg2rad(inc0_3); % Inclination [rad]
% omega0_3_rad = deg2rad(omega0_3); % Arguement of perigee [rad]
% theta0_3_rad = deg2rad(theta0_3); % True Anomaly [rad]
% state_3 = [h0_3,ecc0_3,theta0_3_rad,inc0_3_rad,RAAN0_3_rad,...
% omega0_3_rad]'; % State vector
% 
% % Checking Orbit without perturbation
% % Short
% dt_3 = T_3/2;
% tspan_np_s_3 = [0 t_short + dt_3];
% % tspan_np_s_3 = linspace(0, t_short + dt_3, 40000);
% state_np_3 = [reci0_3; veci0_3]; % State vector
% [tnew_np_s_3, statenew_np_s_3] = ode45(@two_body, tspan_short, state_np_3, options); % ode45
% rnew_np_s_3 = statenew_np_s_3(end,1:3)';
% vnew_np_s_3 = statenew_np_s_3(end,4:6)';
% 
% % Long
% tspan_np_L_3 = [0 t_long + dt_3];
% % tspan_np_L_3 = linspace(0, t_long + dt_3, 40000);
% [tnew_np_L_3, statenew_np_L_3] = ode45(@two_body, tspan_short, state_np_3, options); % ode45
% rnew_np_L_3 = statenew_np_L_3(end,1:3)';
% vnew_np_L_3 = statenew_np_L_3(end,4:6)';
% 
% % Plot of orbit
% figure
% plot3(statenew_np_s_3(:,1), statenew_np_s_3(:,2), statenew_np_s_3(:,3),'Color','r')
% hold on
% plot3(reci0_3(1), reci0_3(2), reci0_3(3), 'o', 'MarkerSize', 5, 'Color', 'b')
% [x,y,z] = sphere;
% globe = surf(x*re,y*re,z*re);
% cdata = imread('earth_base_high_1.jpg');
% sph1 = findobj('Type', 'surface');
% set(globe,'FaceColor','texturemap','CData',cdata,'EdgeColor','none');
% title('No Perturbation RADARSAT 2')
% xlabel('x [km]')
% ylabel('y [km]')
% zlabel('z [km]')
% legend('Orbit', 'Starting Point')
% grid on
% axis equal
% 
% % Short term
% disp('Short Term')
% profile on
% [tnew_s_3, rnew_s_3, vnew_s_3, statenew_s_3, diff_COE_s_3, timer_s_3] = proj2_solver_VoP(@VoP, tspan_short,...
%     state_3, options, date0_3, A_3, m_3, coe0_3, 'Short Term RADARSAT 2');
% fprintf('Time it took to propagate Variation of Parameter: %f [s]\n\n',...
% timer_s_3)
% profile viewer
% [v1_s_3, v2_s_3] = lambert(rnew_s_3(end,:)', rnew_np_s_3, dt_3);
% delta_v_s_1_3 = vnew_s_3(end,:)' - v1_s_3;
% delta_v_s_2_3 = v2_s_3 - vnew_np_s_3;
% delta_v_s_3 = norm(delta_v_s_1_3) +  norm(delta_v_s_2_3);
% fprintf('Delta-V to Correct %f [km/s]\n\n',delta_v_s_3)
% 
% % Long term
% disp('Long Term')
% profile on
% [tnew_L_3, rnew_L_3, vnew_L_3, statenew_L_3, diff_COE_L_3, timer_L_3] = proj2_solver_VoP(@VoP, tspan_long,...
%     state_3, options, date0_3, A_3, m_3, coe0_3, 'Long Term RADARSAT 2');
% fprintf('Time it took to propagate Variation of Parameter: %f [s]\n\n',...
% timer_L_3)
% profile viewer
% [v1_L_3, v2_L_3] = lambert(rnew_L_3(end,:)', rnew_np_L_3, dt_3);
% delta_v_L_1_3 = vnew_L_3(end,:)' - v1_L_3;
% delta_v_L_2_3 = v2_L_3 - vnew_np_s_3;
% delta_v_L_3 = norm(delta_v_L_1_3) +  norm(delta_v_L_2_3);
% fprintf('Delta-V to Correct %f [km/s]\n\n',delta_v_L_3)
% 
