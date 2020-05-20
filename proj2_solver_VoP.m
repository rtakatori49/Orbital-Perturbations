%% Homework 4 Solver for Variation of Parameters
% Solves ode45 and plots
function [t_new, r_new, v_new, state_new, diff_COE, timer] = proj2_solver_VoP(fun,...
    tspan, state, options, date, A, m, initial_COE, title_name)
global mue re
tic % Timer start
date_1 = datetime(date);
jd = juliandate(date);
[t_new, state_new] = ode45(fun, tspan, state, options, date, date_1, jd, A, m); % ode45
timer = toc; % Timer end

% COE states
h_new = state_new(:, 1);
ecc_new = state_new(:, 2);
theta_new = state_new(:, 3);
inc_new = state_new(:, 4);
RAAN_new = state_new(:, 5);
omega_new = state_new(:, 6);

% Pre-allocate
data_COE = zeros(length(h_new), length(initial_COE));
diff_COE = zeros(length(h_new), length(initial_COE));
a_new = zeros(1, length(h_new));
inc_new_deg = zeros(1, length(h_new));
RAAN_new_deg = zeros(1, length(h_new));
omega_new_deg = zeros(1, length(h_new));
theta_new_deg = zeros(1, length(h_new));
r_new = zeros(length(h_new), 3);
v_new = zeros(length(h_new), 3);

% Calculate RV vector
for i = 1:length(h_new)
[r_peri, v_peri, r_new(i,:), v_new(i,:)] = coe2rv_rad(mue, ecc_new(i),...
    h_new(i), inc_new(i), RAAN_new(i), omega_new(i), theta_new(i));
end

% COE difference calculation
for i = 1:length(r_new)
a_new(i) = (h_new(i)^2)/(mue*(1 - ecc_new(i)^2));
inc_new_deg(i) = rad2deg(inc_new(i));
RAAN_new_deg(i) = rad2deg(RAAN_new(i));
omega_new_deg(i) = rad2deg(omega_new(i));
theta_new_deg(i) = rad2deg(theta_new(i));
ra_new = (h_new(i)^2)/(mue*(1 - ecc_new(i))); % Radius of apoapse
rp_new = (2*a_new(i)) - ra_new; % Radius of periapse
current_COE = [rp_new-re, ra_new-re, a_new(i), ecc_new(i),...
    h_new(i), inc_new_deg(i), RAAN_new_deg(i), omega_new_deg(i), theta_new_deg(i)];
data_COE(i,:) = current_COE; % Current COE stored
diff_COE(i,:) = current_COE - initial_COE; % COE difference from inital
end

% Display
figure
plot3(r_new(:, 1), r_new(:, 2), r_new(:, 3),'Color','r');
hold on
[x,y,z] = sphere;
globe = surf(x*re,y*re,z*re);
cdata = imread('earth_base_high_1.jpg');
sph1 = findobj('Type', 'surface');
set(globe,'FaceColor','texturemap','CData',cdata,'EdgeColor','none');
title(strcat({title_name, 'Orbit'}));
xlabel('x [km]')
ylabel('y [km]')
zlabel('z [km]')
grid on
axis equal

figure
plot(t_new/(24*60*60), data_COE(:, 1))
hold on
plot(t_new/(24*60*60), data_COE(:, 2))
title(strcat({title_name, 'Altitude'}));
xlabel('Time [day]')
ylabel('Altitude [km]')
legend('Periapse','Apoapse')
grid on

figure
plot(t_new/(24*60*60), diff_COE(:, 4));
title(strcat({title_name, 'Eccentricity'}));
xlabel('Time [day]')
ylabel('Eccentricity')
grid on

figure
plot(t_new/(24*60*60), diff_COE(:, 6));
title(strcat({title_name, 'Inclination'}));
xlabel('Time [day]')
ylabel('Inclination [^o]')
grid on

figure
plot(t_new/(24*60*60), diff_COE(:, 7));
title(strcat({title_name, 'Right Ascension of Ascending Node'}));
xlabel('Time [day]')
ylabel('\Omega [^o]')
grid on

figure
plot(t_new/(24*60*60), diff_COE(:, 8));
title(strcat({title_name, 'Arguement of Periapse'}));
xlabel('Time [day]')
ylabel('\omega [^o]')
grid on

end

