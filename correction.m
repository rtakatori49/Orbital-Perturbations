function [tnew_c_s, rnew_c_s, vnew_c_s] = correction(state, state_coe, options, date0, A, m, coe0, dt_c_s)
global duration re mue
t_c_s = 7*24*60*60;
tspan_c_s = [0 t_c_s];
if duration == 0
for i = 1:4
    tspan_np = [0 t_c_s + dt_c_s];
    [tnew_np, statenew_np] = ode45(@two_body, tspan_np, state, options); % ode45
        [tnew_c_s{i}, rnew_c_s{i}, vnew_c_s{i}, statenew_c_s{i}, diff_COE_c_s{i}, timer_c_s] = proj2_solver_VoP(@VoP, tspan_c_s,...
        state_coe, options, date0, A, m, coe0, 'BSAT-3C Corrections');
%     [tnew_c_s{i}, statenew_c_s{i}, diff_COE_c_s{i}, timer_c_s] = proj2_solver_ode(@cowell, tspan_c_s,...
%         state, options, date0, A, m, coe0, 'BSAT-3C Corrections');
%     rnew_c_s{i} = statenew_c_s{i}(:,1:3);
%     vnew_c_s{i} = statenew_c_s{i}(:,4:6);
    fprintf('Time it took to propagate VoP: %f [s]\n\n',...
        timer_c_s)
    [v1_c_s, v2_c_s] = lambert(rnew_c_s{i}(end,:)', statenew_np(end,1:3)', dt_c_s);
    delta_v_c_s_1 = vnew_c_s{i}(end,:)' - v1_c_s;
    delta_v_c_s_2 = v2_c_s - statenew_np(end,4:6)';
    delta_v_c_s(i) = norm(delta_v_c_s_1) +  norm(delta_v_c_s_1);
    fprintf('Delta-V to Correct %f [km/s]\n\n',delta_v_c_s(i))
    state_transfer = [rnew_c_s{i}(end,:)'; v1_c_s];
    tspan_transfer = [tnew_c_s{i}(end) tnew_c_s{i}(end) + dt_c_s];
    [tnew_transfer{i}, statenew_transfer{i}] = ode45(@two_body, tspan_transfer, state_transfer, options); % ode45
    rnew_transfer{i} = statenew_transfer{i}(:,1:3);
    vnew_transfer{i} = statenew_transfer{i}(:,4:6);

    [a, ecc, h, inc, RAAN, omega, theta] = coe(rnew_transfer{i}(end,:)', statenew_np(end,4:6)');
    rp = a * (1 - ecc); % Perigee [km]
    ra = (2 * a) - rp; % Apogee [km]
    zp = rp - re; % Perigee altitude [km]
    za = ra - re; % Apogee altitude [km]
    tspan_c_s = [tnew_transfer{i}(end) tnew_transfer{i}(end) + t_c_s];
    state = [rnew_transfer{i}(end,:) statenew_np(end,4:6)];
    coe0 = [zp, za, a, ecc, h, inc, RAAN, omega, theta];
    % In radians
    RAAN_rad = deg2rad(RAAN); % RAAN [rad]
    inc_rad = deg2rad(inc); % Inclination [rad]
    omega_rad = deg2rad(omega); % Arguement of perigee [rad]
    theta_rad = deg2rad(theta); % True Anomaly [rad]
    state_coe = [h,ecc,theta_rad,inc_rad,RAAN_rad,omega_rad]'; % State vector
end

figure
plot3(rnew_c_s{1}(:,1), rnew_c_s{1}(:,2), rnew_c_s{1}(:,3))
hold on
plot3(rnew_transfer{1}(:,1), rnew_transfer{1}(:,2), rnew_transfer{1}(:,3))
plot3(rnew_c_s{2}(:,1), rnew_c_s{2}(:,2), rnew_c_s{2}(:,3))
plot3(rnew_transfer{2}(:,1), rnew_transfer{2}(:,2), rnew_transfer{2}(:,3))
plot3(rnew_c_s{3}(:,1), rnew_c_s{3}(:,2), rnew_c_s{3}(:,3))
plot3(rnew_transfer{3}(:,1), rnew_transfer{3}(:,2), rnew_transfer{3}(:,3))
plot3(rnew_c_s{4}(:,1), rnew_c_s{4}(:,2), rnew_c_s{4}(:,3))
plot3(rnew_transfer{4}(:,1), rnew_transfer{4}(:,2), rnew_transfer{4}(:,3))
[x,y,z] = sphere;
globe = surf(x*re,y*re,z*re);
cdata = imread('earth_base_high_1.jpg');
sph1 = findobj('Type', 'surface');
set(globe,'FaceColor','texturemap','CData',cdata,'EdgeColor','none');
title('Correction Orbit')
xlabel('x [km]')
ylabel('y [km]')
zlabel('z [km]')
legend('Week 1 Orbit', 'Transfer 1','Week 2 Orbit', 'Transfer 2', 'Week 3 Orbit',...
    'Transfer 3', 'Week 4 Orbit', 'Transfer 4')
grid on
axis equal

figure
plot(tnew_c_s{1}/(24*60*60), diff_COE_c_s{1}(:, 1))
hold on
plot(tnew_c_s{1}/(24*60*60), diff_COE_c_s{1}(:, 2))
plot(tnew_c_s{2}/(24*60*60), diff_COE_c_s{2}(:, 1))
plot(tnew_c_s{2}/(24*60*60), diff_COE_c_s{2}(:, 2))
plot(tnew_c_s{3}/(24*60*60), diff_COE_c_s{3}(:, 1))
plot(tnew_c_s{3}/(24*60*60), diff_COE_c_s{3}(:, 2))
plot(tnew_c_s{4}/(24*60*60), diff_COE_c_s{4}(:, 1))
plot(tnew_c_s{4}/(24*60*60), diff_COE_c_s{4}(:, 2))
title(strcat({'BSAT-3C Corrections', 'Altitude'}));
xlabel('Time [day]')
ylabel('Altitude [km]')
legend('Periapse','Apoapse')
legend('Week 1 Periapse', 'Week 1 Apoapse', 'Week 2 Periapse',...
'Week 2 Apoapse', 'Week 3 Periapse', 'Week 3 Apoapse','Week 4 Periapse',...
'Week 4 Apoapse')

grid on

figure
plot(tnew_c_s{1}/(24*60*60), diff_COE_c_s{1}(:, 4));
hold on
plot(tnew_c_s{2}/(24*60*60), diff_COE_c_s{2}(:, 4));
plot(tnew_c_s{3}/(24*60*60), diff_COE_c_s{3}(:, 4));
plot(tnew_c_s{4}/(24*60*60), diff_COE_c_s{4}(:, 4));
title(strcat({'BSAT-3C Corrections', 'Eccentricity'}));
xlabel('Time [day]')
ylabel('Eccentricity')
legend('Week 1', 'Week 2', 'Week 3', 'Week 4')
grid on

figure
plot(tnew_c_s{1}/(24*60*60), diff_COE_c_s{1}(:, 6));
hold on
plot(tnew_c_s{2}/(24*60*60), diff_COE_c_s{2}(:, 6));
plot(tnew_c_s{3}/(24*60*60), diff_COE_c_s{3}(:, 6));
plot(tnew_c_s{4}/(24*60*60), diff_COE_c_s{4}(:, 6));
title(strcat({'BSAT-3C Corrections', 'Inclination'}));
xlabel('Time [day]')
ylabel('Inclination [^o]')
legend('Week 1', 'Week 2', 'Week 3', 'Week 4')
grid on

figure
plot(tnew_c_s{1}/(24*60*60), diff_COE_c_s{1}(:, 7));
hold on
plot(tnew_c_s{2}/(24*60*60), diff_COE_c_s{2}(:, 7));
plot(tnew_c_s{3}/(24*60*60), diff_COE_c_s{3}(:, 7));
plot(tnew_c_s{4}/(24*60*60), diff_COE_c_s{4}(:, 7));
title(strcat({'BSAT-3C Corrections', 'Right Ascension of Ascending Node'}));
xlabel('Time [day]')
ylabel('\Omega [^o]')
legend('Week 1', 'Week 2', 'Week 3', 'Week 4')
grid on

figure
plot(tnew_c_s{1}/(24*60*60), diff_COE_c_s{1}(:, 8));
hold on
plot(tnew_c_s{2}/(24*60*60), diff_COE_c_s{2}(:, 8));
plot(tnew_c_s{3}/(24*60*60), diff_COE_c_s{3}(:, 8));
plot(tnew_c_s{4}/(24*60*60), diff_COE_c_s{4}(:, 8));
title(strcat({'BSAT-3C Corrections', 'Arguement of Periapse'}));
xlabel('Time [day]')
ylabel('\omega [^o]')
legend('Week 1', 'Week 2', 'Week 3', 'Week 4')
grid on

else
end
end

