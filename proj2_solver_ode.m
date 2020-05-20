%% Project 2 Solver for ode45
% Solves ode45 and plots
function [t_new,state_new,diff_COE,timer] = proj2_solver_ode(fun,tspan,...
    state,options,date,A,m,initial_COE,title_name)
global mue re
date_1 = datetime(date);
jd = juliandate(date);
tic % Timer start
[t_new,state_new] = ode45(fun,tspan,state,options,date,date_1,jd,A,m); % ode45
timer = toc; % Timer end

r_new = [state_new(:,1),state_new(:,2),state_new(:,3)]; % Position vectors [km]
v_new = [state_new(:,4),state_new(:,5),state_new(:,6)]; % Velocity vectors [km]

% Pre-allocate
data_COE = zeros(length(r_new),length(initial_COE));
diff_COE = zeros(length(r_new),length(initial_COE));

% COE difference calculation
for i = 1:length(r_new)
[a0,E0,H0,inc0,RAAN0,omega0,theta0] = coe(r_new(i,:),v_new(i,:));
ra0 = (H0^2)/(mue*(1-E0)); % Radius of apoapse
rp0 = (2*a0)-ra0; % Radius of periapse
current_COE = [rp0-re,ra0-re,a0,E0,H0,inc0,RAAN0,omega0,theta0];
data_COE(i,:) = current_COE; % Current COE stored
diff_COE(i,:) = current_COE-initial_COE; % COE difference from inital
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
title(strcat({title_name,'Orbit'}));
xlabel('x [km]')
ylabel('y [km]')
zlabel('z [km]')
grid on
axis equal

figure
plot(t_new/(24*60*60),data_COE(:,1))
hold on
plot(t_new/(24*60*60),data_COE(:,2))
title(strcat({title_name,'Altitude'}));
xlabel('Time [day]')
ylabel('Altitude [km]')
legend('Periapse','Apoapse')
grid on

figure
plot(t_new/(24*60*60),diff_COE(:,4));
title(strcat({title_name,'Eccentricity'}));
xlabel('Time [day]')
ylabel('Eccentricity')
grid on

figure
plot(t_new/(24*60*60),diff_COE(:,6));
title(strcat({title_name,'Inclination'}));
xlabel('Time [day]')
ylabel('Inclination [^o]')
grid on

figure
plot(t_new/(24*60*60),diff_COE(:,7));
title(strcat({title_name,'Right Ascension of Ascending Node'}));
xlabel('Time [day]')
ylabel('\Omega [^o]')
grid on

figure
plot(t_new/(24*60*60),diff_COE(:,8));
title(strcat({title_name,'Arguement of Periapse'}));
xlabel('Time [day]')
ylabel('\omega [^o]')
grid on

end

