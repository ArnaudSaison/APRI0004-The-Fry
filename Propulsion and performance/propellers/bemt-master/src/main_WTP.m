% Customized version of main.m for the wingtip propellers
% 
%

% Initialization
clear variables; close all; clc; warning off;
addpath(genpath('.'));  % Add all functions in subfolders to matlab path

%% Cruise
disp('===================================================================')
disp('WTP Cruise')
disp('===================================================================')
% Additional parameters
WTP_cr.target_thrust = 1650/2; % [N] cruise thrust target
disp(['Target thrust = ', num2str(WTP_cr.target_thrust,'%.2f'), ' N'])

% Configuration
WTP_cr.CONFIG_FILE = 'configurations/FRY_WTP.m';
checkConfig(WTP_cr.CONFIG_FILE);
run(WTP_cr.CONFIG_FILE);

% Parameters specific to flight conditions
Param.Air.AXIAL_VELOC   = kts2mps(170);    % Axial velocity of freestream (i.e. advance speed for propeller, climb velocity for helicopter), [m/s]
Blades.OMEGArpm = 2000;     % Rotational speed, [RPM]

radiusChecker(Blades, 0.7); % checking the Mach number

disp('--------------------------------------')
% Finding the optimal collective pitch
[Results, Elem, Param, Blades] = optiCollective(Param, Blades, WTP_cr.target_thrust, 0, 40);

disp(['Optimal collective pitch = ', num2str(Blades.COLL_PITCHdeg,'%.2f'), '°']);

% Results
printResults

fig.WTP_pitch = figure('Name', 'WTP pitch');
fig.WTP_pitch_ax = axes();
plot(Elem.Ypos, rad2deg(Elem.Pitch), 'DisplayName', 'cruise'); hold on;
xlabel('Radius [m]')
ylabel('Pitch [-]')
grid on
legend();


%% Takeoff
disp('===================================================================')
disp('WTP Takeoff')
disp('===================================================================')
% Additional parameters
WTP_to.HLP_trhust = 240*12; % [N] estimated thrust 
WTP_to.target_thrust = (9900 - WTP_to.HLP_trhust)/2; % [N] takeoff thrust target
disp(['Target thrust = ', num2str(WTP_to.target_thrust,'%.2f'), ' N'])

% Configuration
WTP_to.CONFIG_FILE = 'configurations/FRY_WTP.m';
checkConfig(WTP_to.CONFIG_FILE);
run(WTP_to.CONFIG_FILE);

% Parameters specific to flight conditions
Param.Air.AXIAL_VELOC = 18.5;    % Axial velocity of freestream (i.e. advance speed for propeller, climb velocity for helicopter), [m/s]
Blades.OMEGArpm = 2800;     % Rotational speed, [RPM]
Param.Air.ALTITUDE = 0;             % Altitude, [m]

% Collective pitch
Blades.COLL_PITCHdeg = 0;

% % Finding the optimal collective pitch
% warning on;
% [Results, Elem, Param, Blades] = bemt(Param, Blades);
% warning off;

disp('--------------------------------------')
% Finding the optimal collective pitch
[Results, Elem, Param, Blades] = optiRotationalVel(Param, Blades, WTP_to.target_thrust, 2500, 2900);

disp(['Optimal OMEGArpm = ', num2str(Blades.OMEGArpm,'%.0f'), ' RPM']);
radiusChecker(Blades, 0.8); % checking the Mach number

% Results
printResults

%% Geometry

plot(fig.WTP_pitch_ax, Elem.Ypos, rad2deg(Elem.Pitch), 'DisplayName', 'take-off'); hold on;

fig2pdf(fig.WTP_pitch, 'WTP_pitch', 1, 1.5, 'Figures/')


figure('Name', 'WTP Chord');
plot(Elem.Ypos, Elem.Chord, 'Color', 'red'); hold on;
plot(Elem.Ypos, -Elem.Chord, 'Color', 'red'); hold on;
ylim([-max(Elem.Chord)*1.5 max(Elem.Chord)*1.5])
xlim([0, max(Elem.Ypos)])
axis equal
xlabel('Radius [m]')
ylabel('Chord [m]')
grid on

fig2pdf(gcf, 'WTP_Chord', 1, 1.5, 'Figures/')

%% Takeoff
disp('===================================================================')
disp('WTP throughout take-off')
disp('===================================================================')
mean_to.vel = [];
mean_to.res_T = [];
mean_to.rot_vel = [];
mean_to.eff = [];
mean_to.res_P = [];

for vel = linspace(0,18.5,40)
    try
        % Additional parameters
        WTP_to.HLP_trhust = 240*12; % [N] estimated thrust 
        WTP_to.target_thrust = (9900 - WTP_to.HLP_trhust)/2; % [N] takeoff thrust target
        disp(['Target thrust = ', num2str(WTP_to.target_thrust,'%.2f'), ' N'])

        % Configuration
        WTP_to.CONFIG_FILE = 'configurations/FRY_WTP.m';
        checkConfig(WTP_to.CONFIG_FILE);
        run(WTP_to.CONFIG_FILE);

        % Parameters specific to flight conditions
        Param.Air.AXIAL_VELOC = vel;    % Axial velocity of freestream (i.e. advance speed for propeller, climb velocity for helicopter), [m/s]
        Blades.OMEGArpm = 2800;     % Rotational speed, [RPM]
        Param.Air.ALTITUDE = 0;             % Altitude, [m]

        % Collective pitch
        Blades.COLL_PITCHdeg = 0;

        % % Finding the optimal collective pitch
        % warning on;
        % [Results, Elem, Param, Blades] = bemt(Param, Blades);
        % warning off;

        % Finding the optimal collective pitch
        [Results, Elem, Param, Blades] = optiRotationalVel(Param, Blades, WTP_to.target_thrust, 2500, 2900);

        % results
        mean_to.vel = [mean_to.vel, vel];
        mean_to.res_T = [mean_to.res_T, Results.T];
        mean_to.res_P = [mean_to.res_P, Results.P];
        mean_to.rot_vel = [mean_to.rot_vel, Blades.OMEGArpm];
        mean_to.eff = [mean_to.eff, Results.eff_prop*100];
    catch
        
    end

end

mean_to.mean_T = mean(mean_to.res_T);

disp('--------------------------------------')
disp(['Mean thrust = ', num2str(mean_to.mean_T)])

%
figure('Name', 'take-off velocity vs rotational speed')
plot(mean_to.vel, mean_to.rot_vel);
xlabel('Aircraft velocity [m/s]');
ylabel('Rotational speed [RPM]');
grid on;
xlim([min(mean_to.vel), max(mean_to.vel)]);


%
figure('Name', 'take-off velocity vs eff and power')
yyaxis left
plot(mean_to.vel, mean_to.res_P/1000);
xlabel('Aircraft velocity [m/s]');
ylabel('Power [kW]');

yyaxis right
plot(mean_to.vel, mean_to.eff);
ylabel('Efficiency [%]');

grid on;
xlim([min(mean_to.vel), max(mean_to.vel)]);



