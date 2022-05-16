% Customized version of main.m for the wingtip propellers
% 
%

% Initialization
clear variables; close all; clc; warning off;
addpath(genpath('.'));  % Add all functions in subfolders to matlab path

metric = 0; % 1 = metric units on the graphs / 0 = imperial units
fig_size = 0.75;

%
if metric == 1
    unit.m2in = 1;
    unit.N2lbf = 1;
    unit.mps2ftps = 1;
    unit.mps2kt = 1;
    unit.kW2hp = 1;
    
    unit.m2in_txt = '[m]';
    unit.N2lbf_txt = '[N]';
    unit.mps2ftps_txt = '[m/s]';
    unit.mps2kt_txt = '[m/s]';
    unit.kW2hp_txt = '[kW]';
else
    unit.m2in = 1/0.0254;
    unit.N2lbf = 1/4.448222;
    unit.mps2ftps = 3.28084;
    unit.mps2kt = 0.514444;
    unit.kW2hp = 1e3 / 745.7;
    
    unit.m2in_txt = '[in]';
    unit.N2lbf_txt = '[lbf]';
    unit.mps2ftps_txt = '[ft/s]';
    unit.mps2kt_txt = '[kt]';
    unit.kW2hp_txt = '[kW]';
end

%% Cruise
disp('===================================================================')
disp('WTP Cruise')
disp('===================================================================')
% Additional parameters
WTP_cr.target_thrust = 1707/2; % [N] cruise thrust target
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
plot(Elem.Ypos * unit.m2in, rad2deg(Elem.Pitch), 'DisplayName', 'cruise'); hold on;
xlabel(['Radius ', unit.m2in_txt])
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

plot(fig.WTP_pitch_ax, Elem.Ypos * unit.m2in, rad2deg(Elem.Pitch), 'DisplayName', 'take-off', 'Color', 'red'); hold on;

fig2pdf(fig.WTP_pitch, 'WTP_pitch', fig_size, 1.5, 'Figures_final/')


figure('Name', 'WTP Chord');
plot(Elem.Ypos * unit.m2in, Elem.Chord/2 * unit.m2in, 'Color', 'red'); hold on;
plot(Elem.Ypos * unit.m2in, -Elem.Chord/2 * unit.m2in, 'Color', 'red'); hold on;
ylim([-max(Elem.Chord/2 * unit.m2in)*1.5 max(Elem.Chord/2 * unit.m2in)*1.5])
xlim([0, max(Elem.Ypos * unit.m2in)])
axis equal
xlabel(['Radius ', unit.m2in_txt])
ylabel(['Chord ', unit.m2in_txt])
grid on
yticks(round(sort([0,max(Elem.Chord/2 * unit.m2in), -min(Elem.Chord/2 * unit.m2in), 2*max(Elem.Chord/2 * unit.m2in), -2*max(Elem.Chord/2 * unit.m2in)]), 2))

fig2pdf(gcf, 'WTP_Chord', fig_size, 1.5, 'Figures_final/')

%% Takeoff
disp('===================================================================')
disp('WTP Takeoff 5000 ft')
disp('===================================================================')
% Additional parameters
WTP_to.HLP_trhust = 250*12; % [N] estimated thrust 
WTP_to.target_thrust = (11200 - WTP_to.HLP_trhust)/2; % [N] takeoff thrust target
disp(['Target thrust = ', num2str(WTP_to.target_thrust,'%.2f'), ' N'])

% Configuration
WTP_to.CONFIG_FILE = 'configurations/FRY_WTP.m';
checkConfig(WTP_to.CONFIG_FILE);
run(WTP_to.CONFIG_FILE);

% Parameters specific to flight conditions
Param.Air.AXIAL_VELOC = 20;    % Axial velocity of freestream (i.e. advance speed for propeller, climb velocity for helicopter), [m/s]
Blades.OMEGArpm = 2800;     % Rotational speed, [RPM]
Param.Air.ALTITUDE = 5000/3.2808;             % Altitude, [m]

% Collective pitch
Blades.COLL_PITCHdeg = 7.5;

% % Finding the optimal collective pitch
% warning on;
% [Results, Elem, Param, Blades] = bemt(Param, Blades);
% warning off;

disp('--------------------------------------')
% Finding the optimal collective pitch
[Results, Elem, Param, Blades] = optiRotationalVel(Param, Blades, WTP_to.target_thrust, 2500, 2950);

disp(['Optimal OMEGArpm = ', num2str(Blades.OMEGArpm,'%.0f'), ' RPM']);
radiusChecker(Blades, 0.8); % checking the Mach number

% Results
printResults

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
        [Results, Elem, Param, Blades] = optiRotationalVel(Param, Blades, WTP_to.target_thrust, 2500, 2950);

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
plot(mean_to.vel * unit.mps2ftps, mean_to.rot_vel);
xlabel(['Aircraft velocity ', unit.mps2ftps_txt]);
ylabel('Rotational speed [RPM]');
grid on;
xlim([min(mean_to.vel * unit.mps2ftps), max(mean_to.vel * unit.mps2ftps)]);


%
figure('Name', 'take-off velocity vs eff and power')
yyaxis left
plot(mean_to.vel * unit.mps2ftps, mean_to.res_P/1000 * unit.kW2hp);
xlabel(['Aircraft velocity ', unit.mps2ftps_txt]);
ylabel(['Power ', unit.kW2hp_txt]);

yyaxis right
plot(mean_to.vel * unit.mps2ftps, mean_to.eff);
ylabel('Efficiency [%]');

grid on;
xlim([min(mean_to.vel * unit.mps2ftps), max(mean_to.vel * unit.mps2ftps)]);



