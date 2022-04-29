% Customized version of main.m for the high lift propellers
% 
%

%%
% Initialization
clear variables; close all; clc;
addpath(genpath('.'));  % Add all functions in subfolders to matlab path

% Additional parameters
HLP.target_V_i = 10; % [m/s] target induced velocity
HLP.to_vel = 18.5; % [m/s] takeoff velocity with blown wing effect

% Configuration
CONFIG_FILE = 'configurations/FRY_HLP.m';
checkConfig(CONFIG_FILE);
run(CONFIG_FILE);

Param.Air.AXIAL_VELOC = HLP.to_vel;
Param.Air.ALTITUDE = 0;
Blades.OMEGArpm = 4500;

radiusChecker(Blades, 0.7)

%% Theta tip
% Finds the optimal collective pitch from efficiency standpoint

res_eff = []; % resulting efficiency values
theta_vals = []; % corresponding values of collective pitch
for theta_val = linspace(0, 90, 100)
    try
        Blades.COLL_PITCHdeg = theta_val;
        
        % Run solver based on CONFIG_FILE data
        [Results, Elem, Param, Blades] = bemt(Param, Blades);
        
        res_eff = [res_eff Results.eff_prop];
        theta_vals = [theta_vals theta_val];
    catch
        disp('error')
    end
end

res_eff = coeffCleanup(res_eff);
[~, opt_eff_i] = max(res_eff);

figure('Name', 'Efficiency as a fct of collective pitch')
plot(theta_vals, res_eff); hold on;
plot(theta_vals(opt_eff_i), res_eff(opt_eff_i), 'Marker', '.', 'MarkerSize', 12); hold on;
xlabel('Collective pitch [°]')
ylabel('Efficiency \eta [-]')
grid on;





