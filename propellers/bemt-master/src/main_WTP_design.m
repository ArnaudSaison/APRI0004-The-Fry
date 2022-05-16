% FRY_WTP (mh114) Configuration file for BEMT solver.
% 
%

% Initialization
clear variables; close all; clc;
addpath(genpath('.'));  % Add all functions in subfolders to matlab path


%% Configuration
CONFIG_FILE = 'configurations/FRY_WTP.m';
checkConfig(CONFIG_FILE);
run(CONFIG_FILE);

config = 'c'; % t = takeoff, c = cruise
metric = 0; % 1 = metric units on the graphs / 0 = imperial units
fig_size = 0.75;

%
if config == 't'
    axial_vel = 18.5;
    rpm = 2802;
    thrust = (9900 - 12 * 240)/2;
    altitude = 0;
    
    legend_pos_Collective_pitch = 'northeast';
    
elseif config == '5'
    axial_vel = 20;
    rpm = 2914;
    thrust = (11200 - 12 * 250)/2;
    altitude = 5000/3.28084;
    
    legend_pos_Collective_pitch = 'northeast';
else
    axial_vel = kts2mps(170);
    rpm = 2000;
    thrust = 1707/2;
    altitude = 3962.4;
    
    legend_pos_Collective_pitch = 'southeast';
end

if metric == 1
    unit.m2in = 1;
    unit.N2lbf = 1;
    unit.mps2ftps = 1;
    
    unit.m2in_txt = '[m]';
    unit.N2lbf_txt = '[N]';
    unit.mps2ftps_txt = '[m/s]';
else
    unit.m2in = 0.0254;
    unit.N2lbf = 1/4.448222;
    unit.mps2ftps = 3.28084;
    
    unit.m2in_txt = '[in]';
    unit.N2lbf_txt = '[lbf]';
    unit.mps2ftps_txt = '[ft/s]';
end

Param.Air.AXIAL_VELOC = axial_vel;    % Axial velocity of freestream (i.e. advance speed for propeller, climb velocity for helicopter), [m/s]
Blades.OMEGArpm = rpm;     % Rotational speed, [RPM]
Param.Air.ALTITUDE = altitude;

%% Collective pitch Cruise
% Finds the optimal collective pitch from efficiency standpoint

res_eff = []; % resulting efficiency values
res_thrust = []; % resulting thrussts
col_vals = []; % corresponding values of collective pitch
for col_val = linspace(0, 90, 200)
    try
        Blades.COLL_PITCHdeg = col_val;
        
        % Run solver based on CONFIG_FILE data
        [Results, Elem, Param, Blades] = bemt(Param, Blades);
        
        res_eff = [res_eff Results.eff_prop];
        res_thrust = [res_thrust Results.T];
        col_vals = [col_vals col_val];
    catch
        
    end
end

res_eff = coeffCleanup(res_eff);
[~, opt_eff_i] = max(res_eff);

[~, opt_thrust_i] = max(res_thrust);

% %
% figure('Name', 'Efficiency as a fct of collective pitch')
% plot(col_vals, res_eff); hold on;
% plot(col_vals(opt_eff_i), res_eff(opt_eff_i), 'Marker', 'x', 'MarkerSize', 12); hold on;
% xlabel('Collective pitch [°]')
% ylabel('Efficiency \eta [-]')
% grid on;
% 
% %
% figure('Name',  'Thrust as a fct of collective pitch')
% plot(col_vals, res_thrust); hold on;
% plot(col_vals(opt_thrust_i), res_thrust(opt_thrust_i), 'Marker', 'x', 'MarkerSize', 12); hold on;
% yline(thrust, 'Color', 'red', 'LineWidth', 2);
% xlabel('Collective pitch [°]')
% ylabel('Thrust [N]')
% grid on;

%
figure('Name',  'Collective pitch')

colororder({'b','r'})

yyaxis left
plot(col_vals, res_eff, 'Color', 'b'); hold on;
ylabel('Efficiency \eta [-]')
ylim([0,1]);

xlabel('Collective pitch [°]')

yyaxis right
plot(col_vals, res_thrust * unit.N2lbf, 'Color', 'r'); hold on;
yline(thrust * unit.N2lbf, 'Color', 'r', 'LineWidth', 1, 'LineStyle', '--');
ylabel(['Thrust ', unit.N2lbf_txt])
ylim([0, max(res_thrust * unit.N2lbf)*1.1])

grid on;
xlim([0, max(col_vals)])
legend({'effeciency', 'thrust', 'target thrust'}, 'Location', legend_pos_Collective_pitch);

fig2pdf(gcf, ['WTP_Collective_pitch_', config, '_at_', num2str(round(Param.Air.ALTITUDE))], fig_size, 1.5, 'Figures_final/')


%% Advance Ratio
% Compares different collective pitches

figure('Name', 'advance ratio')

val1_range = linspace(0, 25, 5);

for val1 = val1_range
    res = [];
    res2 = [];
    vals2 = [];
    for val2 = kts2mps(linspace(30, 250, 200))
        try
            Blades.COLL_PITCHdeg = val1;
            Param.Air.AXIAL_VELOC = val2;

            % Run solver based on CONFIG_FILE data
            [Results, Elem, Param, Blades] = bemt(Param, Blades);
            
            res = [res Results.eff_prop];
            res2 = [res2 Results.adv_ratio];
            vals2 = [vals2 val2];
            
            if Results.eff_prop < 0
                break
            end
            
        catch
            
        end
    end
    
    res = coeffCleanup(res);
    plot(res2, res); hold on;
end

xline(axial_vel/(Blades.RADIUS * Blades.OMEGA), 'LineWidth', 2, 'Color', 'red');
grid on;
legend_text = append('\theta = ', strjust(append(string(num2str((val1_range)')), '°'), 'left'));
legend_text(end+1) = [num2str(mps2kts(axial_vel), '%.0f'), ' kts'];
legend(legend_text, 'Location', 'Southeast');
xlabel('Advance ratio J [-]')
xlim([0, max(res2)])
ylabel('Efficiency \eta [-]')

fig2pdf(gcf, ['WTP_Advance_ratio_', config, '_at_', num2str(round(Param.Air.ALTITUDE))], fig_size, 1.5, 'Figures_final/')


