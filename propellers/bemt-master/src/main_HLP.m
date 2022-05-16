% Customized version of main.m for the high lift propellers
% 
%

%%
% Initialization
clear variables; close all; clc;
addpath(genpath('.'));  % Add all functions in subfolders to matlab path

metric = 0; % 1 = metric units on the graphs / 0 = imperial units
fig_size = 0.75;

%
if metric == 1
    unit.m2in = 1;
    unit.N2lbf = 1;
    unit.mps2ftps = 1;
    unit.mps2kt = 1;
    
    unit.m2in_txt = '[m]';
    unit.N2lbf_txt = '[N]';
    unit.mps2ftps_txt = '[m/s]';
    unit.mps2kt_txt = '[m/s]';
else
    unit.m2in = 1/0.0254;
    unit.N2lbf = 1/4.448222;
    unit.mps2ftps = 3.28084;
    unit.mps2kt = 0.514444;
    
    unit.m2in_txt = '[in]';
    unit.N2lbf_txt = '[lbf]';
    unit.mps2ftps_txt = '[ft/s]';
    unit.mps2kt_txt = '[kt]';
    
end


%%
disp('===================================================================')
disp('HLP')
disp('===================================================================')
% Additional parameters
HLP.target_V_i = 10; % [m/s] target induced velocity
%HLP.target_V_i = 10.7; % [m/s] target induced velocity
HLP.to_vel = 18.5; % [m/s] takeoff velocity with blown wing effect
%HLP.to_vel = 20; % [m/s] takeoff velocity with blown wing effect

% Configuration
CONFIG_FILE = 'configurations/FRY_HLP.m';
checkConfig(CONFIG_FILE);
run(CONFIG_FILE);

Param.Air.AXIAL_VELOC = HLP.to_vel;
Param.Air.ALTITUDE = 0;
%Param.Air.ALTITUDE = 5000/3.28084;
Blades.OMEGArpm = 4000;
%Blades.OMEGArpm = 4325;

radiusChecker(Blades, 0.7)

% Axial induction factor
HLP.a = HLP.target_V_i / HLP.to_vel;
disp(['a = ', num2str(HLP.a, '%.2f')]);

% Checking the tangential induction factor (should be maximum 0.5)
HLP.a_p = (1 - sqrt(1-(4*Param.Air.AXIAL_VELOC^2*(1+HLP.a)*HLP.a)/((Blades.OMEGArpm /60*2*pi)^2*Blades.RADIUS^2)))/2;
disp(['a_p = ', num2str(HLP.a_p, '%.2f')]);

% Run solver based on CONFIG_FILE data
warning on;
[Results, Elem, Param, Blades] = bemt(Param, Blades);
warning off;

% Resulting 
HLP.res.a = trapz(Elem.Ypos, Elem.InducedVeloc / HLP.to_vel) / (Elem.Ypos(end) - Elem.Ypos(1));
disp(['resulting a = ', num2str(HLP.res.a, '%.2f')]);

HLP.res.ind_vel = trapz(Elem.Ypos, Elem.InducedVeloc) / (Elem.Ypos(end) - Elem.Ypos(1)); % average
disp(['resulting average induced velocity = ', num2str(HLP.res.ind_vel, '%.2f')]);

% Plots
figure('Name', 'HLP design');
set(gcf,'units','points','position',[1500,0,1280,1000])

subplot(2,2,1)
plot(Elem.Ypos * unit.m2in, Elem.InducedVeloc * unit.mps2ftps); hold on;
title('Name', 'Induced velocity as a fct of radius')
xlabel(['Radius ', unit.m2in_txt])
ylabel(['Induced velocity ', unit.mps2ftps_txt])
grid on;

subplot(2,2,2)
plot(Elem.Ypos * unit.m2in, Elem.InducedVeloc / HLP.to_vel, 'DisplayName', 'Axial'); hold on;
plot(Elem.Ypos * unit.m2in, (1 - (1-(4.*Param.Air.AXIAL_VELOC.^2.*(1+Elem.InducedVeloc).*Elem.InducedVeloc)./((Blades.OMEGArpm /60*2*pi).^2.*Elem.Ypos.^2)).^(1/2))/2, 'DisplayName', 'Tangential'); hold on;
title('Induction factors as a fct of radius')
xlabel(['Radius ', unit.m2in_txt])
ylabel('Induction factor [-]')
grid on;
legend();

subplot(2,2,3)
plot(Elem.Ypos * unit.m2in, Elem.Chord/2 * unit.m2in, 'Color', 'red'); hold on;
plot(Elem.Ypos * unit.m2in, -Elem.Chord/2 * unit.m2in, 'Color', 'red'); hold on;
ylim([-max(Elem.Chord/2 * unit.m2in)*1.5 max(Elem.Chord/2 * unit.m2in)*1.5])
xlim([0, max(Elem.Ypos * unit.m2in)])
axis equal
title('Chord as a fct of radius')
xlabel(['Radius ', unit.m2in_txt])
ylabel(['Chord ', unit.m2in_txt])
grid on

subplot(2,2,4)
plot(Elem.Ypos * unit.m2in, rad2deg(Elem.Pitch)); hold on;
title('Pitch as a fct of radius')
xlabel(['Radius ', unit.m2in_txt])
ylabel('Pitch [-]')
grid on


% Results
printResults


%% printing the other figures
%
figure('Name', 'HLP design');
plot(Elem.Ypos * unit.m2in, Elem.InducedVeloc * unit.mps2ftps, 'Color', 'b', 'DisplayName', 'induced velocity'); hold on;
yline(HLP.target_V_i * unit.mps2ftps, 'Color', 'b', 'LineStyle', '--', 'DisplayName', 'target induced velocity')
xlabel(['Radius ', unit.m2in_txt])
ylabel(['Induced velocity ', unit.mps2ftps_txt])
grid on;
legend({'induced velocity', 'target induced velocity'}, 'Location', 'south');

fig2pdf(gcf, ['HLP_span_velocities', '_at_', num2str(round(Param.Air.ALTITUDE))], fig_size, 1.5, 'Figures_final/')

%
figure('Name', 'HLP Chord');
plot(Elem.Ypos * unit.m2in, Elem.Chord/2 * unit.m2in, 'Color', 'red'); hold on;
plot(Elem.Ypos * unit.m2in, -Elem.Chord/2 * unit.m2in, 'Color', 'red'); hold on;
ylim([-max(Elem.Chord/2 * unit.m2in)*1.5 max(Elem.Chord/2 * unit.m2in)*1.5])
xlim([0, max(Elem.Ypos * unit.m2in)])
axis equal
xlabel(['Radius ', unit.m2in_txt])
ylabel(['Chord ', unit.m2in_txt])
grid on
yticks(round(sort([0, max(Elem.Chord/2 * unit.m2in), -max(Elem.Chord/2 * unit.m2in), 2*max(Elem.Chord/2 * unit.m2in), -2*max(Elem.Chord/2 * unit.m2in)]), 2))

fig2pdf(gcf, 'HLP_Chord', fig_size, 1.5, 'Figures_final/')

%
figure('Name', 'HLP pitch');
plot(Elem.Ypos * unit.m2in, rad2deg(Elem.Pitch)); hold on;
xlabel(['Radius ', unit.m2in_txt])
ylabel('Pitch [-]')
grid on

fig2pdf(gcf, 'HLP_pitch', fig_size, 1.5, 'Figures_final/')


%% thrust, power and induced velocity as a function of the velocity
HLP.to_velocities = linspace(0,88/2,200); % range of velocities to try

HLP.res_vel = [];
HLP.res_T = [];
HLP.res_P = [];
HLP.res_ind_vel = [];

for vel = HLP.to_velocities
    Param.Air.AXIAL_VELOC = vel;
    
    try
        [Results, Elem, Param, Blades] = bemt(Param, Blades);
        
        HLP.res_vel = [HLP.res_vel, vel];
        
        HLP.res_T = [HLP.res_T, Results.T];
        HLP.res_P = [HLP.res_P, Results.P];
        HLP.res_ind_vel = [HLP.res_ind_vel, trapz(Elem.Ypos, Elem.InducedVeloc) / (Elem.Ypos(end) - Elem.Ypos(1))]; % average
        
    catch
        
    end
    
end

%
figure('Name',  'Velocity')

colororder({'b','r'})

yyaxis left
plot(HLP.res_vel * unit.mps2ftps, HLP.res_T * unit.N2lbf, 'Color', 'b'); hold on;
ylabel(['Thrust ', unit.N2lbf_txt])

xlabel(['Aircraft TAS ', unit.mps2ftps_txt])

yyaxis right
plot(HLP.res_vel * unit.mps2ftps, HLP.res_ind_vel * unit.mps2ftps, 'Color', 'r'); hold on;
ylabel(['Induced velocity ', unit.mps2ftps_txt])
yline(HLP.target_V_i * unit.mps2ftps, 'Color', 'r', 'LineWidth', 1, 'LineStyle', '--');

grid on;
xlim([0, max(HLP.res_vel * unit.mps2ftps)])
xline(HLP.to_vel * unit.mps2ftps, 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--');

legend({'thrust', 'induced velocity', 'target induced velocity', 'take-off speed'}, 'Location', 'south');

fig2pdf(gcf, ['HLP_velocities', '_at_', num2str(round(Param.Air.ALTITUDE))], fig_size, 1.5, 'Figures_final/')






