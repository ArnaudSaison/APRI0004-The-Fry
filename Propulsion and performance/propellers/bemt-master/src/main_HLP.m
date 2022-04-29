% Customized version of main.m for the high lift propellers
% 
%

%%
% Initialization
clear variables; close all; clc;
addpath(genpath('.'));  % Add all functions in subfolders to matlab path


%%
disp('===================================================================')
disp('HLP')
disp('===================================================================')
% Additional parameters
HLP.target_V_i = 10; % [m/s] target induced velocity
HLP.to_vel = 18.5; % [m/s] takeoff velocity with blown wing effect

% Configuration
CONFIG_FILE = 'configurations/FRY_HLP.m';
checkConfig(CONFIG_FILE);
run(CONFIG_FILE);

Param.Air.AXIAL_VELOC = HLP.to_vel;
Param.Air.ALTITUDE = 0;
Blades.OMEGArpm = 4000;

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
plot(Elem.Ypos, Elem.InducedVeloc); hold on;
title('Name', 'Induced velocity as a fct of radius')
xlabel('Radius [m]')
ylabel('Induced velocity [m/s]')
grid on;

subplot(2,2,2)
plot(Elem.Ypos, Elem.InducedVeloc / HLP.to_vel, 'DisplayName', 'Axial'); hold on;
plot(Elem.Ypos, (1 - (1-(4.*Param.Air.AXIAL_VELOC.^2.*(1+Elem.InducedVeloc).*Elem.InducedVeloc)./((Blades.OMEGArpm /60*2*pi).^2.*Elem.Ypos.^2)).^(1/2))/2, 'DisplayName', 'Tangential'); hold on;
title('Induction factors as a fct of radius')
xlabel('Radius [m]')
ylabel('Induction factor [-]')
grid on;
legend();

subplot(2,2,3)
plot(Elem.Ypos, Elem.Chord, 'Color', 'red'); hold on;
plot(Elem.Ypos, -Elem.Chord, 'Color', 'red'); hold on;
ylim([-max(Elem.Chord)*1.5 max(Elem.Chord)*1.5])
xlim([0, max(Elem.Ypos)])
axis equal
title('Chord as a fct of radius')
xlabel('Radius [m]')
ylabel('Chord [m]')
grid on

subplot(2,2,4)
plot(Elem.Ypos, rad2deg(Elem.Pitch)); hold on;
title('Pitch as a fct of radius')
xlabel('Radius [m]')
ylabel('Pitch [-]')
grid on


% Results
printResults


%% printing the other figures

figure('Name', 'HLP design');
plot(Elem.Ypos, Elem.InducedVeloc, 'Color', 'b', 'DisplayName', 'induced velocity'); hold on;
yline(HLP.target_V_i, 'Color', 'b', 'LineStyle', '--', 'DisplayName', 'target induced velocity')
xlabel('Radius [m]')
ylabel('Induced velocity [m/s]')
grid on;
legend({'induced velocity', 'target induced velocity'}, 'Location', 'south');

fig2pdf(gcf, 'HLP_span_velocities', 1, 1.5, 'Figures/')

figure('Name', 'HLP Chord');
plot(Elem.Ypos, Elem.Chord, 'Color', 'red'); hold on;
plot(Elem.Ypos, -Elem.Chord, 'Color', 'red'); hold on;
ylim([-max(Elem.Chord)*1.5 max(Elem.Chord)*1.5])
xlim([0, max(Elem.Ypos)])
axis equal
xlabel('Radius [m]')
ylabel('Chord [m]')
grid on

fig2pdf(gcf, 'HLP_Chord', 1, 1.5, 'Figures/')

figure('Name', 'HLP pitch');
plot(Elem.Ypos, rad2deg(Elem.Pitch)); hold on;
xlabel('Radius [m]')
ylabel('Pitch [-]')
grid on

fig2pdf(gcf, 'HLP_pitch', 1, 1.5, 'Figures/')


%% thrust, power and induced velocity as a function of the velocity
HLP.to_velocities = linspace(0,88,100); % range of velocities to try

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
plot(HLP.res_vel, HLP.res_T, 'Color', 'b'); hold on;
ylabel('Thrust [N]')

xlabel('Aircraft TAS [m/s]')

yyaxis right
plot(HLP.res_vel, HLP.res_ind_vel, 'Color', 'r'); hold on;
ylabel('Induced velocity [m/s]')
yline(HLP.target_V_i, 'Color', 'r', 'LineWidth', 1, 'LineStyle', '--');

grid on;
xlim([0, max(HLP.res_vel)])
xline(HLP.to_vel, 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--');

legend({'thrust', 'induced velocity', 'target induced velocity', 'take-off speed'}, 'Location', 'northeast');

fig2pdf(gcf, 'HLP_velocities', 1, 1.5, 'Figures/')






