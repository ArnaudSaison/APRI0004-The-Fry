%% Conceptual Design of a STOL Aircraft 
% (AIAA 2022 Aircraft Deisgn competition)
%
% Code by: Adrien BOURGUIGNON, Arnaud SAISON, Augustin SANDRONT, Maxence 
% CASAGRANDE, Maxime DUMONT, Veysi ASLANCI, Robin VESTRAETE, Tom DETHIER
% 
% Academic year: 2021-2022
% University: Université de Liège - Faculté des Sciences Appliquées
% Master in Aerospace Engineering
% Course: Aerospace Design Project
% 
% All paremeters set in 'parameters.m'
% All outputs defined in 'ouput.m'
% 

clear; close all; clc;

metric = 0;

if metric == 1
    unit.kW2hp = 1;
    unit.kW2hp_txt = '[kW]';
    unit.mps2ftps = 1;
    unit.mps2ftps = '[m/s]';
else
    unit.kW2hp = 1e3/745.7;
    unit.kW2hp_txt = '[hp]';
    unit.mps2ftps = 3.28084;
    unit.mps2ftps_txt = '[ft/s]';
end

%% Additional settings
C_L_blown_wing = linspace(8,5,50);
res_CLmax_to = [];
res_to_dist = [];
res_to_pow = [];
res_ind_vel = [];

%% Loop

par = parameters();


for C_L_i = C_L_blown_wing
    try
        par.takeoff_blown_wing_effect = C_L_i/par.C_L_max;
        res = performance(par);
        
        res_CLmax_to = [res_CLmax_to C_L_i];
        res_to_dist = [res_to_dist par.x_tot];
        res_to_pow = [res_to_pow res.to.P];
        res_ind_vel = [res_ind_vel res.v_i_to];
    catch
    end
end

%% Results

figure('Name', 'Blown wing')

colororder({'b', 'r'});

yyaxis left;
plot(res_CLmax_to, res_to_pow/1000 * unit.kW2hp); hold on;
xlabel('C_{L max} at take-off [-]')
xlim([min(res_CLmax_to), max(res_CLmax_to)])
ylabel(['Mean power required ', unit.kW2hp_txt])
ylim([min(res_to_pow/1000 * unit.kW2hp), max(res_to_pow/1000 * unit.kW2hp)])
grid on;

yyaxis right;
plot(res_CLmax_to, res_ind_vel * unit.mps2ftps); hold on;
ylabel(['Induced velocity required ', unit.mps2ftps_txt])
ylim([min(res_ind_vel * unit.mps2ftps), max(res_ind_vel * unit.mps2ftps)])
grid on;

make_fig(gcf, 'Figures/blown_wing', 0.75, 0.75, 1.5)



