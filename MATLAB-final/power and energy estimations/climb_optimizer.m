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

clc; clear; close all;

%% 
%==========================================================================
% Initialization of the parameters (set by user in parameters.m)
%==========================================================================

par = parameters();

climb_C_L = linspace(par.C_L_minD, par.C_L_max, 50);
climb_rho = [par.rho_0, par.rho_cruise];
climb_res_TAS = zeros(2, 50);
climb_res_power = zeros(2, 50);
climb_res_angle = zeros(2, 50);
climb_res_time = zeros(2, 50);

%==========================================================================
% Solving the performance problem
%==========================================================================
for j = 1:2
for i = 1:1:50
    par.C_L_climb = climb_C_L(i);
    par.rho_climb = climb_rho(j);
    
    res = performance(par);
    
    climb_res_TAS(j, i) = res.climb.V_TAS;
    climb_res_power(j, i) = res.climb.P_engine;
    climb_res_angle(j, i) = res.climb.angle;
    climb_res_time(j, i) = res.climb.time;
end
end

% Solution
[~, sol.pos] = min(climb_res_power(1, :));
sol.climb_C_L = climb_C_L(sol.pos);
sol.TAS = climb_res_TAS(1, sol.pos);
sol.power = climb_res_power(1, sol.pos);
sol.angle = climb_res_angle(1, sol.pos);
sol.time = climb_res_time(1, sol.pos);

[~, sol.pos2] = min(climb_res_power(2, :));
sol.climb_C_L2 = climb_C_L(sol.pos2);

%==========================================================================
% Figures
%==========================================================================
%
fig.Power = figure('Name', 'Power');

plot(climb_C_L, climb_res_power(1, :) / 1e3, '-b', 'DisplayName', 'takeoff altitude'); hold on;
plot(climb_C_L, climb_res_power(2, :) / 1e3, '-r', 'DisplayName', 'cruise altitude');
xlabel('C_{L climb}');
grid('on');
ylabel('Power [kW]');
legend();

make_fig(fig.Power, 'Figures/climb_Power', 1, 1, 1.5);

%
fig.Power_hp = figure('Name', 'Power hp');

plot(climb_C_L, climb_res_power(1, :) / 746, '-b', 'DisplayName', 'takeoff altitude'); hold on;
plot(climb_C_L, climb_res_power(2, :) / 746, '-r', 'DisplayName', 'cruise altitude');
xlabel('C_{L climb}');
grid('on');
ylabel('Power [hp]');
xticks([par.C_L_minD, sol.climb_C_L, par.C_L_max]);
legend();

make_fig(fig.Power_hp, 'Figures/climb_Power_hp', 1, 1, 1.5);

%
fig.V_TAS = figure('Name', 'V_TAS');

plot(climb_C_L, climb_res_TAS(1, :), '-b', 'DisplayName', 'takeoff altitude'); hold on;
plot(climb_C_L, climb_res_TAS(2, :), '-r', 'DisplayName', 'cruise altitude');
xlabel('C_{L climb}');
grid('on');
ylabel('V_{TAS} [m/s]');
legend();

make_fig(fig.V_TAS, 'Figures/climb_V_TAS', 1, 1, 1.5);

%
fig.angle = figure('Name', 'angle');

plot(climb_C_L, climb_res_angle(1, :), '-b', 'DisplayName', 'takeoff altitude'); hold on;
plot(climb_C_L, climb_res_angle(2, :), '-r', 'DisplayName', 'cruise altitude');
xlabel('C_{L climb}');
grid('on');
ylabel('angle [°]');
legend();

make_fig(fig.angle, 'Figures/climb_angle', 1, 1, 1.5);

%==========================================================================
% Command window
%==========================================================================
disp(['<strong>Power of the engine = </strong>', ...
num2str(sol.power / 1000), ' [kW]  /  ', ...
num2str(sol.power / 746), ' [imp. hp]']);

disp(['<strong>TAS during climb = </strong>', ...
num2str(sol.TAS), ' [m/s]']);

disp(['<strong>Climb angle = </strong>', ...
num2str(sol.angle), ' [°]']);

disp(['<strong>Climb C_L = </strong>', ...
num2str(sol.climb_C_L)]);

disp(['<strong>Climb C_L2 = </strong>', ...
num2str(sol.climb_C_L2)]);

disp(['<strong>Time to climb = </strong>', ...
num2str(sol.time/60), ' [min]']);




