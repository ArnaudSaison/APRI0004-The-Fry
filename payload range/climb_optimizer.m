function [climb_C_L_sol] = climb_optimizer()

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


%% 
%==========================================================================
% Initialization of the parameters (set by user in parameters.m)
%==========================================================================

par = parameters();
nb = 20;

climb_C_L = linspace(par.C_L_minD, par.C_L_max, nb);
climb_rho = par.rho_0;
climb_res_power = zeros(1,nb);

%==========================================================================
% Solving the performance problem
%==========================================================================

for i = 1:1:nb
    par.C_L_climb = climb_C_L(i);
    par.rho_climb = climb_rho;
    
    res = performance_drag_study(par);
    
    climb_res_power(i) = res.climb.P_engine;
end


% Solution
[~, sol.pos] = min(climb_res_power);
climb_C_L_sol = climb_C_L(sol.pos);
end


