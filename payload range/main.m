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
% User accessible codes are:
%  • 'main.m' solves the general problem
%  • 'climb_optimizer' finds the C_L for minimal power (to be changed 
%    afterwards in 'parameters.m')
%  • 'final_battery' allows to use a fixed battery mass 
%    /!\ does not check if total mass is consistent. The value must thus
%    result from a thourough analysis and be used to check that the actual
%    weight of battery chosen doesn't change other variables.
% 

clc; clear; close all;
disp('************************************************************')
disp('Program logs')
disp('************************************************************')

%% 
%==========================================================================
% Initialization of the parameters (set by user in parameters.m)
%==========================================================================
par = parameters();

%==========================================================================
% Calculating the different phases of the flight
%==========================================================================
res = performance(par);

%==========================================================================
% Estimating the energy needs and mass repartition
%==========================================================================
% Initialization of parameters
batteries_iter = 1;
res.en.mass_battery = par.mass_battery;
res.en.mass_left_after_fuel = par.mass_battery_margin;

% Iteration loop to find optimal baterries
while (batteries_iter < par.mass_battery_max_iter) && ...
      ((res.en.mass_left_after_fuel >= par.mass_battery_margin) || ...
      (res.en.mass_left_after_fuel < 0))
    
    % computing the battery
	res = batteries(res, par);
    
    % correcting the battery mass
    res.en.mass_battery = res.en.mass_battery + res.en.mass_left_after_fuel;
    
    % account for double precision so the loop won't be stuck at -1e16
    if res.en.mass_left_after_fuel < 0
        res.en.mass_battery = res.en.mass_battery - par.mass_battery_margin/20;
    end
    
    % display current state of the loop
    disp(['Battery mass at iteration ', ...
          num2str(batteries_iter), ...
          ' is ', num2str(res.en.mass_battery), ...
          ' [kg] for initial mass of ', num2str(par.mass_battery), ...
          ' [kg] and ', num2str(res.en.mass_left_after_fuel), ...
          ' [kg] still available'])
    
      % increment the loop
    batteries_iter = batteries_iter + 1;
end

% display result of the iteration
if batteries_iter < par.mass_battery_max_iter
    disp('Battery mass converged within maximum number of iterations.')
else 
    disp('Battery mass did not converge.')
end

%==========================================================================
% Printing the output
%==========================================================================
output(par, res)

