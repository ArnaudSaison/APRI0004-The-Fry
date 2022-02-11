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
% Initialising flight phases at maximum range
%==========================================================================
res = performance(par);
res = batteries(res, par);

%==========================================================================
% Estimating the energy needs and mass repartition
%==========================================================================

% set of ranges to find the max payload for
% upper limit is a guess that must be adjusted according to the results in
% order not to either cut the graph or lose computation time solving too
% far
% below the design range, the payload is the max payload
p_r.nb_ranges = 50;
p_r.ranges = linspace(par.range, par.range*1.5, p_r.nb_ranges)

p_r.total_masses = zeros(1, p_r.nb_ranges)
p_r.fuel_masses = zeros(1, p_r.nb_ranges)
p_r.payload_masses = zeros(1, p_r.nb_ranges)

p_r.iter = 1:1:p_r.nb_ranges;

% loop to test every range
for p_r_i = p_r.iter
    
    % setting the current range
    par.range = p_r.ranges(p_r_i);

    % Initialization of parameters
    iter.i = 1;
    iter.max = 20;
    iter.mass_margin = 1; % [kg]
    iter.max_ICE_power = 200 * 746; % [W]
    iter.max_fuel = 75; % [l]

    iter.

    % imposed range
    % we are looking for takeoff mass
    % we start looking from 100 kg below the last value 

    % conditions 
    % 1) < max itaterations for the loop not to get stuck -> is notified as
    %    error
    % 2) takeoff mass is consistant with total mass
    % 3) ICE is powerful enough
    % 4) above 0 payload
    % 5) fuel below max vollume
    %
    % at each iteration, max takeoff mass is changed to reflect changes on
    % objective: maximize total mass
    
    iter.last_res = res;
    iter.last_par = par;
    
    while (iter.i < iter.max) && ...
          ((res.en.mass_left_after_fuel >= par.mass_battery_margin) || ...
          (res.en.mass_left_after_fuel < 0))

        % to res-estimate power requirements
        res = performance(par);

        % to res-estimate fuel mass
        res = batteries(res, par);

        % 
    end




    % Iteration loop to find optimal baterries
    while (iter.i < par.mass_battery_max_iter) && ...
          ((res.en.mass_left_after_fuel >= par.mass_battery_margin) || ...
          (res.en.mass_left_after_fuel < 0))

        res = performance(par);

        % computing the battery
        res = batteries(res, par);

        % correcting the battery mass
        res.en.mass_battery = res.en.mass_battery + res.en.mass_left_after_fuel;

        % account for double precision so the loop won't be stuck at -1e16
        if res.en.mass_left_after_fuel < 0
            res.en.mass_battery = res.en.mass_battery - par.mass_battery_margin/20;
        end

        % increment the loop
        iter.i = iter.i + 1;
    end

    % display result of the iteration
    if iter.i < par.mass_battery_max_iter
        disp('Battery mass converged within maximum number of iterations.')
    else 
        disp('Battery mass did not converge.')
    end
end
%==========================================================================
% Printing the output
%==========================================================================







