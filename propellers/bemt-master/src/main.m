% MAIN Blade Element Momentum Theory wrapper for single isolated rotors.
%
% Finds the Thrust, Power and Torque of a specific rotor. This code DOES NOT
% accommodate varying planform, but takes into account linear twist and taper of
% the blades. It can also calculate the ideal twist and taper and take into
% account the blade tip losses.
%
% Three solvers are implemented (user decides which one to use, see configuration/config_template.m)
% 1. approx_heli:
%   Solves a simple system with a few assumptions (small angles, no swirl,
%   linear cl,...) as described in "Leishman, Principles of Helicopter
%   Aerodynamics, Cambridge University Press, 2006".
%
% 2. fulliter:
%   Solves a more complete set of equations iteratively, as commonly done for propellers
%   (using inflow factors). The convergence is not always guaranteed in some edge cases.
%
% 3. fullsolve:
%   Solves a more complete set of equations based on the methodology proposed
%   in "Stahlhut, C. and Leishman, J.G. Aerodyanamic Design
%   Optimization of Proprotors for Convertible-Rotor Concepts, American Helicopter
%   Society 68th Annual Forum, 1-3 May 2012".
%
% See more details on the code wiki page: https://gitlab.uliege.be/thlamb/bemt/-/wikis/home
%
% USAGE: Create a configuration file based on configurations/config_template
%        and indicate it in the firsts line of the code. Then run main.m
% -----
%
% Synopsis: [-] = MAIN(-)
%
% Input:    none
%
% Output:   none
%
% Calls:    bemt                : BEMT solver
%           checkConfig         : Validate user configuration input file
%           plotResults         : Plots main results of the solver
%           plotRotor           : Simplified 3D drawing of the blade and rotor assembly
%
% See also: BEMT, CONFIG_TEMPLATE, CHECKCONFIG, PLOTRESULTS, PLOTROTOR.

% -----
% This m-file is a wrapper for the BEMT solver.
% Author: Thomas Lambert <t.lambert@uliege.be>
% ULiege - Aeroelasticity and Experimental Aerodynamics
% MIT License
% https://gitlab.uliege.be/thlamb/bemt
% ------------------------------------------------------------------------------

clear variables; close all; clc;
addpath(genpath('.'));  % Add all functions in subfolders to matlab path


% ==============================================================================
% =                         USER INPUT REQUIRED HERE                           =
% ==============================================================================

% Specify the configuration file to be loaded here
CONFIG_FILE = 'configurations/FRY_WTP.m';

% ==============================================================================
% =                                                                            =
% ==============================================================================


%% READS CONFIG AND RUN THE BEMT SOLVER

% Check the configuration validity and load the file if no errors were raised
checkConfig(CONFIG_FILE);
run(CONFIG_FILE);

% Run solver based on CONFIG_FILE data
[Results, Elem, Param, Blades] = bemt(Param, Blades);


%% OUTPUT

% Print results in console
if Param.Simul.PRINT_CONSOLE == 1
    fprintf('Total Thrust: %0.2f N\n', Results.T)
    fprintf('Total Torque: %0.2f Nm\n', Results.Q)
    fprintf('Total Power: %0.2f kW (%0.1f hp)\n', Results.P/1000, Results.P/1000/0.746)
    disp(' ');
    fprintf('Advance ratio: %0.2f \n', Results.adv_ratio)
    fprintf('Propulsive efficiency: %0.2f %%\n', Results.eff_prop*100)
    fprintf('Rotor CT: %0.6f\n', Results.CT)
    fprintf('Rotor CQ: %0.6f\n', Results.CQ)
    disp('----------------------------------')
    disp(' ')
end

% Plot the results
if Param.Simul.SHOW_GRAPHS == 1
    plotResultsCustom(Elem);
end
if Param.Simul.SHOW_3DVIEW == 1
    plotRotor(Blades, Elem);
end

% Saves the results
if Param.Simul.SAVE_RESULTS == 1
    % Create dir if required then save all solver outputs
    if ~exist(Param.Simul.SAVE_PATH, 'dir')
        mkdir(Param.Simul.SAVE_PATH);
    end
    save([Param.Simul.SAVE_PATH, Param.Simul.SAVE_FILENAME, '.mat'], 'Results', 'Elem', 'Param', 'Blades');
end
