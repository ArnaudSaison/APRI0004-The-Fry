% 

%% Initialization
close all; clear; clc;
addpath(genpath('.'));  % Add all functions in subfolders to matlab path

%% Base parameters
% Specify the configuration file to be loaded here
CONFIG_FILE = 'configurations/config_custom.m';

checkConfig(CONFIG_FILE);
run(CONFIG_FILE);

%% Additional parameters
Param.Add.VEL_RANGE = linspace(100,170,5) /3.6*1.852;  % [kts]->[m/s]   velocities range
Param.Add.DRAG = 1805;          % [N]   drag during level flight = thrust
Param.Add.DRAG_MAX_I = 20;      % [#]   max iterations
Param.Add.DRAG_TOL = 1e-3;      % [N]   tolerance
Param.Add.INIT_COLL = [0,20];   % [N]   tolerance


%% Finding optimal collective pitches for range of velocities
Add.eta_P = []

for vel = Param.Add.VEL_RANGE
    
    % changing the velocity
    Param.Air.AXIAL_VELOC = vel;

    % starting values of the bisection method
    iter.a = Param.Add.INIT_COLL(1) /180*pi;
    iter.b = Param.Add.INIT_COLL(2) /180*pi;

    % iterator
    i = 1;

    while Param.Add.DRAG_MAX_I >= i && (iter.b - iter.a) > Param.Add.DRAG_TOL
        % using middle value
        Blades.COLL_PITCHdeg = (iter.a + iter.b) / 2;
        
        try
            [Results, Elem, Param, Blades] = bemt(Param, Blades);
            
            % choosing new boundaries
            if Param.Add.DRAG < Results.T
                iter.b = Blades.COLL_PITCHdeg;

            else % greater than searched value
                iter.a = Blades.COLL_PITCHdeg;

            end
        catch
            iter.a = iter.a - 5;
        end
        
        % iteration increment
        i = i+1;
    end

    % warning if no convergence
    if Param.Add.DRAG_MAX_I >= i
        disp(['converged after ', num2str(i), ' iterations']);
        
        Add.eta_P = [Add.eta_P Results.eff_prop];
    else
        warning('did not converge')
        
        Add.eta_P = 0;
    end
    clear i;
end

%plot(Param.Add.VEL_RANGE, Results.eff_prop)



