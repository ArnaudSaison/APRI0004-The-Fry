% CONFIG_CARADO Configuration file for BEMT solver based on Caradonna and Tung rotor.
%
% This file specifies the blade gemoetry and kinematics, the simulation parameters
% and the freestream constants.
% -----
% Synopsis: [-] = CONFIG_CARADO(-)
%
% Input:    none
%
% Output:   N/A
%
% Calls:    none
%
% See also: BEMT.


% -----
% Caradonna and Tung, Experimental and Analytical Studies of a Model Helicopter Rotor in Hover. 1981. [NASA TM-81232].
%
% -----
% Example config file for the Caradonna and Tung rotor.
% Author: Thomas Lambert <t.lambert@uliege.be>
% ULiege - Aeroelasticity and Experimental Aerodynamics
% MIT License
% https://gitlab.uliege.be/thlamb/bemt
% -------------------------------------------------------------------------


% -------- Simulation options --------
Param.Simul.SAVE_RESULTS    = 0;    % Save results in mat files
Param.Simul.SAVE_PATH       = '../validation_data/';  % Folder where to save the results (put './' for current folder)
Param.Simul.SAVE_FILENAME   = 'carado';    % Name of the saved results file
Param.Simul.SHOW_GRAPHS     = 1;    % Show all graphs (forces, angles, speed, ...)
Param.Simul.SHOW_3DVIEW     = 1;    % Show the 3D view of the whole rotor
Param.Simul.PRINT_CONSOLE   = 1;    % Print the final results in console

Param.Simul.SYSTEM_SOLVER = 'approx_heli';  % System of equations to solve ('approx_heli', 'fulliter', 'fullsolve')

Param.Simul.LOSSES      = 1;    % Include losses using Prandtl formula (0 = Ignore, 1 = Tip only, 2 = Tip and Hub)
Param.Simul.TWIST       = 0;    % Twist (0 = Linear or None, 1 = Ideal, 2 = Custom)
Param.Simul.TAPER       = 0;    % Taper (0 = Linear or Constant, 1 = Ideal, 2 = Custom)
Param.Simul.TRIM_IDEAL  = 1;    % Trim IDEAL twist and taper to realistic bound values set in Blade.TAPER_MAX and Blade.TWIST_MAX (0 = No trim, 1 = Trim)
Param.Simul.CONV_CRIT   = 1e-4; % Convergence criterion value


% ------------ Freestream ------------
Param.Air.ALTITUDE      = 0;                % Altitude, [m]
Param.Air.AXIAL_VELOC   = 10;               % Axial velocity of freestream (i.e. advance speed for propeller, climb velocity for helicopter), [m/s]
Param.Air.DYN_VISCOSITY = 1.7332.*10.^-5;   % Dynamic viscosity of air (at 0 C), [Pa s]


% -------------- Blades --------------
Blades.PROFILE_FILE = 'naca0012.dat';       % Profile coordinates for 3D view of the final blade
Blades.COEFF_FILE   = 'naca0012polar.mat';  % Lift and drag polars to use for cl and cd calculation

Blades.OMEGArpm = 1250;     % Rotational speed, [RPM]
Blades.nELEM    = 100;      % Number of blade elements

Blades.nBLADES          = 2 ;       % Number of blades on the rotor
Blades.RADIUS           = 1.143 ;   % Radius (from HUB center to tip), [m]
Blades.ROOT_CHORD       = 0.1905 ;  % Root chord (will only be used if linear/constant taper), [m]
Blades.TIP_CHORD        = 0.1905 ;  % Tip chord, [m]
Blades.ROOT_CUTOUT      = 0.1905 ;  % Root cutout, [m]`
Blades.THETA_TIPdeg     = 0;        % Twist angle at tip, [deg]
Blades.COLL_PITCHdeg    = 5;        % Collective pitch, [deg]


% -------- Optional parameters --------
% If Param.Simul.TRIM_IDEAL = 1 (trim idal twist or taper to avoid infinite values near root)
Blades.TAPER_MAX = 0.5;     % Maximal taper ratio (to avoid infinite taper at root), [-]
Blades.TWIST_MAXdeg = 40;   % Maximal twist (to avoid infinite twist at root), [deg]

% If Param.Simul.TWIST = 1 (Linear or constant twist)
Blades.TWISTdeg = -20;      % Blade twist (0 if untwisted blades) [deg]

% If Param.Simul.TWIST = 2 (Custom twist)
Blades.CUSTOM_TWIST = []; % Twist of each individual element [deg] (! vector of [Blades.nELEM x 1])

% If Param.Simul.TAPER = 2 (Custom taper)
Blades.CUSTOM_CHORD = []; % Chord of each individual element [m] (! vector of [Blades.nELEM x 1])
