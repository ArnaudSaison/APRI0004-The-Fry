function checkConfig(config_file)
% CHECKCONFIG Check configuration file validity.
% -----
%
% Synopsis: CHECKCONFIG(config_file)
%
% Input:    config_file = (required) configuration file to be checked
%
% Output:   No output if everything is fine.
%           Error message with detail about the issue otherwise.
%
% Calls:    none
%
% See also: BEMT.


% ------
% This function is an optional part of the BEMT code.
% Author: Thomas Lambert <t.lambert@uliege.be>
% ULiege - Aeroelasticity and Experimental Aerodynamics
% MIT License
% https://gitlab.uliege.be/thlamb/bemt
% -------------------------------------------------------------------------


% Load config file
run(config_file);

% Warning if running the default config file
if strcmpi(config_file, 'configurations/config_template.m')
    warning('Config:useTemplate', 'You are currently running the template configuration.\nCreate a new config for your system based on config_template and indicate it in main.m (l.54).')
    disp(' ')
end

% Check binary variables (all following var must be either 0 or 1)
validateattributes(Param.Simul.SAVE_RESULTS, {'double'},{'binary'}, '', 'Param.Simul.SAVE_RESULTS')
validateattributes(Param.Simul.SHOW_GRAPHS, {'double'},{'binary'}, '', 'Param.Simul.SHOW_GRAPHS')
validateattributes(Param.Simul.SHOW_3DVIEW, {'double'},{'binary'}, '', 'Param.Simul.SHOW_3DVIEW')
validateattributes(Param.Simul.PRINT_CONSOLE, {'double'},{'binary'}, '', 'Param.Simul.PRINT_CONSOLE')
validateattributes(Param.Simul.TRIM_IDEAL, {'double'},{'binary'}, '', 'Param.Simul.TRIM_IDEAL')

% Check if int variables are allowed (all following var must be either 0, 1 or 2)
validateattributes(Param.Simul.TWIST, {'double'},{'scalar','integer','>=',0,'<=',2}, '', 'Param.Simul.TWIST')
validateattributes(Param.Simul.TAPER, {'double'},{'scalar','integer','>=',0,'<=',2}, '', 'Param.Simul.TAPER')
validateattributes(Param.Simul.LOSSES, {'double'},{'scalar','integer','>=',0,'<=',2}, '', 'Param.Simul.LOSSES')

% Check if correct number element for custom twist and taper (should be equal to number of blade elements)
if Param.Simul.TWIST == 2
    validateattributes(Blades.CUSTOM_TWIST, {'double'},{'numel',Blades.nELEM}, '', 'Elem.CUSTOM_TWIST')
end
if Param.Simul.TAPER == 2
    validateattributes(Blades.CUSTOM_CHORD, {'double'},{'numel',Blades.nELEM}, '', 'Elem.CUSTOM_CHORD')
end

% Check tip speed and raise warnings accordingly (approx a_sound = 343)
tipMach = Blades.RADIUS * Blades.OMEGArpm * (2*pi/60) / 343;
if  (tipMach > 0.75) && (tipMach < 1)
    warning('Config:transonicTip', 'Transonic flow at the tip (Mach number is %0.2f).', tipMach)
    disp(' ')
end
if  (tipMach >= 1)
    warning('Config:supersonicTip', 'SUPERSONIC flow at the tip (Mach number is %0.2f).', tipMach)
    disp(' ')
end

% Check if solver name is correct
validatestring(Param.Simul.SYSTEM_SOLVER, ["approx_heli", "fulliter", "fullsolve"], '', 'Param.Simul.SYSTEM_SOLVER');

% Check if there is an axial speed with fulliter
if strcmpi(Param.Simul.SYSTEM_SOLVER,"fulliter") &&  Param.Air.AXIAL_VELOC == 0
    error("fulliter solver requires an AXIAL_VELOC > 0");
end

end