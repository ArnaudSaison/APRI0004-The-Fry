function [Rotor, Elem, Param, Blades] = bemt(Param, Blades)
% BEMT Blade Element Momentum Theory solver for single isolated rotors.
%
% Finds the Thrust, Power and Torque of a specific rotor. This code DOES NOT
% accommodate varying planform, but takes into account linear twist and taper of
% the blades. It can also calculate the ideal twist and taper and take into
% account the blade tip (and hub) losses.
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
% -----
%
% Synopsis: [Rotor, Elem, Param, Blades] = BEMT(Param, Blades)
%
% Input:    Param  = (required) General parameters defined in config file
%           Blades = (required) Blade parameters defined in config file
%
% Output:   Rotor   = Results for the rotor computed with BEMT
%           Elem    = Variables computed for any individual blade element
%           Param   = General parameters updated with flight conditions
%           Blades  = Blade parameters updated with unit conversion and spanwise variation of forces, inflow and angles
%
% Calls:    approx_heli         : Solver for apporximate equation (Leishman)
%           fulliter            : Solver for full system equations based on iteration (UNSTABLE)
%           fullsolve           : Solver for full system equations (Stahlhut)
%           calcPolarChar       : Calculates special angles (zero-lift, stall) and slope of polar curves
%           interpElemPolarRe   : Interpolate polar data for each element (based on Reynolds)
%           prandtlLoss         : Calculate losses with Prandtl formula
%           getCoeffValue       : Get cl and cd from airfoil data
%           getCoeffAngle       : Get specific angles from cl polar of the airfoil
%
%
% See also: CONFIG_TEMPLATE, APPROX_HELI, FULLSYST, PRANDTLLOSS, GETCOEFFVALUE, GETCOEFFANGLE.

% -----
% This function is the main solver of the BEMT implementation.
% Author: Thomas Lambert <t.lambert@uliege.be>
% ULiege - Aeroelasticity and Experimental Aerodynamics
% MIT License
% https://gitlab.uliege.be/thlamb/bemt
% -------------------------------------------------------------------------


%% INITIALIZATION

% Convert angles to radians and RPM to rad/s
Blades.COLL_PITCH = Blades.COLL_PITCHdeg * (pi/180);
Blades.TWIST_MAX = Blades.TWIST_MAXdeg * (pi/180);
Blades.THETA_TIP = Blades.THETA_TIPdeg * (pi/180);
Blades.TWIST = Blades.TWISTdeg * (pi/180);
Blades.OMEGA = Blades.OMEGArpm * (2*pi/60);

% Get atmospheric conditions
if exist('atmosisa','file') == 2
    [Param.Air.temperature, Param.Air.vSound, Param.Air.pressure, Param.Air.rho] = atmosisa(Param.Air.ALTITUDE);
else
    error('Function ''%s'' is missing. You can download an alternative version from: https://github.com/lentzi90/octave-atmosisa\nPlace it this project folder and restart this computation.', 'atmosisa');
end

% Calculate air kinematic viscosity (used for Reynolds) with Sutherland model (https://www.cfd-online.com/Wiki/Sutherland's_law)
T0 = 273.15;  % Reference temperature, [K]
MU0 = 1.716e-5;  % Reference viscosity, [kg/ms]
SUTH = 110.4;  % Sutherland temperature, [K]
Param.Air.mu = MU0 * (Param.Air.temperature/T0).^(3/2) .* (T0 + SUTH) ./ (Param.Air.temperature + SUTH);  % Dynamic viscosity
Param.Air.nu = Param.Air.mu / Param.Air.rho;  % Kinematic viscosity

% Load airfoil polars from file
tmpPolar = load(['./airfoils_data/', Blades.COEFF_FILE]);
tmpField = fieldnames(tmpPolar);
Blades.Polar = tmpPolar.(tmpField{1});
Blades.Polar = calcPolarChar(Blades.Polar);
clear tmpPolar tmpField


%% GEOMETRY

% Rotor parameters
tipSpeed = Blades.RADIUS * Blades.OMEGA;

% Blade discretization
dy = (Blades.RADIUS - Blades.ROOT_CUTOUT) / Blades.nELEM;  % Span of one blade element
Elem.Ypos = Blades.ROOT_CUTOUT + dy/2:dy:Blades.RADIUS-dy/2;  % Absolute position of the blade elements
r = Elem.Ypos ./ Blades.RADIUS;  % Relative position of the blade elements
Elem.r = r;
Elem.r0 = Blades.ROOT_CUTOUT / Blades.RADIUS;  % Relative position of the root cutout
Elem.dy = dy;

% Geometric twist of each element
if Param.Simul.TWIST == 0  % Linear twist
    Elem.Twist = Blades.THETA_TIP + (r-1) * Blades.TWIST;
    
elseif Param.Simul.TWIST == 1  % Ideal twist
    Elem.Twist = Blades.THETA_TIP ./ r;
    if Param.Simul.TRIM_IDEAL == 1  % Trim excessive values
        ind = find(Elem.Twist > Elem.Twist(end) + Blades.TWIST_MAX);
        if ~isempty(ind)
            Elem.Twist(ind) = Elem.Twist(ind(end)+1) * ones(1,length(ind));
        end
    end
    
elseif Param.Simul.TWIST == 2  % Custom twist
    Elem.Twist = Blades.CUSTOM_TWIST * (pi/180);
end

Elem.geomPitch = Blades.COLL_PITCH + Elem.Twist;  % Geometric pitch angle (collective pitch + geometric twist)

% Chord of each blade element (taper)
if Param.Simul.TAPER == 0  % Linear taper
    Elem.Chord = Blades.ROOT_CHORD - (Blades.ROOT_CHORD - Blades.TIP_CHORD) * r;
    
elseif Param.Simul.TAPER == 1  % Ideal taper
    Elem.Chord = Blades.TIP_CHORD./r;
    if Param.Simul.TRIM_IDEAL == 1  % Trim excessive values
        ind = find(Elem.Chord > 1/Blades.TAPER_MAX * Elem.Chord(end));
        if ~isempty(ind)
            Elem.Chord(ind) = Elem.Chord(ind(end)+1) * ones(1,length(ind));
        end
    end
    
elseif Param.Simul.TAPER == 2  % Custom taper
    Elem.Chord = Blades.CUSTOM_CHORD;
end

% Solidity
Elem.Solidity = (Blades.nBLADES * Elem.Chord) / (pi * Blades.RADIUS);
rotorSolidity = Blades.nBLADES * sum(Elem.Chord * dy) / (pi * Blades.RADIUS.^2);

% Reynolds number estimate for each element
Elem.Re = sqrt((Blades.OMEGA .* Elem.Ypos).^2 + Param.Air.AXIAL_VELOC^2) .* Elem.Chord / Param.Air.nu;

% Calculate the zero-lift aoa, stall aoa and lift curve slope for each element
Elem = interpElemPolarRe(Blades.Polar, Elem);
% Define the element true pitch based on zero-lift line
Elem.Pitch = Elem.geomPitch - Elem.alpha_0;


%% INFLOW

% Calculate inflow velocities and angles
if strcmpi(Param.Simul.SYSTEM_SOLVER, 'approx_heli')  % Simple solver
    % No swirl, small angles, linear cl
    
    % Calculate first estimation of inflow without losses
    Elem = approx_heli(Elem, Blades, Param, 1);
    
    % Include tip (and hub) losses
    if Param.Simul.LOSSES == 1 || Param.Simul.LOSSES == 2
        
        % Iterate until inflow ratio is converged
        lambda_converged = false;
        loopCount = 1;
        
        while ~lambda_converged && loopCount <= 500
            lambda_old = Elem.Lambda;
            
            % Calculate loss factor and update inflow
            lossFact = prandtlLoss(Blades.nBLADES, Elem.r, Elem.r0, Elem.InflowAngle,Param.Simul.LOSSES);
            Elem = approx_heli(Elem, Blades, Param, lossFact);
            
            % Convergence criterions
            lambda_converged = isConverged(Elem.Lambda, lambda_old, Param.Simul.CONV_CRIT);
            if loopCount == 500
                warning('Solver:PrandtlNotConverged','Solution NOT CONVERGED after reaching maximum number of iterations for losses calculation (%d).',loopCount);
            end
            loopCount = loopCount +1;
        end
        
    end
    
elseif strcmpi(Param.Simul.SYSTEM_SOLVER, 'fulliter')  % Complete system of equation by iteration (UNSTABLE)
    % No particular aditional assumption
    
    % Calculate first estimation of inflow without losses
    Elem = fulliter(Elem, Blades, Param, 1);
    
    % Include tip losses
    if Param.Simul.LOSSES == 1 || Param.Simul.LOSSES == 2
        
        % Iterate until inflow ratio is converged
        lambda_converged = false;
        loopCount = 1;
        
        while ~lambda_converged && loopCount <= 500
            lambda_old = Elem.Lambda;
            
            % Include tip (and hub) losses
            lossFact = prandtlLoss(Blades.nBLADES, Elem.r, Elem.r0, Elem.InflowAngle, Param.Simul.LOSSES);
            Elem = fulliter(Elem, Blades, Param, lossFact);
            
            % Convergence criterions
            lambda_converged = isConverged(Elem.Lambda, lambda_old, Param.Simul.CONV_CRIT);
            if loopCount == 500
                warning('Solver:PrandtlNotConverged','Solution NOT CONVERGED after reaching maximum number of iterations for losses calculation (%d).',loopCount);
            end
            loopCount = loopCount +1;
        end
        
    end
    
elseif strcmpi(Param.Simul.SYSTEM_SOLVER, 'fullsolve')   % Solve full system of equation
    % Solve full system of equation (Include swirl, large angles and tip (and hub) losses)
    Elem = fullsolve(Elem, Blades, Param);
end


% Check if stalled parts or non lifting parts
if any(Elem.AOA > Elem.alpha_stall)
    warning('Inflow:Stall', 'Possible stall: Some sections have an higher AOA than the profile stall angle.')
end
if any(Elem.AOA < Elem.alpha_0)
    warning('Inflow:NonLift','Non lifting parts: Some sections have an AOA smaller than the zero-lift angle.')
end


%% AERODYNAMIC COEFFICIENTS AND FORCES

% Aerodynamic coefficient of each blade element
Elem.cl = getCoeffValue('cl', Blades.Polar, Elem.AOA, Elem.Re);
Elem.cd = getCoeffValue('cd', Blades.Polar, Elem.AOA, Elem.Re);

% Lift and drag of each blade element
Elem.dL = 0.5 * Param.Air.rho * Elem.TotalVeloc.^2 .* Elem.Chord .* Elem.cl * dy;
Elem.dD = 0.5 * Param.Air.rho * Elem.TotalVeloc.^2 .* Elem.Chord .* Elem.cd * dy;

% Thrust and Torque of each blade element
Elem.dT = Blades.nBLADES * (Elem.dL .* cos(Elem.InflowAngle) - Elem.dD .* sin(Elem.InflowAngle));
Elem.dQ = Blades.nBLADES * ((Elem.dL .* sin(Elem.InflowAngle) + Elem.dD .* cos(Elem.InflowAngle)) .* Elem.Ypos);
Elem.dCT = Elem.dT / (Param.Air.rho * pi * (Blades.RADIUS)^2 * (tipSpeed)^2);
Elem.dCQ = Elem.dQ / (Param.Air.rho * pi * (Blades.RADIUS)^2 * (tipSpeed)^2 * Blades.RADIUS);

% Induced, profile and total power requirements of each blade element
Elem.dP_i = Blades.nBLADES * (Elem.dL .* sin(Elem.InflowAngle) .* Elem.Ypos * Blades.OMEGA);
Elem.dP_p = Blades.nBLADES * (Elem.dD .* cos(Elem.InflowAngle) .* Elem.Ypos * Blades.OMEGA);
Elem.dP = Elem.dP_i + Elem.dP_p;
Elem.dCP = Elem.dP / (Param.Air.rho * pi * (Blades.RADIUS)^2 * (tipSpeed)^3);

% Forces and coefficients for the whole rotor
Rotor.T = sum(Elem.dT);
Rotor.Q = sum(Elem.dQ);
Rotor.P = sum(Elem.dP);
Rotor.CT = sum(Elem.dCT);
Rotor.CQ = sum(Elem.dCQ);
Rotor.CP = sum(Elem.dCP);

% Other useful outputs
Rotor.adv_ratio = Param.Air.AXIAL_VELOC / tipSpeed; % Advance Ratio
Rotor.eff_prop = Rotor.T * Param.Air.AXIAL_VELOC / Rotor.P;  % Rotor propulsive efficiency
Rotor.diskLoading = Rotor.CT * Param.Air.rho * (tipSpeed)^2;  % Disk loading
Rotor.CT_Sigma = Rotor.CT / rotorSolidity;  % Thrust coefficient over solidity

Rotor.FM = 1/sqrt(2) * Rotor.CT^(3/2) / Rotor.CP;  % Figure of Merit


end
