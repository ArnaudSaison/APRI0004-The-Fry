function [Elem] = fulliter(Elem, Blades, Param, lossFact)
% FULLITER Solves BEMT classical formulation for propellers by iteration
%
% Calculates the inflow velocity ratio with the following assumptions:
%   - Small angles
%
% IMPORTANT NOTE:
%   - This formulation requires a non-zero axial velocity as the inflow factor is devined as (a=v_ind/V_ax -1)
%   - This solver should not be used for hovering rotors (prefer the approx_heli in tht situation).
%
% -----
%
% Synopsis: [Elem] = FULLITER(Elem, Blades, Param, lossFact)
%
% Input:    Elem     = (required) Variables computed for any individual blade element
%           Blades   = (required) Blade parameters defined in config file
%           Param    = (required) General parameters updated with flight conditions
%           lossFact = (required) Prandtl loss factor value
%
% Output:   Elem     = Updated Elem struct with inflow angles and velocities
%
% Calls:    none
%
% See also: BEMT.
%

% -----
%
% -----
% This function is one set of equations that can be used with the BEMT.
% Author: Thomas Lambert <t.lambert@uliege.be>
% ULiege - Aeroelasticity and Experimental Aerodynamics
% MIT License
% https://gitlab.uliege.be/thlamb/bemt
% -------------------------------------------------------------------------


RELAX=0.1; % Relaxation coefficient when updating inflow factors (facilitate convergence for nonlinear systems)

%% Abbreviations

tipSpeed = Blades.RADIUS*Blades.OMEGA;


%% Solve iteratively for the inflow factors

% Initialize inflow ratios using the analytical solution
tempElem = approx_heli(Elem, Blades, Param, lossFact);
a_ax = tempElem.InducedVeloc ./ Param.Air.AXIAL_VELOC;  % Take the inflow factor from the analytical solution
a_ang = zeros(1,Blades.nELEM);  % No swirl initially

% Initialize loop
inflow_converged = false;
loopCount = 0;

while ~inflow_converged && loopCount <= 500

    a_ax_OLD = a_ax;
    a_ang_OLD = a_ang;

    % Axial and angular velocities
    V_ax = (1+a_ax) * Param.Air.AXIAL_VELOC;
    V_ang = (1-a_ang) .* Blades.OMEGA .*Elem.Ypos;
    V_tot = sqrt(V_ax.^2 + V_ang.^2);

    % Angles
    phi = atan2(V_ax,V_ang);  % Inflow
    alpha = Elem.Pitch - phi;  % Angle of attack

    % Get new estimates for the coefficients
    cl = getCoeffValue('cl', Blades.Polar, alpha, Elem.Re);
    cd = getCoeffValue('cd', Blades.Polar, alpha, Elem.Re);

    % Analytical solution to the momentum and blade element equations.
    tmp = ((V_tot.^2 * Blades.nBLADES .* Elem.Chord) .* (cl.*cos(phi) - cd.*sin(phi))) ./ (8 * pi * Elem.Ypos .* Param.Air.AXIAL_VELOC.^2 .* lossFact);
    a_ax = (- 1 + sqrt(1 + 4 * tmp))/2;
    a_ang = ((V_tot.^2 * Blades.nBLADES .* Elem.Chord) .* (cl.*sin(phi) + cd.*cos(phi))) ./ (8 * pi * Elem.Ypos.^2 .* Param.Air.AXIAL_VELOC .* (1 + a_ax) .* Blades.OMEGA .* lossFact);

    % Use relaxation to faciliate convergence of nonlinear system
    a_ax = a_ax_OLD * (1-RELAX) + a_ax * RELAX;
    a_ang = a_ang_OLD * (1-RELAX) + a_ang * RELAX;

    % Convergence criterions
    axial_converged = isConverged(a_ax, a_ax_OLD, Param.Simul.CONV_CRIT);
    swirl_converged = isConverged(a_ang, a_ang_OLD, Param.Simul.CONV_CRIT);
    inflow_converged = axial_converged && swirl_converged;
    if loopCount == 500
        warning('Solver:InflowFactorsNotConverged','Solution NOT CONVERGED after reaching maximum number of iterations for inflow factors (%d).',loopCount);
    end
    loopCount = loopCount +1;
end

%% Output important results
Elem.Lambda = V_ax ./ tipSpeed;
Elem.Xi = V_ang ./ tipSpeed;

Elem.InducedVeloc = V_ax - Param.Air.AXIAL_VELOC;
Elem.TotalVeloc = sqrt(V_ax.^2 + V_ang.^2);

Elem.InflowAngle = phi;
Elem.AOA = alpha;


end
