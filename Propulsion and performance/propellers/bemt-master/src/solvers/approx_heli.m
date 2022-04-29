function [Elem] = approx_heli(Elem, Blades, Param, lossFact)
% APPROX_HELI Bemt solver following Leishmann's Principle of Helicopter Aerodynamics.
%
% Calculates the inflow velocity ratio with the following assumptions:
%   - Small angles
%   - No swirl
%   - Linear lift coefficient
%   - cd << cl
% -----
%
% Synopsis: [Elem] = APPROX(Elem, Blades, Param, lossFact)
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
% Leishman, "Principles of Helicopter Aerodynamics", Cambridge University Press, 2006.
% -----
% This function is one set of equations that can be used with the BEMT.
% Author: Thomas Lambert <t.lambert@uliege.be>
% ULiege - Aeroelasticity and Experimental Aerodynamics
% MIT License
% https://gitlab.uliege.be/thlamb/bemt
% -------------------------------------------------------------------------

% Abbreviations
tipSpeed = Blades.RADIUS*Blades.OMEGA;
lambda_c = Param.Air.AXIAL_VELOC/tipSpeed;

cla = Elem.cla;

% Direct calculation of the inflow ratios
Elem.Lambda = sqrt(((Elem.Solidity.*cla)./(16*lossFact)-lambda_c/2).^2+(Elem.Solidity.*cla.*Elem.Pitch.*Elem.r)./(8*lossFact))-((Elem.Solidity.*cla)./(16*lossFact)-lambda_c/2);
Elem.Xi = (Blades.OMEGA*Elem.Ypos)./tipSpeed ;  % Swirl velocity is neglected -> (u_i = 0)

Elem.InducedVeloc = Elem.Lambda*tipSpeed - Param.Air.AXIAL_VELOC;
Elem.TotalVeloc = sqrt((Elem.Lambda*tipSpeed).^2+(Elem.Ypos*Blades.OMEGA).^2);

Elem.InflowAngle = atan2((Elem.Lambda*tipSpeed),(Blades.OMEGA*Elem.Ypos));
Elem.AOA = Elem.Pitch-Elem.InflowAngle;
end
