function [lossFact] = prandtlLoss(nBl, r, r0, phi, typeloss)
% PRANDTLLOSS Calculate Prandtl Loss Factor.
% -----
%
% Synopsis: [lossFact] = PRANDTLLOSS(nBl, r, phi, typeloss)
%
% Input:    nBl     = (required) Number of blades [scalar]
%           r       = (required) Non-dimensional radial position of the elements [vect]
%           r0      = (required) Non-dimensional radial position of root cutout
%           phi     = (required) Inflow angle of the elements [vect]
%           typeloss= (required) Type of loss to take into account [scalar]
%
% Output:   lossFact = Value of the loss factor
%
% Calls:    none
%
% See also: BEMT.

% -----
% This function is a mandatory part of the BEMT implementation.
% Author: Thomas Lambert <t.lambert@uliege.be>
% ULiege - Aeroelasticity and Experimental Aerodynamics
% MIT License
% https://gitlab.uliege.be/thlamb/bemt
% -------------------------------------------------------------------------

% Tip losses
f_tip = (nBl/2) * ((1-r)./(r.*sin(phi)));
F_tip = calcF(f_tip);

% Root losses
if typeloss == 2
    f_root = (nBl/2) * (r-r0)./(r.*sin(phi));
    F_root = calcF(f_root);
else
    F_root = 1;
end

% Overall loss factor
lossFact = F_tip .* F_root;

end



function [lossF] = calcF (f)
% CALCF Prandtl formula for loss factor
lossF = (2/pi) * acos(exp(-f));
end