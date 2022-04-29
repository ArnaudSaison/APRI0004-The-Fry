function [converged] = isConverged(newArray, oldArray, convCrit)
% ISCONVERGED Determine convergence status of an array
% -----
%
% Synopsis: [converged] = ISCONVERGED(newArray, oldArray, convCrit)
%
% Input:    newArray = (required) New values of the Array
%           oldArray = (required) Previous value of the Array
%           convCrit = (required) Convergence criterion (scalar)
%
% Output:   converged = True/False (bool).
%
% Calls:    none

% ------
% Author: Thomas Lambert <t.lambert@uliege.be>
% ULiege - Aeroelasticity and Experimental Aerodynamics
% MIT License
% https://gitlab.uliege.be/thlamb/bemt
% -------------------------------------------------------------------------


% Relative difference between new and old values for every element of array
converged = all(abs(newArray-oldArray)./oldArray < convCrit);

end
