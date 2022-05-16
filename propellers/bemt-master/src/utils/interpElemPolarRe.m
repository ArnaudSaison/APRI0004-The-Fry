function [Elem] = interpElemPolarRe(Polar, Elem)
% INTERPELEMPOLARRE Get the value of special angles (zero-lift, stall) and the lift slope based on the polars for the right Reynolds, for each blade element.
% -----
%
% Synopsis: [Elem] = CALCELEMSPECIALANGLES(Polar, reyn)
%
% Input:    Polar   = (required) Structure with the polars
%           Elem    = (required) Variables computed for any individual blade element
%
% Output:   Elem    = (required) Variables computed for any individual blade element
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


if numel(Polar.reynolds) == 1
    % Interpolate special angles for each element from the polars
    Elem.alpha_0 = Polar.alpha_0;
    Elem.alpha_stall = Polar.alpha_stall;

    % Interpolate lift curve slope for each element from the polars
    Elem.cla = Polar.cla;
else
    % Interpolate special angles for each element from the polars
    Elem.alpha_0 = interp1(Polar.reynolds, Polar.alpha_0, Elem.Re, 'linear','extrap');
    Elem.alpha_stall = interp1(Polar.reynolds, Polar.alpha_stall, Elem.Re, 'linear','extrap');

    % Interpolate lift curve slope for each element from the polars
    Elem.cla = interp1(Polar.reynolds, Polar.cla, Elem.Re, 'linear','extrap');
end


end
