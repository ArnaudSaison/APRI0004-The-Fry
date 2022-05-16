function [coeff] = getCoeffValue(coeffType, Polar, aoa_vect, reyn_vect)
% GETCOEFFVALUE Find the aerodynamic coefficients for a given AOA, based on the polars.
% -----
%
% Synopsis: [coeff] = GETCOEFFVALUE(coeffType, profileCoeff, Elem)
%
% Input:    coeffType    = (required) Type of coefficient to output ['cla','cl',cd']
%           Polar        = (required) Polars data
%           aoa_vect     = (required) AOA where the coefficients need to be evaluated
%           reyn_vect    = (required) Reynolds number for every element
%
% Output:   coeff = Value of the coefficients
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

% Check if enough aoas on the polars
isok_aoa = (aoa_vect(1)>Polar.alpha(1) && aoa_vect(end)<Polar.alpha(end));

if numel(Polar.reynolds) == 1
    if ~isok_aoa
        warning('Not enough angles of attacks on the polars. The range of AOA should be at least [%f - %f] deg.\nPoints for missing AOAs were given values of cl and cd corresponding to the last AOA available.',aoa_vect(1)*180/pi,aoa_vect(end)*180/pi);
    end
    % Simple interpolation
    switch coeffType
        case 'cl'
            coeff = interp1(Polar.alpha,Polar.cl',aoa_vect,'linear',Polar.cl(end));
        case 'cd'
            coeff = interp1(Polar.alpha,Polar.cd',aoa_vect,'linear',Polar.cd(end));
        otherwise
            error('Input ''coeffType''  needs to be either, ''cl'' or ''cd''. Found ''%s'' instead',coeffType)
    end

else

    % Tweak for extrapolation outside Reynolds range
    betterPolar.cl = [Polar.cl(:,1), Polar.cl, Polar.cl(:,end)];
    betterPolar.cd = [Polar.cd(:,1), Polar.cd, Polar.cd(:,end)];
    betterPolar.reynolds = [0, Polar.reynolds, 1e15];

    % Create meshgrid with the polar data points
    [X, Y] = meshgrid(Polar.alpha, betterPolar.reynolds');

    % Create meshgrid with the actual values for each element
    [Xq, Yq] = meshgrid(aoa_vect, reyn_vect');

    % Check if enough data on the polars for the interpolation
    isok_reyn = (reyn_vect(1)>betterPolar.reynolds(1) && reyn_vect(end)<betterPolar.reynolds(end));
    if ~isok_aoa || ~isok_reyn
        error('Not enough data on the Polars. The range of AOA should be at least [%f - %f] deg and the range of Reynolds should be at least [%f - %f]e5.',aoa_vect(1)*180/pi,aoa_vect(end)*180/pi,reyn_vect(1)/1e5,reyn_vect(end)/1e5);
    end

    % 2D interpolation
    switch coeffType
        case 'cl'
            Vq = interp2(X,Y,betterPolar.cl',Xq,Yq,'linear',0)';
        case 'cd'
            Vq = interp2(X,Y,betterPolar.cd',Xq,Yq,'linear',0)';
        otherwise
            error('Input ''coeffType''  needs to be either, ''cl'' or ''cd''. Found ''%s'' instead',coeffType)
    end
    coeff = diag(Vq)';
end


end
