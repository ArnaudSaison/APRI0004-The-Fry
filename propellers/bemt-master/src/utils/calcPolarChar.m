function [Polar] = calcPolarChar(Polar)
%  CALCPOLARCHAR Calculate Polar Characteristic values (zero-lift angle, stall angle, lift slope).
% -----
%
% Synopsis: [Polar] = CALCPOLARCHAR(Polar)
%
% Input:    Polar   = (required) Structure with the polars
%
% Output:   Polar   = Updated structure with extra fields for the angles and the lift curve slope
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

% Abbreviations
nPolars = length(Polar.reynolds); % Number of polars

% Initialize
Polar.alpha_0 = zeros(size(Polar.reynolds));
Polar.alpha_stall = zeros(size(Polar.reynolds));

% Calculate angles for each polar in the input structure
for iPol = 1:nPolars
    
    % Zero-lift angles
    ind_firstClpos = find(Polar.cl(2:end,iPol)>0,1);  % Find first cl > 0
    if isempty(ind_firstClpos) % If no zero lift angle, throw an error
        error('Not enough data on polars (impossible to determine zero lift angle). Please regenerate the polar with a greater range of AOA.');
    end
    if ind_firstClpos == 1
        ind_firstClpos = 2;
    end
    cl0(:,1) = Polar.alpha(ind_firstClpos-1:ind_firstClpos); % AOA of the two points around cl = 0
    cl0(:,2) = Polar.cl(ind_firstClpos-1:ind_firstClpos,iPol); % cl of the two points around cl = 0
    Polar.alpha_0(iPol) = cl0(1,1) - cl0(1,2) * (cl0(2,1)-cl0(1,1))/(cl0(2,2)-cl0(1,2));  % Line equation to find the angle
    
    % Stall angle (not perfect, but only used to warn users)
    maxFound = false;
    minWidth = 0;
    while ~maxFound % Iterate in order to find only one peak (largest one should indicate the stall)
        [~,locs] = findpeaks(Polar.cl(:,iPol),'MinPeakWidth',minWidth);
        if numel(locs)==1
            maxFound = true;
        else
            minWidth = minWidth + 1;
        end
    end
    Polar.alpha_stall(iPol) = Polar.alpha(locs);
    
    
    % Lift curve linear range
    linRangeInd = linearRange(Polar.alpha, Polar.cl(:,iPol));
    Polar.linRangeAOA(:,iPol) = Polar.alpha(linRangeInd);
    
    % Lift curve slope
    p = polyfit(Polar.alpha(linRangeInd),Polar.cl(linRangeInd,iPol),1); % First-order plynomial fit
    Polar.cla(iPol)=p(1); % Slope
end

end


function [linRangeInd] = linearRange(aoa, cl)
% LINEARRANGE Finds the linear range of the lift curve
% -----
%
% Synopsis: [Polar] = LINEARRANGE(aoa, cl)
%
% Input:    aoa   = (required) angle of attacks
%           cl    = (required) lift coefficient
%
% Output:   linRangeInd = range of indexes between which Cl curve is (quasi) linear
%
% Calls:    none
%
% See also: BEMT, CALCPOLARCHAR.

%% CONSTANTS
% May require some user tweaking depending on input Polar (low Re, very high AOA, etc)

cla_ideal = 2*pi;   % Ideal Cl curve, [rad]
dclTol = 0.50;      % Tolerance on the slope (dCl should be at least  > dclTol*cla_ideal in order to be counted in the linear zone)
ddClTol = 0.01;     % Tolerance on the curvature (ddClTol > ddCl > -ddClTol to be considered linear)

forceLinResults = 1; % Throws an error if valid linear range bounds were not founds (0 to deactivate)

%% CALCULATION OF LINEAR RANGE

% First and second derivatives
daoa = diff(aoa);
dCl = diff(cl)./daoa(1);
ddCl = diff(cl,2)./daoa(1)^2;

% Adapt tolerance to angle of attack
ddClTol = ddClTol/daoa(1)^2;

% Initialize loop
maxFound = false;
minFound = false;

% Initialize linear range start and end (useful for debugging)
aoaMinInd = 1;
aoaMaxInd = length(aoa);

for i = 1:(length(aoa))-2
    % Linear range is when
    %   1. The slope is larger than the tolerance on cla_ideal
    %   2. The curvature is almost zero (i.e. between its tolerance bounds).
    
    % Beginning of linear range
    if (dCl(i) > dclTol*cla_ideal) && ((ddCl(i) < ddClTol) && (ddCl(i) > -ddClTol)) && ~minFound
        minFound = true;
        aoaMinInd = i;
    end
    
    % End of the linear range (same criterion as before, but starts from the other end of the array)
    if (dCl(end+1-i) > dclTol*cla_ideal) && ((ddCl(end+1-i) < ddClTol) && (ddCl(end+1-i) > -ddClTol)) && ~maxFound
        maxFound = true;
        aoaMaxInd = length(dCl) - i;
    end
    
    % Stop when both min and max were found
    if minFound && maxFound
        break;
    end
end

% Ensure valid linear range was found
if ~minFound || ~maxFound
    if forceLinResults == 1
        % Display issue
        figure('Name','Debug:LinearRange Polar')
        hold on
        plot(aoa,cl,'--');
        plot(aoa(aoaMinInd:aoaMaxInd),cl(aoaMinInd:aoaMaxInd), 'color', 'r', 'lineWidth', 1.2);
        if ~minFound
            plot(aoa(aoaMinInd),cl(aoaMinInd),'s','MarkerFaceColor','b','Markersize',12)
        end
        if ~maxFound
            plot(aoa(aoaMaxInd),cl(aoaMaxInd),'s','MarkerFaceColor','b','Markersize',12)
        end
        grid on
        legend('Polar','Linear range','Missing point(s)')
        
        % Throw error
        error('Impossible to determine proper linear range for the given polar (see figure). Check CONSTANTS section to tweak the function.')
    else
        warning('LinearRange:BoundNotFound', 'Not able to calculate at least one bound of the lift linear range.\nAssuming linear range missing bound to be equal to min/max aoa available.');
    end
end

% Output indexes of linear range
linRangeInd = [aoaMinInd, aoaMaxInd];

end
