function radiusChecker(Blades, maxMach)
%RADIUSCHECKER checks if the maximum mach number is not exceeded
%   Detailed explanation goes here
% Blade radius
Blades.maxMach = maxMach;       % Maximum Mach number (0.6 for wooden, 0.75-0.8 for metal and composite)
Blades.RADIUS_MAX = Blades.maxMach * 343 / (Blades.OMEGArpm /60*2*pi);    % Radius (from HUB center to tip), [m]
Blades.Mach = Blades.RADIUS * (Blades.OMEGArpm /60*2*pi) / 343;

disp(['Maximum blade radius (for max Mach of ', num2str(Blades.maxMach) ,') = ', num2str(Blades.RADIUS_MAX, '%.2f'), ' m']);
disp(['Set radius = ', num2str(Blades.RADIUS,'%.2f')]);
disp(['Resulting Mach = ', num2str(Blades.Mach,'%.2f')]);
end

