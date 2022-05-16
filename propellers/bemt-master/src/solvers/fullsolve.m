function [Elem] = fullsolve(Elem, Blades, Param)
% FULLSOLVE Bemt solver following Stahlhut equations (solved with fsolve).
%
% Calculates the inflow velocity ratios without assumptions on the angles,
% the velocities or the forces.
% -----
%
% Synopsis: [Elem] = FULLSOLVE(Elem, Blades, Param)
%
% Input:    Elem     = (required) Variables computed for any individual blade element
%           Blades   = (required) Blade parameters defined in config file
%           Param    = (required) General parameters updated with flight conditions
%
% Output:   Elem     = Updated Elem struct with inflow angles and velocities
%
% Calls:    prandtlLoss   = Calculate tip losses with Prandtl formula
%           getCoeffValue = Get cl and cd from airfoil data
%
% See also: BEMT, PRANDTLLOSS, GETCOEFFVALUE.
%

% -----
% Stahlut, C. W., & Leishman, J. G. (2012). "Aerodynamic design optimization of
% proprotors for convertible-rotor concepts". American Helicopter Society
% 68th Annual Forum. Fort Worth, Texas, USA.
%
% -----
% This function is one set of equations that can be used with the BEMT.
% Author: Thomas Lambert <t.lambert@uliege.be>
% ULiege - Aeroelasticity and Experimental Aerodynamics
% MIT License
% https://gitlab.uliege.be/thlamb/bemt
% -------------------------------------------------------------------------

% Initialize Elem struct in the correct order
Elem.Lambda =0;
Elem.Xi =0;
Elem.InducedVeloc =0;
Elem.TotalVeloc =0;
Elem.InflowAngle =0;

% Abbreviations
tipSpeed = Blades.RADIUS*Blades.OMEGA;

% Calculate for each section one at a time (due to fsolve limitations)
for i = 1:Blades.nELEM
    % Solve equation system
    try
        % Lower and upper bounds for the induced angle
        unk0(1) = 1e-12;    % Avoid singularity (g(0)=0);
        unk0(2)= (Elem.Pitch(i)-Blades.Polar.alpha(2)); % Bound the inflow in order to ensure that AOA is on the polars
        
        [phi,~] = fzero(@(unk) gfun(unk, Blades.nBLADES, Elem.Solidity(i), Blades.OMEGA, Param.Air.AXIAL_VELOC, Blades.Polar, Elem.Pitch(i), Elem.Re(i), Elem.Ypos(i), Elem.r(i), Elem.r0, Param.Simul.LOSSES), unk0);
        
        % Output values
        Elem.InflowAngle(i) = phi;
    catch ME
        % If issue with the solver, calculate the value of g(phi) for this section
        phi_test = -pi/2:0.01:pi/2;
        gval_test = zeros(size(phi_test));
        for j = 1:length(phi_test)
            gval_test(j) = gfun(phi_test(j), Blades.nBLADES, Elem.Solidity(i), Blades.OMEGA, Param.Air.AXIAL_VELOC, Blades.Polar, Elem.Pitch(i), Elem.Re(i), Elem.Ypos(i), Elem.r(i), Elem.r0, Param.Simul.LOSSES);
        end
        
        fprintf(2,'************************************************* \n');
        fprintf(2,'Can not find solution for g(phi) at element %d \n', i');
        fprintf(2,'************************************************* \n');
        
        % Plot g(phi) for the section that caused the issue
        figname = ['g(phi) for Element',num2str(i)];
        figname_latex = ['$g(\phi)$ for Element',num2str(i)];
        figure('Name',figname)
        plot(rad2deg(phi_test),gval_test)
        grid on
        hXLabel=xlabel('Inflow angle, $\phi$ [deg]');
        hYLabel=ylabel('$g(\phi)$');
        hTitle=title(figname_latex);
        set(gca, 'FontName','Helvetica');
        set(gca, 'Fontsize', 14);
        set([hXLabel, hYLabel, hTitle], 'interpreter','latex', 'Fontsize', 14);
        set(gcf, 'PaperPositionMode', 'auto');
        
        % Rethrow excepetion to stop the process
        rethrow(ME)
    end
end

aoa = Elem.Pitch-Elem.InflowAngle;

% Get new estimates for the coefficients
cl = getCoeffValue('cl', Blades.Polar, aoa, Elem.Re);
cd = getCoeffValue('cd', Blades.Polar, aoa, Elem.Re);

gamma = atan(cd/cl);

% Calculate inflow ratios and AOA
b2phi = cos(Elem.InflowAngle) + 1./(8*Elem.r) .* Elem.Solidity .*cl .*sec(gamma) .*csc(abs(Elem.InflowAngle)).*sin(Elem.InflowAngle+gamma);
Elem.Xi = Elem.r .* cos(Elem.InflowAngle)./b2phi;
Elem.Lambda = Elem.Xi .* tan(Elem.InflowAngle);
Elem.AOA = Elem.Pitch-Elem.InflowAngle;

Elem.InducedVeloc = Elem.Lambda*tipSpeed - Param.Air.AXIAL_VELOC;

Elem.TotalVeloc = sqrt((Elem.Lambda*tipSpeed).^2+(Elem.Xi*tipSpeed).^2);

end


function gphi = gfun(phi, nBlades, sigma, omega, v_inf, Polar, pitch, reyn, y_i, r_i, r0, inclLoss)
% GFUN Equation for g(phi) where phi = Inflow angle

% Loss factor (separate swirl and axial components)
if inclLoss == 1 || inclLoss == 2
    lossFact = prandtlLoss(nBlades, r_i, r0, phi, inclLoss);
    K_T = 1-(1-lossFact).*cos(phi);
    K_P = 1-(1-lossFact).*sin(phi);
else
    K_T = 1;
    K_P = 1;
end

% Aerodynamic coefficients
cl = getCoeffValue('cl', Polar, pitch-phi, reyn);
cd = getCoeffValue('cd', Polar, pitch-phi, reyn);

% Equation to solve (eq 32 in Stahlhut paper)
if cl == 0
    gphi = omega*y_i*sin(phi) - v_inf .*cos(phi) ;
else
    gamma = atan(cd/cl);
    gphi = (omega*y_i.*sin(phi) - v_inf .*cos(phi)) .*sin(phi) ...
           - sign(phi) * 1/(8*r_i) * sigma * cl *sec(gamma) ...
           *(omega * y_i ./K_T .* cos(phi+gamma) + v_inf ./K_P .* sin(phi+gamma)) ;
end

end
