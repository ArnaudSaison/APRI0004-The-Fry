function [par] = parameters()
%**************************************************************************
%**************************************************************************
% Parameters
%**************************************************************************
%**************************************************************************
% This file contains all the parameters used during the conceptual design
% phase. Some of them are initialization values for iterative processes.
% Only the code outputs should be used as a base for other design phases.
%**************************************************************************
%%
%==========================================================================
% Performance
%==========================================================================
%--------------------------------------------------------------------------
% Runway
%--------------------------------------------------------------------------
par.runway = 300;                       % [ft]
par.runway = 300 / 3.281;               % [m] conversion to meters
par.x_tot = par.runway * 0.95;          % safety margin (around 85 m)
par.h = 50 * 0.3048;                    % obtacle height (feet to meters conversion)
par.to_alt = 0 / 3.281;

%--------------------------------------------------------------------------
% Physical constants
%--------------------------------------------------------------------------
par.g = 9.81;                           % gravity
par.rho_0 = 1.225;                      % air density at ground
[~, ~, ~, par.rho_0] = atmosisa(par.to_alt);

%--------------------------------------------------------------------------
% Aircraft characteristics and geometry
%--------------------------------------------------------------------------
par.M_to = 1600;                        % takeoff mass
par.W_to = par.M_to * par.g;            % takeoff weight
par.S = 15;                             % [m^2] wing surface

par.AR = 8;                             % between 7 and 9
par.mu = 0.02;                          % ground friction coeff at takeoff on concrete
par.sweep_angle = 0;                    % for speed range

par.C_L_max = 2.8;                      % before sweep
par.takeoff_blown_wing_effect = 6/par.C_L_max;  % coefficient from blown lift effet
par.C_L_blown = par.C_L_max * par.takeoff_blown_wing_effect; % before sweep

%--------------------------------------------------------------------------
% Assumptions from estimated shape and size of aircraft
%--------------------------------------------------------------------------
par.C_D0_to = 0.02;                     % drag coeff (assumption)
par.e = 0.75;
par.k = 1 / (par.e * pi * par.AR);

%--------------------------------------------------------------------------
% Engines (must be refined using the BEMT code)
%--------------------------------------------------------------------------
par.prop_efficiency_to = 0.53;           % at take-off
par.prop_efficiency_cruise = 0.97;       % 

%--------------------------------------------------------------------------
% Propeller (must be refined using the BEMT code)
%--------------------------------------------------------------------------
par.prop_diameter = 1;                  % [m]
par.prop_surface = pi * par.prop_diameter^2; % [m^2]

par.surface_blown = 0.75;               % [-] percentage of the wing that is actually blown by the propellers

%--------------------------------------------------------------------------
% Cruise
%--------------------------------------------------------------------------
par.cruise_altitude = 13000 * 0.3048;       % [m]
[~, ~, ~, par.rho_cruise] = atmosisa(par.cruise_altitude);
par.C_L_minD = 0.33;
%par.C_D0_min = 0.029;
par.C_D0_min = 0.0318;
par.C_D0_min = 0.0318;
%par.C_D0_min = 0.00619;
par.range = 300 * 1.852e3;              % [m] total range of the aircraft
par.speed = 170;                        % [knots] max speed of aircraft

%--------------------------------------------------------------------------
% Climb constraints
%--------------------------------------------------------------------------
%par.cruise_altitude = 6800;
par.climb_rate = 1500;                      % [ft/min] imposed climb rate = V_z
par.climb_rate = par.climb_rate / 3.28084 / 60; % [m/s] conversion

[~, ~, ~, par.rho_climb_1] = atmosisa(0);
[~, ~, ~, par.rho_climb_2] = atmosisa(par.cruise_altitude);

par.climb_velocity_fact = 1;                % for test puposes only

% drag and thus power can be optimized by using the appropriate lift
% coefficient
[~, ~, ~, par.rho_climb] = atmosisa(0);
par.climb_blown_wing_effect = 1;          % []
%par.C_L_climb = par.C_L_max;
par.C_L_climb = 1.439;                        % [] based on climb_optimizer.m
par.C_D0_climb = par.C_D0_min * 1.2;        % []

%--------------------------------------------------------------------------
% Landing 
%--------------------------------------------------------------------------
par.mean_deceleration = 0.7*par.g; % tabulated value for light aircraft with breaks
%par.approach_angle = 12.28*pi/180; % same as climb angle


%%
%==========================================================================
% Energy
%==========================================================================
%--------------------------------------------------------------------------
% Battery and fuel parameters
%--------------------------------------------------------------------------
par.battery_max_SOC = 0.80;                 % [-] max state of charge of the battery (min 15%, max 95% = 80% useful capacity)

% energy density
par.en_dens_Li_ion = 700 * 3600;            % [J/kg] assumption based on doubling in energy density before 2031
par.en_dens_kerosene = 43.2e6;              % [J/kg]

% mass density
par.mass_dens_li_ion = 500;                 % [kg/m^3]
par.mass_dens_kerosene = 821;               % [kg/m^3]

% CO2 emissions
par.kerosene_CO2 = 3 / par.en_dens_kerosene;% kg of CO2 / J of fuel energy
par.electricity_CO2 = 0.012 / 3.3e6;        % kg of CO2 / J of electric energy


%--------------------------------------------------------------------------
% energy storage mass (used as initial value for iteration)
%--------------------------------------------------------------------------
par.mass_battery = 250;                     % [kg]
par.mass_battery_margin = 1;                % [kg] stop iterating when less than x kg left available
par.mass_battery_max_iter = 50;             % stop iterating if doesn't converge

%--------------------------------------------------------------------------
% structure mass
%--------------------------------------------------------------------------
par.mass_dry = 600;                         % [kg]
par.mass_payload = 400;                     % [kg] (190 lbs passenger or pilot + 30 lbs luggage) * 4

%--------------------------------------------------------------------------
% Power assumptions
%--------------------------------------------------------------------------
%par.mass_power_ratio_elec = 30e3;           % [W / kg] based on 2030 engines
%par.mass_power_ratio_ICE = 2.2e3;           % [W / kg]

par.nb_HLP_engines = 12;                     % number of high lift propeller engines
par.mass_power_elec = par.nb_HLP_engines * 10 + 2 * 50;% [kg] 8 * 10
par.mass_power_ICE = 150;                    % [kg] 2 * 50

%(157 * 745.6) / ((0.3*157)*1e-3 * par.mass_dens_kerosene * par.en_dens_kerosene / 3600)
par.peak_power_elec = 1.3;                   % [-] coefficient multiplying electric engine to have peak power
par.generator_fuel_efficiency = 0.40;        % [-] from primary energy to electric energy at cruise altitude
%par.generator_fuel_efficiency = (200 * 745.6) / (33e-3 * par.mass_dens_kerosene * par.en_dens_kerosene / 3600)
par.charging_efficiency = 0.98;              % [-] charging efficiency between ICE and battery
par.electric_engines_efficiency = 0.96;      % [-]
par.ICE_margin = 1.1;                        % [-] slected engine at least x times more powerfull than theoretically needed

%--------------------------------------------------------------------------
% problem maximum mass
%--------------------------------------------------------------------------
par.total_max_mass = par.M_to;


end