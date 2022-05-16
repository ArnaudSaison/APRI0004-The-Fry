function [res] = batteries(res, par)
%==========================================================================
% Energy required for flight
%==========================================================================
%% ------------------------------------------------------------------------
% Time of each flight part
%--------------------------------------------------------------------------
% estimated time during which power is used
res.en.to_time = 120; %[s]
res.en.climb_time = res.climb.time; % [s] higher estimate to take into account taxi, takeoff and climb


%% ------------------------------------------------------------------------
% Energy
%--------------------------------------------------------------------------
res.en.tot_to_power = res.P_tot_BWE_to + res.to.P_eng;

% energy consumption at end of motor
res.en.tot_to_energy = res.en.to_time * res.en.tot_to_power;
res.en.tot_climb_energy = res.en.climb_time * res.climb.P_available; % climb and takeoff (which is about the )
res.en.tot_cruise_energy = res.cruise.time * res.cruise.P_engine * 1.1; % 10% margin required by AIAA

res.en.tot_energy = res.en.tot_to_energy + res.en.tot_climb_energy + res.en.tot_cruise_energy;


%% ========================================================================
% Repartition of mass for fuel and batteries
%==========================================================================
%% ------------------------------------------------------------------------
% Mass for power requirements (engines)
%--------------------------------------------------------------------------
% mass of electric engines based on cruise requirements
% this is assuming the electric engines can momentarily provide takeoff
% power
%res.en.mass_elec_engines = res.to.P_eng / par.mass_power_ratio_elec / par.peak_power_elec;
%res.en.mass_fuel_engine = res.cruise.P_engine / par.mass_power_ratio_ICE;

res.en.mass_elec_engines = par.mass_power_elec;
res.en.mass_fuel_engine = par.mass_power_ICE;

res.en.mass_engines_tot = res.en.mass_elec_engines + res.en.mass_fuel_engine;

% Actual power of engines
%res.en.P_nom_electric_engines = res.to.P_eng / par.peak_power_elec; % nominal power of electric engines
%res.en.P_nom_ICE_engine = res.cruise.P_engine;

% mass available for fuel after power requirements
res.en.fuel_available_mass = par.total_max_mass - res.en.mass_engines_tot - par.mass_dry - par.mass_payload - res.en.mass_battery; % mass available for fuel


%% ------------------------------------------------------------------------
% Mass for energy requirements (fuel)
%--------------------------------------------------------------------------
% battery taking into account the efficiency of the engines -> total useful
% energy the battery can deliver = final energy
res.en.battery_tot_energy_final = res.en.mass_battery * par.en_dens_Li_ion * par.battery_max_SOC * par.electric_engines_efficiency;

% fuel
res.en.fuel_tot_energy_final = res.en.tot_energy - res.en.battery_tot_energy_final; % energy needed that must be fuel
res.en.fuel_tot_energy_primary = res.en.fuel_tot_energy_final / par.generator_fuel_efficiency / par.charging_efficiency / par.electric_engines_efficiency;

% mass needed for fuel considering energy requirements
res.en.fuel_req_mass = res.en.fuel_tot_energy_primary / par.en_dens_kerosene;

% throw warning in case the battery already contains all energy
if res.en.fuel_req_mass <= 0
    warning('Battery already meets all energy requirements, no fuel needed :)')
end

% mass left after fuel has been added
% (negative mass indicates the plane can't meet requirements and a new iteration must be run)
res.en.mass_left_after_fuel = res.en.fuel_available_mass - res.en.fuel_req_mass; % [kg]

% throw warning in case plane can't fly :(
if res.en.mass_left_after_fuel <= 0
    warning("The plane is too heavy, it won't fly :(")
end

% Ratio between battery energy and tot energy
res.en.tot_ernergy_battery_ratio = res.en.battery_tot_energy_final / res.en.tot_energy * 100; % [%]

% ernergy contained in the battery = primary energy
res.en.battery_tot_energy_primary = res.en.mass_battery * par.en_dens_Li_ion; % [J]


%% ------------------------------------------------------------------------
% Volume
%--------------------------------------------------------------------------
res.en.volume_battery = res.en.mass_battery / par.mass_dens_li_ion;
res.en.volume_fuel = res.en.fuel_req_mass / par.mass_dens_kerosene;


%% ------------------------------------------------------------------------
% Carbon emission
%--------------------------------------------------------------------------
res.en.CO2_battery = res.en.battery_tot_energy_primary * par.electricity_CO2 * par.battery_max_SOC;
res.en.CO2_kerosene = res.en.fuel_tot_energy_primary * par.kerosene_CO2;


%% ------------------------------------------------------------------------
% ICE engine actual power
%--------------------------------------------------------------------------
% total final energy the ICE must provide is res.en.fuel_tot_energy_final
% amount of time the ICE must run is the cruise time
res.en.ICE_time = res.cruise.time;
res.en.ICE_actual_power = res.en.fuel_tot_energy_final / res.en.ICE_time;



end