function [] = output(par, res)
%%
disp(' ')
disp(' ')
disp(' ')
disp('============================================================')
disp('Performance')
disp('============================================================')
disp('------------------------------------------------------------')
disp('Takeoff')
disp('------------------------------------------------------------')

disp(['<strong>Takeoff distance = </strong>', ...
num2str(par.x_tot), ' [m] = ', ...
num2str(res.to.x_g), ' [m] (ground) + ', ...
num2str(res.to.x_a), ' [m] (obstacle clearance)']);

disp(['<strong>Power required = </strong>', ...
num2str(res.to.P), ' [W]  /  ', ...
num2str(res.to.P/1000), ' [kW]  /  ', ...
num2str(res.to.P/746), ' [imp. hp]']);

disp(['<strong>Power of the engine = </strong>', ...
num2str(res.to.P_eng), ' [W]  /  ', ...
num2str(res.to.P_eng/1000), ' [kW]  /  ', ...
num2str(res.to.P_eng/746), ' [imp. hp]']);

disp(['<strong>Mean takeoff thrust = </strong>', ...
num2str(res.to.T_mean), ' [N]  /  ', ...
num2str(res.to.T_mean/1000), ' [kN]']);

disp(['<strong>TAS at takeoff = 1.1 * V_stall = </strong>', ...
num2str(res.to.V_LOF), ' [m/s]  /  ', ...
num2str(res.to.V_LOF * 3.6), ' [km/h]']);

disp(['<strong>Maximum lift coefficient at takeoff = </strong>', ...
num2str(res.to.C_L_max * par.takeoff_blown_wing_effect)]);

%
disp(['<strong>Blown wing effect: </strong>'])
disp(['<strong>  • Induced velocity = </strong>', ...
num2str(res.v_i_to), ' [m/s]']);

disp(['<strong>  • HLP power = </strong>', ...
num2str(res.P_BWE_to / 1000), ' [kW]  /  ', ...
num2str(res.P_BWE_to / 746), ' [imp. hp]']);

%
disp(['<strong>Total takeoff power = </strong>', ...
num2str(res.en.tot_to_power), ' [W]  /  ', ...
num2str(res.en.tot_to_power/1000), ' [kW]  /  ', ...
num2str(res.en.tot_to_power/746), ' [imp. hp]']);


disp(' ')
disp('------------------------------------------------------------')
disp('Climb')
disp('------------------------------------------------------------')

disp(['<strong>Power of the engine = </strong>', ...
num2str(res.climb.P_engine / 1000), ' [kW]  /  ', ...
num2str(res.climb.P_engine / 746), ' [imp. hp]']);

disp(['<strong>TAS during climb = </strong>', ...
num2str(res.climb.V_TAS), ' [m/s]']);

disp(['<strong>Climb rate = </strong>', ...
num2str(par.climb_rate), ' [m/s]']);

disp(['<strong>Climb angle = </strong>', ...
num2str(res.climb.angle), ' [°]']);

disp(['<strong>Climb C_L = </strong>', ...
num2str(par.C_L_climb), ' [°]']);

disp(['<strong>Time to climb = </strong>', ...
num2str(res.climb.time/60), ' [min]']);

% %
% disp(['<strong>Blown wing effect: </strong>'])
% disp(['<strong>  • Induced velocity = </strong>', ...
% num2str(res.v_i_climb), ' [m/s]']);
% 
% disp(['<strong>  • HLP power = </strong>', ...
% num2str(res.P_BWE_climb / 1000), ' [kW]  /  ', ...
% num2str(res.P_BWE_climb / 746), ' [imp. hp]']);


disp(' ')
disp('------------------------------------------------------------')
disp('Cruise')
disp('------------------------------------------------------------')

disp(['<strong>Speed = </strong>', ...
num2str(par.speed), ' [knots]  /  ', ...
num2str(res.cruise.V), ' [m/s]']);

disp(['<strong>Range = </strong>', ...
num2str(par.range), ' [nmi]  /  ', ...
num2str(res.cruise.range), ' [m]  /  ', ...
num2str(res.cruise.range/1000), ' [km]']);

disp(['<strong>Time = </strong>', ...
num2str(res.cruise.time / 3600), ' [h]']);

disp(['<strong>Thrust = </strong>', ...
num2str(res.cruise.T), ' [N]  /  ', ...
num2str(res.cruise.T / 1000), ' [kN] ']);

disp(['<strong>Thrust to weight ratio = </strong>', ...
num2str(res.cruise.T_to_W)]);

disp(['<strong>Power = </strong>', ...
num2str(res.cruise.P), ' [W]  /  ', ...
num2str(res.cruise.P / 1000), ' [kW]  /  ', ...
num2str(res.cruise.P / 746), ' [imp. hp]']);

disp(['<strong>Power of the engine = </strong>', ...
num2str(res.cruise.P_engine), ' [W]  /  ', ...
num2str(res.cruise.P_engine / 1000), ' [kW]  /  ', ...
num2str(res.cruise.P_engine / 746), ' [imp. hp]']);

disp(['<strong>Power of the generator = </strong>', ...
num2str(res.en.ICE_actual_power), ' [W]  /  ', ...
num2str(res.en.ICE_actual_power / 1000), ' [kW]  /  ', ...
num2str(res.en.ICE_actual_power / 746), ' [imp. hp]']);


disp(' ')
disp('------------------------------------------------------------')
disp('Landing')
disp('------------------------------------------------------------')

disp(['<strong>landing distance = </strong>', ...
num2str(res.landing.ground_distance), ' [m]']);

%%
disp(' ')
disp(' ')
disp(' ')
disp('============================================================')
disp('Batteries')
disp('============================================================')
disp('------------------------------------------------------------')
disp('Energy and power')
disp('------------------------------------------------------------')

%
% disp(['<strong>Electric engine nominal power = </strong>', ...
% num2str(res.en.P_nom_electric_engines / 1000), ' [kW]  /  ', ...
% num2str(res.en.P_nom_electric_engines / 746), ' [imp. hp]']);
% 
% disp(['<strong>ICE nominal power = </strong>', ...
% num2str(res.en.P_nom_ICE_engine / 1000), ' [kW]  /  ', ...
% num2str(res.en.P_nom_ICE_engine / 746), ' [imp. hp]']);

disp(['<strong>Total energy required = </strong>', ...
num2str(res.en.tot_energy / 1e6), ' [MJ]  /  ', ...
num2str(res.en.tot_energy / 3.6e6), ' [kWh]']);

%
disp(['<strong>Fuel total energy: </strong>'])
disp(['<strong>  • primary energy = </strong>', ...
num2str(res.en.fuel_tot_energy_primary / 1e6), ' [MJ]  /  ', ...
num2str(res.en.fuel_tot_energy_primary / 3.6e6), ' [kWh]']);

disp(['<strong>  • final energy = </strong>', ...
num2str(res.en.fuel_tot_energy_final / 1e6), ' [MJ]  /  ', ...
num2str(res.en.fuel_tot_energy_final / 3.6e6), ' [kWh]']);

disp(['<strong>  • efficiency = </strong>', ...
num2str(res.en.fuel_tot_energy_final / res.en.fuel_tot_energy_primary * 100), ' [%]']);

%
disp(['<strong>Battery total energy: </strong>'])
disp(['<strong>  • primary energy = </strong>', ...
num2str(res.en.battery_tot_energy_primary / 1e6), ' [MJ]  /  ', ...
num2str(res.en.battery_tot_energy_primary / 3.6e6), ' [kWh]']);

disp(['<strong>  • final energy = </strong>', ...
num2str(res.en.battery_tot_energy_final / 1e6), ' [MJ]  /  ', ...
num2str(res.en.battery_tot_energy_final / 3.6e6), ' [kWh]']);

disp(['<strong>  • efficiency = </strong>', ...
num2str(res.en.battery_tot_energy_final / res.en.battery_tot_energy_primary * 100), ' [%]']);

%
disp(['<strong>Battery final energy to total energy ratio = </strong>', ...
num2str(res.en.battery_tot_energy_final / res.en.tot_energy * 100), ' [%]']);

%
disp(['<strong>Flight phases energy consumption: </strong>'])
disp(['<strong>  • Takeoff = </strong>', ...
num2str(res.en.tot_to_energy / res.en.tot_energy * 100), ' [%]']);

disp(['<strong>  • Climb = </strong>', ...
num2str(res.en.tot_climb_energy / res.en.tot_energy * 100), ' [%]']);

disp(['<strong>  • Cruise and approach = </strong>', ...
num2str(res.en.tot_cruise_energy / res.en.tot_energy * 100), ' [%]']);


disp(' ')
disp('------------------------------------------------------------')
disp('Fuel and battery mass')
disp('------------------------------------------------------------')
%
disp(['<strong>Battery: </strong>'])
disp(['<strong>  • mass (initial input) = </strong>', ...
num2str(par.mass_battery), ' [kg]']);

disp(['<strong>  • mass (converged to) = </strong>', ...
num2str(res.en.mass_battery), ' [kg]']);

disp(['<strong>  • volume = </strong>', ...
num2str(res.en.volume_battery * 1000), ' [l]']);

disp(['<strong>  • CO2 emissions = </strong>', ...
num2str(res.en.CO2_battery), ' [kg CO2]']);

%
disp(['<strong>Fuel: </strong>'])
disp(['<strong>  • Fuel required mass = </strong>', ...
num2str(res.en.fuel_req_mass), ' [kg]']);

disp(['<strong>  • volume = </strong>', ...
num2str(res.en.volume_fuel * 1000), ' [l]']);

disp(['<strong>  • CO2 emissions = </strong>', ...
num2str(res.en.CO2_kerosene), ' [kg CO2]']);

%
disp(['<strong>Total CO2 emissions = </strong>', ...
num2str(res.en.CO2_battery + res.en.CO2_kerosene), ' [kg CO2]']);



%%
disp(' ')
disp(' ')
disp(' ')
disp(' ')




end