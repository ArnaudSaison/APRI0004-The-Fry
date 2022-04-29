function [res] = performance(par)
%==========================================================================
% FLIGHT PERFORMANCE
%==========================================================================
%% ------------------------------------------------------------------------
% TAKEOFF
%--------------------------------------------------------------------------
% System of symbolic equations
syms x_a x_g T_mean mu_p P V_LOF C_L_max P_eng

runway_length = par.x_tot == x_a + x_g; % total length of runway
air_distance_cond = x_a > 0;
air_distance = x_a == V_LOF^2 / (par.g*sqrt(2)) + par.h / (0.9 * T_mean / par.W_to - 0.3 / sqrt(par.AR));

ground_distance_cond = x_g > 0;
ground_distance = x_g == V_LOF^2 / (2 * par.g) / (T_mean / par.W_to - mu_p);

friction_coeff = mu_p == par.mu + 0.72 * par.C_D0_to / C_L_max;

power = P == T_mean * V_LOF / sqrt(2);
power_eng = P_eng == P / par.prop_efficiency_to;

liftoff_speed = V_LOF == 1.1 * sqrt(2 * par.W_to / (par.rho_0 * par.S * C_L_max * par.takeoff_blown_wing_effect));
lift_coeff = C_L_max == par.C_L_max * cos(par.sweep_angle);

% Solving the system of equations
sol = solve([runway_length, ...
             air_distance, ...
             air_distance_cond, ...
             ground_distance, ...
             ground_distance_cond, ...
             friction_coeff, ...
             power, ...
             liftoff_speed, ...
             lift_coeff, ...
             power_eng], ...
             [x_a, x_g, T_mean, mu_p, P, V_LOF, C_L_max, P_eng]);

if sol.x_a <= 0
    warning('Plane will not take off. Please change parameters')
end

% results
res.to.x_a = double(sol.x_a);
res.to.x_g = double(sol.x_g);
res.to.T_mean = double(sol.T_mean);
res.to.P = double(sol.P);
res.to.V_LOF = double(sol.V_LOF);
res.to.C_L_max = double(sol.C_L_max);
res.to.C_L_blown = res.to.C_L_max * par.takeoff_blown_wing_effect;
res.to.P_eng = double(sol.P_eng);


%% ------------------------------------------------------------------------
% CLIMB
%--------------------------------------------------------------------------
% P_a = P_r + W * V_z
% where:
% • P_r = cruise power at same airspeed for level flight
% • V_z = rate of climb = V * sin gamma
% • P_a = availeble power for climb
% -> we are looking for V_TAS_climb

res.climb.time = par.cruise_altitude / par.climb_rate; % [s] time to reach cruise altitude

% the loop solves the system of equations (usually takes 4 to 5 iterations)

res.climb.gamma_climb = 1;
climb_las_iter = 0;
climb_iter_stop_criterion = 1e-4;
climb_iter = 1;
climb_max_iter = 10;

while ((cos(res.climb.gamma_climb) - climb_las_iter) > climb_iter_stop_criterion) && climb_iter < climb_max_iter
    climb_las_iter = cos(res.climb.gamma_climb);
    
    res.climb.V_TAS = 1 * par.climb_velocity_fact * sqrt(2 * par.W_to * cos(res.climb.gamma_climb) / (par.rho_climb * par.S * par.C_L_climb * par.climb_blown_wing_effect));
    res.climb.C_D = par.C_D0_min + (par.C_L_climb)^2 * par.k;
    
    res.climb.D_1 = 1/2 * par.rho_climb_1 * res.climb.V_TAS^2 * res.climb.C_D * par.S;
    res.climb.D_2 = 1/2 * par.rho_climb_2 * res.climb.V_TAS^2 * res.climb.C_D * par.S;
  
    
    res.climb.T_1 = res.climb.D_1 + par.W_to * sin(res.climb.gamma_climb);
    res.climb.T_2 = res.climb.D_2 + par.W_to * sin(res.climb.gamma_climb);
    
    % we can compute gamma
    res.climb.gamma_climb = asin(par.climb_rate / res.climb.V_TAS); % [rad]
    
    
    climb_iter = climb_iter + 1;
end

if climb_iter >= climb_max_iter
    warning('Impossible to find a climb speed.');
end

% required power for climb
res.climb.P_1 = res.climb.T_1*res.climb.V_TAS;
res.climb.P_2 = res.climb.T_2*res.climb.V_TAS;

res.climb.P_r = 1/2 * par.rho_climb * res.climb.V_TAS^3 * res.climb.C_D * par.S;
res.climb.P_available = res.climb.P_r + par.W_to * res.climb.V_TAS * sin(res.climb.gamma_climb);
res.climb.angle = double(res.climb.gamma_climb / pi * 180); % [°] climb angle

res.climb.P_engine = res.climb.P_available / par.prop_efficiency_cruise;

% "ground" distance covered during climb
res.climb.distance = res.climb.time * cos(res.climb.angle) * res.climb.V_TAS;


%% ------------------------------------------------------------------------
% CRUISE
%--------------------------------------------------------------------------
res.cruise.V = par.speed * 1.852 / 3.6;   % [m/s] cruise speed
res.cruise.range = par.range - res.climb.distance;   % [m] range of the cruise
res.cruise.time = res.cruise.range / res.cruise.V; % [s]

%Cruise power requirements
%res.cruise.C_D = par.C_D0_min + par.C_L_minD^2 * par.k; % total drag
res.cruise.C_D = 0.0318;

res.cruise.T = 1/2 * par.rho_cruise * res.cruise.V^2 * res.cruise.C_D * par.S; % because T = D at equilibrium
res.cruise.T_to_W = res.cruise.T / par.W_to; % thrust to weight ratio
res.cruise.P = res.cruise.T * res.cruise.V; % usable power delivered by the engine after propeller efficiency
res.cruise.P_engine = res.cruise.P / par.prop_efficiency_cruise;
res.cruise.P_generator = res.cruise.P / par.prop_efficiency_cruise / par.charging_efficiency / par.electric_engines_efficiency;


%% ------------------------------------------------------------------------
% LANDING
%--------------------------------------------------------------------------
% considering a touchdown velocity of 1.1 V_stall = V_LOF
% res.landing.ground_distance = res.to.V_LOF^2 / 2 * par.mean_deceleration;
% 
% if res.landing.ground_distance > par.x_tot
%     warning("Runway not long enough for landing")
% end
res.landing.V_stall = sqrt((2*(par.W_to))/(par.rho_0*par.S*par.C_L_blown));
res.landing.V_a = 1.3*res.landing.V_stall;
res.landing.V_TD = 1.1*res.landing.V_stall;
res.landing.dist_g = res.landing.V_TD^2/(2*par.mean_deceleration);

res.landing.C_D = par.C_D0_min + (par.C_L_blown)^2 * par.k;
res.landing.D = 1/2 * par.rho_0 * res.landing.V_a^2 * res.landing.C_D * par.S;
res.landing.L = 1/2 * par.rho_0 * res.landing.V_a^2 * par.C_L_blown * par.S;
res.landing.approach_angle = atan2(res.landing.D,par.W_to);

R = 0.1512*res.landing.V_stall^2;
h_f = R*(1-cos(res.landing.approach_angle));                           % flare height 
res.landing.dist_a = (par.h - h_f)/tan(res.landing.approach_angle);      % approach distance

res.landing.dist_f = R*sin(res.landing.approach_angle);      %flare distance

% res.landing.dist_fr = 0.9*res.landing.V_TD;         % free-roll distance
res.landing.dist_fr = 0.9*res.landing.V_TD;         % free-roll distance

res.landing.dist_tot = res.landing.dist_g + res.landing.dist_a + ...
res.landing.dist_f + res.landing.dist_fr; % total landing distance
res.landing.approach_angle = res.landing.approach_angle * 180/pi;
if res.landing.dist_tot > par.x_tot
warning("Runway not long enough for landing")
end


%% ------------------------------------------------------------------------
% Blown wing effect
%--------------------------------------------------------------------------
res.v_i_to = sqrt(par.W_to * 2 / (par.rho_0 * res.to.C_L_max * par.S * par.surface_blown)) - res.to.V_LOF;

% res.P_BWE_to = 2 * par.rho_0 * par.prop_surface * (res.to.V_LOF + res.v_i_to) * res.v_i_to^2;
% res.T_BWE_to = 2 * par.rho_0 * par.prop_surface * (res.to.V_LOF + res.v_i_to) * res.v_i_to;
% 
% res.P_tot_BWE_to = res.P_BWE_to * par.nb_HLP_engines;
% res.T_tot_BWE_to = res.T_BWE_to * par.nb_HLP_engines;

res.P_BWE_to = 8.15e3;
res.T_BWE_to = 247;

res.P_tot_BWE_to = res.P_BWE_to * par.nb_HLP_engines;
res.T_tot_BWE_to = res.T_BWE_to * par.nb_HLP_engines;

% res.v_i_climb = sqrt(par.W_to * 2 / (par.rho_0 * par.C_L_climb * par.S * par.surface_blown)) - res.climb.V_TAS;
% res.P_BWE_climb = 2 * par.rho_0 * par.prop_surface * (res.climb.V_TAS + res.v_i_climb) * res.v_i_climb^2;


end