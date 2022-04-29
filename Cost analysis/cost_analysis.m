%% Cost analysis 
close all; clear; clc;

w_tot = 1600; %[kg]
w_luggages = 40.8; %[kg]
w_people = 344.73; %[kg]
w_fuel = 53; %[kg]
w_batteries = 175; %[kg]

%% Development cost (Gudmundsson)

% Given 

W_empty = w_tot - w_luggages - w_people - w_fuel - w_batteries; %[kg]
W_airframe = 600 * 2.20462262; %0.65 * W_empty * 2.20462262; %[lbs]

V_H = 170; %[kts]
F_CF = 1; %For simple flaps system = 1
F_CERT = 1; % =1 for the case of the FAR23
f_comp_standard = 0.9; % Amount of the airframe made of composite
F_COMP_standard = f_comp_standard + 1; 
F_COMP_MFG = 0.25*f_comp_standard + 1; 
F_COMP_DEV = 0.5*f_comp_standard + 1; 
F_PRESS = 1; %Unpressurized aircraft
F_TAPER = 1; %For a tapered wing = 1
CPI_2012 = 1.24; %Consumer price index relative to the year 2012
N_p = 5; %Number of prototypes
N_pp = 1; %Number of engines (generator)
N_ee = 15; %Number of electric engines
P_BHP = 200; %Rated brake-horsepower 
eta_p = 0.95; %Efficiency of electic motors
P_b = 282; %Power of the big motors [hp]
P_s = 40; %Power of the small motors [hp]
P_SHP_b = P_b*eta_p; %Rated shaft-horsepower big
P_SHP_s = P_s*eta_p; %Small
D_p_b = 1.8 * 3.2808399; %Diameter of the big propellers [ft]
D_p_s = 0.7 * 3.2808399; %Diameter of the small propellers [ft]
%inflation = 1.24; %https://www.bls.gov/data/inflation_calculator.htm
Q_m = 4:1:10; %Estimated production rate in number of aircraft per month

for i = 1:length(Q_m)
    
    N = Q_m(i) * 12 * 5; %Number of planned aircraft to be produced over 5 years period
    
    %Hours
    %Engineering man-hours
    dev.H_eng = 0.0396 * W_airframe.^(0.791) .* V_H.^(1.526) .* N^(0.183) * F_CERT * F_CF * F_COMP_standard * F_PRESS;
    
    %Tooling man-hours
    dev.H_tool = 1.0032 * W_airframe.^(0.764) .* V_H.^(0.899) .* N^(0.178) * Q_m(i)^(0.066) * F_TAPER * F_CF * F_COMP_standard * F_PRESS;
    
    %Manufacturing labor man-hours
    dev.H_mfg = 9.6613 *  W_airframe.^(0.74) .* V_H.^(0.543) .* N^(0.524) * F_CERT * F_CF * F_COMP_MFG;
    
    %Cost per hour
    cost_h_eng = 92; %[$/h]
    cost_h_tool = 61; %[$/h]
    cost_h_mfg = 53; %[$/h]
    
    %Costs
    dev.c_eng = 2.0969 * dev.H_eng .* cost_h_eng * CPI_2012 ; %[$]
    dev.c_tool = 2.0969 * dev.H_tool .* cost_h_tool * CPI_2012 ; %[$]
    dev.c_mfg = 2.0969 * dev.H_mfg .* cost_h_mfg * CPI_2012 ; %[$]
    
    %Other costs
    %Development support
    dev.dev_s = 0.06458 * W_airframe.^(0.873) .* V_H.^(1.89) .* N_p^(0.346) * CPI_2012 * F_CERT * F_CF * F_COMP_DEV * F_PRESS ;
    
    %Flight test operations
    dev.flight_test = 0.009646 * W_airframe.^(1.16) .* V_H.^(1.3718) .* N_p^(1.281) * CPI_2012 * F_CERT ;
    
    %Quality control
    dev.quality_control = 0.13 * dev.c_mfg * F_CERT * F_COMP_DEV;
    
    %Materials 
    dev.materials = 24.896 * W_airframe.^(0.689) .* V_H.^(0.624) .* N^(0.792) * CPI_2012 * F_CERT * F_CF * F_PRESS;
    
    %Certification
    dev.certification(i) = dev.c_eng + dev.dev_s + dev.flight_test + dev.c_tool;
    
    %Avionics
    dev.avionics = 15000;
    
    %Landing gears 
    dev.landing_gears = -7500;
    
    %Power plants
    %generator = 174 * N_pp .* P_BHP .* CPI_2012 * inflation;
    generator = 35000;
    batteries = 12250;
    electric_motors = 3 * 13000 + 12 * 3100;
    %electric_motors = (2 * 174 * P_SHP_b * CPI_2012 * inflation + 12 * 174 * P_SHP_s * CPI_2012 * inflation)/3;
    big_propellers = 3145 * 2 .* CPI_2012 ; %Constant speed propeller 
    small_propellers = 3145 * 4 .* CPI_2012 ; % Fixed pitch propeller

    dev.power = big_propellers + small_propellers + electric_motors + generator + batteries;
    
    dev.tot(i) = dev.certification(i) + dev.c_mfg + dev.materials + dev.quality_control + N*dev.avionics + N*dev.power + N*dev.landing_gears; %Total development cost
    
    %Liability 
    dev.liability = 0.13 * dev.tot(i);
    
    dev.tot(i) = dev.tot(i) + dev.liability;
    
    dev.cost_per_aircraft(i) = dev.tot(i)/N;
    
    dev.variable_cost(i) = dev.c_mfg + dev.quality_control + dev.materials + N*dev.landing_gears + N*dev.power + N*dev.avionics + dev.liability;
end

% figure
% plot(Q_m,dev.cost_per_aircraft * 10^(-6),'LineWidth',1)
% grid on
% xlabel('Number of units produced per month [-]','fontsize',16,'interpreter','Latex')
% ylabel('Unit production cost [millions \$]','fontsize',16,'interpreter','Latex')

%% Selling price

for i = 1:length(Q_m)
    
    N = Q_m(i) * 12 * 5;
    selling_price(i) = dev.cost_per_aircraft(i)*1.15;
    
    N_BE(i) = dev.certification(i)/(selling_price(i) - (dev.variable_cost(i)/N));
end

% figure
% plot(Q_m,selling_price * 10^(-6),'LineWidth',1)
% grid on
% xlabel('Number of units selled [-]','fontsize',16,'interpreter','Latex')
% ylabel('Unit selling price [millions \$]','fontsize',16,'interpreter','Latex')

%% Operational cost (8 flights/day, 1h10 per flight (mean))

BHP_cruise = 160; %[hp]
R_STOR = 250; %[$/month] Assumption of 250$/month of storage
R_crew = 70; %[$/h] 
R_AP = 53 * CPI_2012; % Hourly rate for a certified Airframe and Powerplant (A&P) mechanic (typ. $53or67 per hr)
Q_FLGT = 2800; %[h] Number of flight hours per year
route_day = 8; %[-] Number of route per day
oper_days = 300; %[-] Number of operationnal days
energy_cons_battery = 0.6 * 441000000; %[J] Energy consumption of the batteries per route
price_kWh = 0.15; %[$/kWh] (2020) https://www.statista.com/statistics/263492/electricity-prices-in-selected-countries/
Thrust = 1805 * 0.2248; %[lbf]
cons = 30 * 2.20462; %[lb/h] 
SFC_cruise = cons / Thrust; %Typical specific fuel consumption during cruise
R_fuel = 5; %[$/gallon] %Price of fuel (04/04/22, mean, AVGAS 100ll) http://100ll.com/ https://www.av-fuel.com/map/#

% Maintenance cost (per year)

F1 = 0; 
F2 = 0;
F3 = 0;
F4 = 0.02;
F5 = 0.04;
F6 = 0.01;
F7 = 0; 
F8 = 0;

F_MF = 0.3 + F1 + F2 + F3 + F4 + F5 + F6 + F7 + F8; % Ratio of maintenance man-hours to flight hours
oper.maintenance = F_MF * R_AP * Q_FLGT;

% Storage cost (per year)

oper.storage = 12 * R_STOR;

% Cost of fuel (per year)

oper.fuel = BHP_cruise .* SFC_cruise * Q_FLGT * R_fuel ./ 6.5;

% Insurance cost (per year)
C_AC = dev.cost_per_aircraft(7);  % Insured value of the aircraft (10 aircraft/month)
oper.insurance = 500 + 0.015 * C_AC;

% Inspection cost (per year)

oper.inspection = 500;

% Engine overhaul fund (per year)

oper.fund = 5 * N_pp .* Q_FLGT;

% Electricity cost (per year)

oper.electricity = energy_cons_battery/3.6e6 * price_kWh * route_day * oper_days;

% Batteries replacement cost (per year)

oper.batteries = batteries * route_day * oper_days /1000; % Change every 1000 cycles (routes)

% Crew cost (per year)

oper.crew = R_crew * Q_FLGT;

% Total cost per year

oper.tot = oper.maintenance + oper.storage + oper.fuel + oper.insurance + oper.inspection + oper.fund + oper.electricity + oper.crew + oper.batteries;

% Cost per flight hour

oper.flight_hour = oper.tot/Q_FLGT;

%% Example LA -> San Fransisco (300 nmi)

flight_duration = 1.84; %[h]
charge = 1/3; %[h] Charge of batteries + boarding

LA_SF = flight_duration * oper.flight_hour; %[$]

%% LA -> San Diego

f_d = 0.65; %[h]

LA_SD = f_d * oper.flight_hour; %[$]

%% NY -> Washington

f_d_2 = 1.13; %[h]

NY_W = f_d_2 * oper.flight_hour; %[$]

%% Seattle -> Portland

f_d_3 = 0.83; %[h]

S_P = f_d_3 * oper.flight_hour; %[$]

%% Left to do
% Compare with others tranports (flexibus, airlines, train,...)(distance travelled, price, time)
% Make a sheet with all the values 
% Non-recuring and recuring costs -> production cost, operationnal cost
% Break-even graph
% Pie graph of each contribution

%% Break-even analysis 
%variable_costs(1) = dev.variable_cost(7)/600;
sell_2 = 900000;
%sell_3 = 1.1e6;
for i=1:600
    variable_costs(i) = i*dev.variable_cost(7)/600;
    profit(i) = i*selling_price(7) - dev.certification(7) - i*(dev.variable_cost(7)/600);
    profit_2(i) = i*sell_2 - dev.certification(7) - i*(dev.variable_cost(7)/600);
    %profit_3(i) = i*sell_3;
end

total_cost = variable_costs + dev.certification(7);
units_produced = 1:1:600;

figure
plot(units_produced,profit * 10^(-6),units_produced,profit_2 *10^(-6),[0 600],[0 0] * 10^(-6),'-.','LineWidth',1)
hold on 
plot(326.65,0,'o','LineWidth',1)
text(326.65,0,' Break-even (327 units)','VerticalAlignment','top','HorizontalAlignment','left')
hold on
plot(218.5,0,'o','LineWidth',1)
text(218.5,0,'Break-even (219 units)   ','VerticalAlignment','bottom','HorizontalAlignment','right')
grid on
%ylim([0 0.3])
xlabel('Number of units produced [N]','fontsize',14,'interpreter','Latex')
ylabel('Profit [\$M]','fontsize',14,'interpreter','Latex')
lg=legend('$\mathrm{Selling~price} = \$800,000$','$\mathrm{Selling~price} = \$900,000$');
set(lg,'interpreter','Latex')
set(lg,'fontsize',12)
set(lg,'location','northwest')