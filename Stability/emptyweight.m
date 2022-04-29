W_dg=2000;% gross weight lbf
W_ev=zeros(1,100);%lbf
W_ev(1)=W_dg;
%% Units:
kg2lb = 2.2046226218;%Passage de kg à lb
lb2kg = 1/kg2lb;%Passage de lbf en kg
m2ft = 3.28084;%Passage de m en feet
ft2m = 1/m2ft;%Passage de pied en m
m_s = 1.94384;%Passage de m/s à KTAS
Pa = 0.000145038;%Passage de pascal à psi

for i=2:length(W_ev)
%% important parameters
S_w_feet=17.955*10.7639;%ft
W_fw=102.40; %fuel weight in wing lbf
Tot_W=0;
%% Wing parameters

b=12.3*m2ft;%wingspan
A=b^2/S_w_feet; %aspect Ratio
Lambda=0.55;%wing sweep at 25% MGC
lambda=0.55;%taper ratio
c_mean=b/A;%mean chord
%c_root=(3*c_mean/2)/((1+lambda+lambda^2)/(1+lambda));%chord at root
c_root=1.77*3.28;
c_tip=c_root*lambda;
density_aircruise=0.049;% lb/ft^3 0 feet
aircraft_speed=288.714;%feet/sec
q=87.27;%dynamic pressure lbf/ft^2
max_t=0.74; % maximum thickness [feet]
t_over_c=max_t/c_mean;% thickness to cord ratio
n=6.06;%load factor
N_z=1.5*n; %ultimate load factor

%% wing equations
W_wing=0.036*S_w_feet^(0.758)*1*(A/(cos(Lambda)^2))^(0.6)*q^(0.006)*lambda^(0.04)*(100*t_over_c/(cos(Lambda)))^(-0.3)*(N_z*W_ev(i-1))^(0.49);
Tot_W=Tot_W+W_wing;

%% Horizontal tail parameters
% N_z=N_z;
% W_dg=W_dg;
% q=q;
% Lambda=Lambda;
b_ht = 4.380*m2ft;
S_ht=34.74;%horizontal tail area ft^2
A_ht=b_ht^2/S_ht ;%!!!!!! to check: is A the tail one or the main wing one?!!!!!
t_over_c_ht=0.18;%!!!!!! to check: is t/c the tail one or the main wing one?!!!!!
Lambda_htail=0;%wing sweep horizontal tail
lambda_h=0.55;%taper ratio horizontal tail
%% horizontal tail equation
W_horizontal_tail=0.016*(N_z*W_ev(i-1))^(0.414)*q^(0.168)*S_ht^(0.896)*(100*t_over_c_ht/cos(Lambda_htail))^(-0.12)*(A/(cos(Lambda_htail)^2))^(0.043)*lambda_h^(-0.02);
Tot_w=Tot_W+W_horizontal_tail;

%% Vertical tail parameters
% N_z=N_z;
% W_dg=W_dg;
% q=q;
A_ht=2;%!!!!!! to check: is A the tail one or the main wing one?!!!!!
t_over_c_vt=0.18;%!!!!!! to check: is t/c the tail one or the main wing one?!!!!!
ht_over_hv=0;%0 for conventional tail, 1 for T tail
S_vt=16.21;% vertical tail area ft^2
Lambda_vt=0;%wing sweep vertical tail
lambda_vt=1.77/3.23;%tper ratio vertical tail (If lambda_vt is less than 0.2, use 0.2)
%% vertical tail equations
W_vertical_tail=0.073*(1+0.2*ht_over_hv)*(N_z*W_ev(i-1))^(0.376)*q^(0.122)*S_vt^(0.873)*(100*t_over_c_vt/cos(Lambda_vt))^(-0.49)*(A/(cos(Lambda_vt)^(2)))^(0.357)*lambda_vt^(0.039);
Tot_W=Tot_W+W_vertical_tail;

%% fuselage parameters
% N_z=N_z;
% W_dg=W_dg;
% q=q;
l_f=29.53;%fuselage length, ft
w_f=4.95;%fuselage width, ft
h_f=4.95;%fuselage height, ft
d_barre=sqrt(h_f*w_f);%diameter, ft
S_f=pi*d_barre*(l_f-1.3*d_barre);%fuselage wetted area, ft^2
L_t=20;%tail length; wing quarter-MAC to tail quarter-MAC, ft (Hypothesis)
t_f=0.003333;%thickness of the fuselage, ft 
L_over_D=l_f/t_f;%Length over structural depth 
W_press=0;%weight penalty due to pressurization
%% fuselage equation
W_fuselage=0.052*S_f^(1.086)*(N_z*W_ev(i-1))^(0.177)*L_t^(-0.051)*(L_over_D)^(-0.072)*q^(0.241)+W_press;
Tot_W=Tot_W+W_fuselage;

% % main landing gear parameters
% N_l=;%ulitmate landing gear factor =1.5*N_gear
% W_l=0.04*W_dg;%landing design gross weight, lb  (4% of total weight)
% L_m=;%extended length of main landing gear, in.
% % main landing gear equation
% W_main_lgear=0.095*(N_l*W_l)^(0.768)*(L_m/12)^(0.409);
% % nose landing gear parameters
% N_l=;%!!!! same as main lgear?!!!!
% W_l=;% !!!! same as main lgear?!!!!
% L_m=;%!!!! same as main lgear?!!!!
% % nose landing gear equation
% W_nose_lgear=0.125*(N_l*W_l)^(0.566)*(L_m/12)^(0.845);%(reduce total landing gear weight by 1.4% of TOGW if nonretractable)
% Tot_W=Tot_W+W_main_lgear+W_nose_lgear;
Tot_W=Tot_W+0.04*W_ev(i-1); %4% of total weight for landing gear)

%% installed engine parameters

W_installed_engine=220*kg2lb;%lbf
%% installed engine equation
Tot_W=Tot_W+W_installed_engine;

%% Fuel system parameters
V_t=55*1.2346*0.26417;% total fuel volume, gal
V_i=20;% integral tanks volume, gal 
N_t=1;% number of fuel tanks 
N_en=9;
%N_en=N_en;
%% Fuel system equation
W_fuel_system=2.49*V_t^(0.726)*(1/(1+(V_i/V_t)))^(0.363)*N_t^(0.242)*N_en^(0.157);
Tot_W=Tot_W+W_fuel_system;

%% flight control parameters
% N_z=N_z;
% W_dg=W_dg;
L=l_f;%fuselage structural length, ft (excludes radome cowling, tail cap)
B_w=b;%wing span, f
%% flight control equation
W_flight_control=0.053*L^(1.536)*B_w^(0.371)*(N_z*W_ev(i-1)*10e-4)^0.8;
Tot_W=Tot_W+W_flight_control;

%% hydraulics parameters
% W_dg=W_dg;
M_cruise=0.238;
K_h=0.11;%0.05 for low subsonic with hydraulics for brakes and retracts only; 0.11 for medium subsonic with hy?raulics for flaps; = 0.12 for high subsonic with hydraulic flight controls; = 0.013 for light 
M_design=1.06*M_cruise;%Mach number (design maximum) 
%% hydraulics equations
W_hydraulics=K_h*W_ev(i-1)^(0.8)*M_design^(0.5);
Tot_W=Tot_W+W_hydraulics;

%% avionics parameters
W_uav=30;%uninstalled avionics weight, lb (typically = 800-1400 lb) 
%% avionics equation
W_avionics=2.117*W_uav^(0.933);
Tot_W=Tot_W+W_avionics;

%% electrical parameters
%% electrical equation
W_electrical=12.57*(W_fuel_system+W_avionics)^(0.51);
Tot_W=Tot_W+W_electrical;

%% air conditionning and anti-ice parameters
% W_dg=W_dg;
% M=M;
N_p=4;%number of personnel onboard (crew and passengers
%% air conditionning and anti-ice equation
W_acaai=0.265*W_ev(i-1)^(0.52)*N_p^(0.68)*W_avionics^(0.17)*M_design^(0.08);
Tot_W=Tot_W+W_acaai;

%% furnishings equations
W_furnishings=120 %0.0582*W_ev(i-1)-65;
Tot_W=Tot_W+W_furnishings;



W_ev(i)=Tot_W; %lbf

end
W_kg=W_ev(end)*lb2kg

kg2lb = 2.2046226218;%Passage de kg à lb
lb2kg = 1/kg2lb;%Passage de lbf en kg
m2ft = 3.28084;%Passage de m en feet
ft2m = 1/m2ft;%Passage de pied en m
m_s = 1.94384;%Passage de m/s à KTAS
Pa = 0.000145038;%Passage de pascal à psi

W_wing = lb2kg*W_wing
W_flight_control = lb2kg*W_flight_control
W_horizontal_tail = lb2kg*W_horizontal_tail
W_hydraulics = lb2kg*W_hydraulics
W_vertical_tail = lb2kg*W_vertical_tail
W_avionics = lb2kg*W_avionics 
W_electrical = lb2kg*W_electrical 
W_fuselage = lb2kg*W_fuselage
W_fuel_system = lb2kg*W_fuel_system
W_acaai =lb2kg*W_acaai
W_furnishings = lb2kg*W_furnishings
W_installed_engine = lb2kg*W_installed_engine
W_landing_gear = lb2kg*W_landing_gear













