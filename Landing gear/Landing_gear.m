%% Landing gear
close all; clear; clc;

g = 9.81; %[m/s^2]
w_tot = 1600; %[kg]
cg_dist = 3.11; %[m]
cg_height = 0.5; %[m]
fuselage_length = 9.45; %[m]
tail_height = 0.85; %[m]
fuselage_D = 1.5; %[m] diameter of the fuselage 

W_main = (0.9*w_tot)/2; %[kg] 90% of the mass on the main landing gear

W_nose = 0.1*w_tot; %[kg]

%% wheel diameter
A=5.1;
B=0.349;

d_main_1= A*W_main^B*10^(-2); %[m]
d_nose_1= A*W_nose^B*10^(-2); %[m]

rad_main_1 = d_main_1/2;
rad_nose_1 = d_nose_1/2;

%% wheel width
A=2.3;
B=0.312;

wi_main_1=A*W_main^B*10^(-2); %[m]

wi_nose_1=A*W_nose^B*10^(-2); %[m]

%% New dimensions from Raymer table : type 3, 7-8

% First assumption
wi_main=0.185; %[m]
d_main=0.53; %[m]
d_main_wheel=0.2; %[m]
rol_rad_main=0.2109; %[m] rolling radius

wi_nose= 0.7 * wi_main; %[m]
d_nose= 0.7 * d_main; %[m]

% Choice in the
% http://www.goodyearaviation.com/resources/pdf/2021%20Data%20Section.pdf
% in order to support the load
% Main: 606C81B1

wi_main = 0.16; %[m]
d_main = 0.44; %[m]
rol_rad_main = 0.175; %[m]
d_rim_main = 0.152; %[m]

% Nose: 505C41-4

wi_nose = 0.125; %[m]
d_nose = 0.35; %[m]
rol_rad_nose = 0.145; %[m]
d_rim_nose = 0.127; %[m]

%% Stroke determination

V_vert=3; %[m/s]
S_T= (d_main/2 - rol_rad_main); %[m] tire stroke
eta=0.7; %[-] shock absorber efficiency
eta_T=0.47;
N_gear=3; %[m/s^2] vertical deccelaration rate

S_1 =  (V_vert^2)/(2*g*eta*N_gear)- (eta_T/eta)*S_T; %[m] stroke dimension
S_1=S_1+0.03; %[m] security

%% OLEO

static_dist=0.5*S_1;
oleo_length_main = 2*static_dist;
oleo_length_nose = oleo_length_main+(d_main - d_nose);

length_tot_gear_main= oleo_length_main+d_main; %[m] total length of main gear
length_tot_gear_nose= oleo_length_nose+d_nose; %[m] total length of main gear

%% Tip back angle

Mf_B=0.2; %[-]

H_1= length_tot_gear_main + cg_height; %[m] cg to ground

h_1= length_tot_gear_main + tail_height; %[m] tail to ground

alpha=14.5; %[°] angle imposed tipback angle

l_1=h_1/(tan(alpha*pi/180)); %[m] horizontal distance btw main gear and tail

M_f= fuselage_length - cg_dist - l_1; %[m] horizontal distance btw cg and main gear

theta_1 = atan(M_f/H_1); %[rad] angle btw cg and main wheel >tipback angle

theta_1= theta_1 * 180/pi; %[°]

B_1= M_f/Mf_B;
N_f= B_1 - M_f; %[m] horiz distance btw nose gear and cg

fus_nose_gear_nose= cg_dist - N_f;
fus_nose_gear_main= cg_dist + M_f;

max_stat_load_main= w_tot * 0.9/2; %static load on main gear
max_stat_load_nose= w_tot * 0.1; %static load on nose tire

decc=4; %[m/s^2] dynamic load on nose gear

dynamic_break_load_nose= 5* H_1* w_tot/(g*B_1);

%% OLEO next

L_oleo_main = max_stat_load_main / 0.45359237; %[lb] load carry by main oleo
L_oleo_nose = (max_stat_load_nose+dynamic_break_load_nose) / 0.45359237; %[lb] '' nose oleo

P=1800; %[psi]

D_oleo_main= 1.3 * sqrt((4*L_oleo_main)/(P*pi)); %[in]
D_oleo_nose= 1.3 * sqrt((4*L_oleo_nose)/(P*pi)); %[in]

D_oleo_main = D_oleo_main * 2.54 * 10^(-2); %[m]
D_oleo_nose = D_oleo_nose* 2.54 * 10^(-2); %[m]

%% Overturn angle 
Wheel_base = N_f+M_f;
omega = 60* pi/180; %[rad]
alpha_track = asin(H_1/(N_f*tan(omega)));
Wheel_track = 2 * Wheel_base * tan(alpha_track);

%% 5° roll 

prop_R = 0.9; %[m] radius of the wingtip propeller 
span = 12.3/2; %[m] semi-span
h_wing = length_tot_gear_main + fuselage_D; %[m] height of the wing
h_tip = h_wing - span * sin(5*pi/180); %[m] wing tip height > 0.95 because of the propeller

%% Results

res.wheel_base = Wheel_base * 3.2808; %[ft]
res.wheel_track = Wheel_track * 3.2808; %[ft]
res.H_cg = H_1 * 3.2808; %[ft]
res.M_f = M_f * 3.2808; %[ft]
res.N_f = N_f * 3.2808; %[ft]
res.heigth = length_tot_gear_main * 3.2808; %[ft]
res.alpha = alpha; %[°]
res.theta = theta_1; %[°]
res.alpha_track = alpha_track * 180/pi; %[°]
res.overturn_angle = omega * 180/pi; %[°]
res.length_gear = length_tot_gear_main * 3.2808; %[ft]
res.d_main = 0.44 * 39.37; %[in]
res.w_main = 0.16 * 39.37; %[in]
res.d_rim_main = 0.152 * 39.37; %[in]
res.d_nose = 0.35 * 39.37; %[in]
res.w_nose = 0.125 * 39.37; %[in]
res.d_rim_nose = 0.127 * 39.37; %[in]
res.static_load_main = max_stat_load_main * 2.2046; %[lbf]
res.static_load_nose = max_stat_load_nose * 2.2046; %[lbf]

