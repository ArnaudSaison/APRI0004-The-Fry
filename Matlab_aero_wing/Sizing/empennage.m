clear
clc
close all
% NACA 0012
%% Wing parameter 
S_W = 164.69; %trapezoidal wing surface ft^2
AR_W = 8; %aspect Ratio
b_W = sqrt(AR_W*S_W); %wingspan
Lambda = 0; %wing sweep at 25% MAC
lambda = 0.55; %taper ratio
C_W = b_W/AR_W; %mean chord
AR_VT = 2;
AR_HT = 4;
Sweep_vertical = 20; %°

%% Tail volume coefficient 
L_HT = 20; %ft
L_VT = L_HT;
V_VT = 0.07;
V_HT = 0.6;

%% Conventional
S_VT_conv = V_VT*b_W*S_W/L_VT;
b_VT_conv = sqrt(AR_VT * S_VT_conv);
c_root_VT_conv = b_VT_conv/AR_VT;
c_tip_VT_conv = c_root_VT_conv*lambda; 
S_VT = (c_root_VT_conv+c_tip_VT_conv)*b_VT_conv/2;

S_HT_conv = V_HT*C_W*S_W/L_HT;
b_HT_conv = sqrt(AR_W * S_HT_conv);
c_root_HT_conv = b_HT_conv/AR_HT;
c_tip_HT_conv = c_root_HT_conv*lambda;
S_HT = (c_root_HT_conv+c_tip_HT_conv)*b_HT_conv/2;

%% Dorsal Fin

S_df = 0.18*S_VT;
phi_df = 60; % [°]
phi_VT = 20; % [°]
h_df = sqrt(2*S_df/(tand(phi_df)-tand(phi_VT)));

L_df = h_df/sind(90-phi_df);
c_df = L_df - h_df/sind(90-phi_VT);


%%
dx = 2;
CF = 2*c_root_VT_conv;
x = [0 c_df L_df 0];
y = [0 0 h_df 0];

x2 = [c_df c_df+c_root_VT_conv c_df+4.1525 c_df+2.3737 c_df];
y2 = [0 0 b_VT_conv b_VT_conv 0];

plot(x,y)
hold on
plot(x2,y2)
