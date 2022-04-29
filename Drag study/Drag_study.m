function [CD,CD0,CDi_flow,delta_trim_CD,CL,CD_twist,CD_Vfuselage,CD_HT,CD_skin,CD_p,CD_tail,CD_fus,CD_landing_gear] = Drag_study(AoA,M)
% Torenbeek - annexe F et G
%% Units:
kg2lb = 2.2046226218;%Passage de kg à lb
lb2kg = 1/kg2lb;%Passage de lbf en kg
m2ft = 3.28084;%Passage de m en feet
ft2m = 1/m2ft;%Passage de pied en m
m_s = 1.94384;%Passage de m/s à KTAS
Pa = 0.000145038;%Passage de pascal à psi

[T_alti, a, P, rho_alt] = atmosisa(0*m2ft^(-1))
v_cruise = M*a*m_s;%KTAS calculated in propulsion section
q = 0.5*rho_alt*(1.688*v_cruise)^2;% %page 63
mu = 3.17 * 10^(-11) * T_alti^(1.5) * (734.7/(T_alti+216)); %lb * s / ft²
    
% Wing
Sw = 17.95*m2ft*m2ft; %ft²
bw = 12.3*m2ft;
ARw = bw^2/Sw;
sweep_angle = 2.3; %deg
lambda_w = 0.55;
c_root1 = 1.77*m2ft;
c_tip1 = c_root1*lambda_w;
t_c_w = 0.12;

% VT
c_VT_root = 0.984*m2ft;
c_VT_tip = c_VT_root*lambda_w;
c_VT = (c_VT_tip+c_VT_root)/2;
S_VT = 1.5058*m2ft*m2ft;%ft2
t_c_VT = 0.12;
sweep_angle_VT = 38;%deg

% HT
c_HT_root = 1.02*m2ft;
c_HT_tip = c_HT_root*lambda_w;
c_HT = (c_HT_tip+c_HT_root)/2;
S_HT = 3.227*m2ft*m2ft;%ft2
t_c_HT = 0.12;
sweep_angle_HT = 2;%deg
 
% Fuselage
L_F = 9.455 * m2ft;%ft
w_fus = 1.5 * m2ft;%ft
l_A = 2.958*m2ft ; %length tail
l_N = 0.9777*m2ft ; %length nose

%landing gear
Wheel_track = 1.9271*m2ft ;
Wheel_base = 3.513*m2ft ;
wheel_heigth = 0.959*m2ft ;
main_wheel_diameter = 0.44*m2ft ;
nose_wheel_diameter = 0.35*m2ft ;
width_main_wheels = 0.16*m2ft ;
width_nose_wheel = 0.125*m2ft ;


%% Data
CL0 = 0.015;
Cl_alpha_deg = 0.0812;%deg^-1
Cl_alpha_rad = Cl_alpha_deg*180/pi;%rad^-1
CL_alpha_deg = 0.1217;%deg^-1
CL_alpha_rad = CL_alpha_deg * 180/pi;%rad^-1
%CL = CL_alpha_deg * AoA + CL0;
%CL = 0.3397 ; %Lift cruise 
CL = 6 ; %Take off  
e = 1-(0.25/(cosd(sweep_angle))^2);
b_2 = bw/2;
v_cd = v_cruise %* 1.68781;%ft/s

%Exposed planform area
S_fus = 28.79*m2ft^2;
S_wing = 15.3*m2ft^2;
AR_wing = bw^2/S_wing;

% Wetted surface
S_wing_wetted = 33.27*m2ft^2;
S_fus_wetted = S_fus;

%% Profil drag: (page 543 to 545)
cf_c = 0.49/1.4103; %flaps choice
K_L = 1/(pi*ARw*e);
K_23 = 1.15;

%bf0_down = 56.27;
%bf0_top = 40.19;
bfi = 0.5*m2ft;
%bf_down = bf0_down-bfi;
%bf_top = bf0_top-bfi;
bf = 2*(2.69+2.19)*m2ft ;
%bf = (bf_top+bf_down)/2;
%Swf_down = 69.79;
%Swf_top = 35.96;
%Swf = Swf_top+Swf_down;
Swf = 0.8*S_wing;
Sf = 2.39*m2ft^2;
%Sf_top = 179.81;%ft²
%Sf_down = 348.98;%ft²

delta_f = 0; %degré
theta_f_Cl = acos(2*(cf_c)-1);
alpha_d_Cl = 1-((theta_f_Cl - sin(theta_f_Cl))/pi);

theta_f_CL = acos(2*(Sf/Sw)-1);
alpha_d_CL = 1-((theta_f_CL - sin(theta_f_CL))/pi);

K_b = (2/pi)*((bf/bw)*sqrt(1-(bf/bw)^2)+asin(bf/bw));
F_delta = 0; %fig G21 page545

delta_Cd_p0 = 0.55 * (cf_c) * (cf_c/t_c_w^1.5)^(2/9)*F_delta;
Cd_p0 = delta_Cd_p0/2; 
%delta_Cd_p0 = kd * Cl_alpha * alpha_d * (cf_c) * delta_f * sind(delta_f) + Cd_p0 * ((c_prime/c)-1);

%Swf_S_wing = ((bf0_down-bfi)/b_2)*(1+((1-lambda_w)/(1+lambda_w)) * (1-((bf0_down-bfi)/b_2))); % page 559
delta_f_Cl0 = alpha_d_Cl * Cl_alpha_rad * delta_f*pi/180;
delta_f_CL0 = delta_f_Cl0 * (CL_alpha_rad/Cl_alpha_rad) * (alpha_d_CL/alpha_d_Cl) * K_b;
%G-48
delta_CD_p = K_23 * (Swf/Sw)* delta_Cd_p0 * cosd(sweep_angle) - ...
    K_L * Cd_p0 * delta_f_CL0 * (CL - (CL0+0.25*delta_f_CL0));

%% Vortex-induced drag: (page 545)
Kff = 1/3;
z = (0.07/(1+lambda_w)) * (1-Kff)^2 * (bfi/b_2);
delta_f_Cl = delta_f_Cl0 * (sin(theta_f_Cl)/(pi-(theta_f_Cl-sin(theta_f_Cl))));
w_down = (0.003+0.0025)/2; %fig G-22 page 546
v_down = 0.0005;
w_top = (0.0085+0.009)/2;
v_top = 0.003;
w=w_top+w_down;
v = v_top+v_down;

delta_CD_v = (w+z) * delta_f_Cl^2 + v * CL * delta_f_Cl;

%% Trim drag: (page 545 à 547)
sweep_h = -sweep_angle;
e_h = 1-(0.25/(cosd(sweep_h))^2);
CL_HT = CL/((S_HT*5.729)/(S_wing*0.48));
b_h = bw;
AR_h = b_h^2/S_HT;

delta_trim_CD = ((CL_HT^2/(pi*AR_h*e_h))+CL_HT*sin(e_h)-2*CL*CL_HT/(pi*ARw))*S_HT/Sw; 

%% CD0 (annexe F page 487

% Drag due to twist (p.493)
epsilon_t = -3;%degré
CD_twist = 3.7 * 10^(-5) * epsilon_t^2;

% Vortex induced by fuselage (p.496)
Vf = 86.09*m2ft^3;%ft3
alpha_f = (CL-CL0)/CL_alpha_rad;%rad
CD_Vfuselage = (0.15 * alpha_f^2 * Vf ^(2/3))/Sw;

% HT (page 496)
AR_HT = bw^2/S_HT;
CD_HT = 1.02*(S_HT/Sw)*CL_HT^2/(pi*AR_HT); %F-26 

% CD0 skin friction
phi = 2.7*t_c_w+100*t_c_w^4;
c_wing = (c_root1 + c_tip1)/2;
Re_wing = rho_alt * v_cd * c_wing / mu
Cf_wing = 0.001479; %graph F4 or aerotoolbox
CD_skin = Cf_wing * (1+phi) * S_wing_wetted/S_wing;
CD_skin = CD_skin * 1.09 ; %roughness

% Profil (p.500)
Cl0 = 0.5;
Cl3 = 0.8;
ks = 0.07;
Cl_max = 1.3894;%rad^-1
phi_w = 2.7 * t_c_w + 100 * t_c_w^4; %F-30
Cd_pmin = 2 * Cf_wing * (1+phi_w); %F-29

% Cli = @(y) (Cl3 - (Cl3-Cl0)*y/b_2);
%if Re>10E7
delta_l_Cd_p_ref = 0.01 * Cl_max - 0.0046 *(1+2.75*t_c_w); %(F-32)
% delta_l_Cd_p = @(y) (0.75 * delta_l_Cd_p_ref * ((Cl3 - (Cl3 - (Cl3-Cl0)*y/b_2))/(Cl_max - (Cl3 - (Cl3-Cl0)*y/b_2)))^2);
% Cd_p = @(y) (Cd_pmin + 0.75 * delta_l_Cd_p_ref * ((Cl3 - (Cl3 - (Cl3-Cl0)*y/b_2))/(Cl_max - (Cl3 - (Cl3-Cl0)*y/b_2)))^2);

c = @(y) ((((Sw/2)/((1+lambda_w)*(b_2)))*((1-((1-lambda_w)/(b_2))*y))))...
    *(Cd_pmin + 0.75 * delta_l_Cd_p_ref * ((Cl3 - (Cl3 - (Cl3-Cl0)*y/b_2))...
    /(Cl_max - (Cl3 - (Cl3-Cl0)*y/b_2)))^2);

CD_p = (2/Sw) * (integral(c,w_fus/2,b_2)); %F-34

% Fuselage and tail boom (p.501)
A_c = 2.405*m2ft^2;%coupe transversale du fusalge
Df_eff = sqrt((4/pi)*A_c);
lambda_eff1 = L_F/Df_eff;
lambda_eff2 = ((l_N+l_A)/Df_eff) + 2;
if lambda_eff1 < lambda_eff2
    lambda_eff = lambda_eff1;
else
    lambda_eff = lambda_eff2;
end
phi_f = (2.2/lambda_eff^1.5) + (3.8 / (lambda_eff^3));
Re_fus = rho_alt * v_cd * L_F / mu
Cf_fus = 0.0012; %F-4aerotoolbox
CD_S_basic = Cf_fus * S_fus_wetted * (1+phi_f);
L_cyl = L_F-l_A-l_N ;
L_front = l_N ;
L_back = l_A ;
A_I = L_cyl * w_fus + L_front * w_fus/2; %ft² voir fig F-13 page 505
A_II = L_back * w_fus/2; %ft² voir fig F-13 page 505
beta_fus = 15*pi/180; %boat tail angle
delta_ab_CD_S = A_I * abs((sin(alpha_f))^3) + A_II * abs((sin(alpha_f-beta_fus))^3)/cos(beta_fus);
delta_ab_CD_S = 0 ;
CD_fus = (CD_S_basic + delta_ab_CD_S)/S_fus 

% Tail drag:
Re_VT = rho_alt * v_cd * c_VT / mu 
Cf_VT = 0.0014;
CD_VT = 2*Cf_VT * (1+2.75*t_c_VT*(cosd(sweep_angle_VT))^2)*S_VT/Sw;

Re_HT = rho_alt * v_cd * c_HT / mu 
Cf_HT = 0.0014;
CD_HT_basic = 2*Cf_HT*(1+2.75*t_c_w*(cosd(sweep_angle))^2)*S_HT/Sw;
delta_l_CD_HT = 0.33*CL_HT^2/((cos(sweep_angle))^2*pi*AR_HT)*S_HT/Sw;
CD_tail = CD_VT + CD_HT_basic + delta_l_CD_HT;

% CD0
CD0 = CD_twist + CD_Vfuselage + CD_HT + CD_skin + CD_p +  + CD_tail + CD_fus;

%% Interference correction
% Wing/fuselage
nf = w_fus/bw;
% Vortex induced (p.509)
delta_i_CDV = (0.55*nf)/(1+lambda_w) * (2-pi*nf)*CL0^2/(pi*ARw); %F-65

% Viscous interference (p.509-510)
Cci = 4.5*c_root1;
tr = t_c_w*c_root1;
delta_i_CDviscous = 1.5 * Cf_wing * tr * Cci * (cosd(sweep_angle))^2/Sw; %F-66

% Nacelle: %F-70 (p.510)
r_nacelle = 0.257*m2ft/2; %petites nacelles
delta_i_nacelle = 12*0.008*r_nacelle^2*pi/Sw;
r_nacelle = 0.406*m2ft*2; %grosses nacelles 
delta_i_nacelle = delta_i_nacelle + 2*0.008*r_nacelle^2*pi/Sw;
%landing gear 
CD_pi_main = 0.51 ;
CD_pi_nose = 0.42 ;
main_bD = 2*main_wheel_diameter*width_main_wheels ; 
bD = nose_wheel_diameter*width_nose_wheel ;

CD_landing_gear = CD_pi_main/(rho_alt*main_bD*(v_cruise^2)/2)+CD_pi_nose/(rho_alt*bD*(v_cruise^2)/2) ;
% CD interference
CD_interf = delta_i_CDV + delta_i_CDviscous + delta_i_nacelle + CD_landing_gear;


%% CDi:
CDi = CL^2/(pi*ARw*e);
CDi_flow = 0.0082;
 
%% Total drag:
CD = CD0 + delta_CD_p + delta_CD_v + delta_trim_CD + CDi+CD_interf; %p.543 G-41

D = 0.5 * rho_alt * CD * Sw * v_cd^2;

%% Breakdown:
Vortex_drag = delta_CD_v + CD_twist + CD_Vfuselage+CD_HT+delta_trim_CD;
perc_vortex = Vortex_drag/CD
Profil_drag = delta_CD_p + CD_skin + CD_p + CD_fus + CD_tail;
perc_profil = Profil_drag/CD
Induced_drag = CDi;
perc_induced = Induced_drag/CD
perc_interf = CD_interf/CD
perc_tot = perc_vortex + perc_profil + perc_induced + perc_interf;
end