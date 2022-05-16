%% drag study
%% This code comes largely from the Aerialeater based their drag study on the Gundmundsson, 
%% thanks to them for allowing us to make a better drag study at the AAIA 
% Call coefficients before
% clear all

%% 

% sectional lift, chord & angle of attack distribution, length and surface
% wing
% t/c & x/c & sweep angle for wing and tails
% Wetted and projected surfaces
% Cruise speed

%% Conversion
InToFtsq = 0.08333333^2;
knot2m = 0.51444444;
Density_kgmToLbsFt = 0.062427960576145;

kg2lb = 2.2046226218;%Passage de kg à lb
lb2kg = 1/kg2lb;%Passage de lbf en kg
m2ft = 3.28084;%Passage de m en feet
inch2m = 0.0254;
ft2m = 1/m2ft;%Passage de pied en m
m_s = 1.94384;%Passage de m/s à KTAS
Pa = 0.000145038;%Passage de pascal à psi

%% Données
% Weight
W = 1600; 
cl = [0.3411];
ClEx = [0 1];
DragToLiftRatio = cl ;% linspace(ClEx(1), ClEx(2),10000);
% Speeds 
alttarget = 13000*ft2m;
Vl = linspace(86,86,1);
altl = linspace(alttarget, alttarget,1);

% Coefficients
Tr = 371;                                 % Thrust thermal                [lbf] 

% Need AOA and Cl for landing and takeoff

for i=1:length(Vl)
for j=1:length(altl)
for m=1:length(DragToLiftRatio)
[T,c,p,rho] = atmosisa(altl(j));
V = Vl(i);
Mach = V/c;
Cl = DragToLiftRatio(m);
AR = 8;
%Cl = 0.38;
%% Wetted and reference surfaces :
% 
% tail v
S_wetted_tailv = 2*1.5058*m2ft*m2ft;        % [ft²]
S_projected_tailv = 16.21;      % [ft²] might need to be computed from side

% tail h
S_wetted_tailh = 2*3.227*m2ft*m2ft;               % [ft²]
S_projected_tailh = 34.74;      % [ft²] might need to be computed from side

% wing:
S_wetted_wing =  33.27*m2ft^2;                 % [ft²]
Sref = 17.95*m2ft*m2ft;                  % [ft²] 

% fuselage:
S_wetted_fus = 28.79*m2ft^2;                   % [ft²]
S_projected_fus = 2.1*m2ft^2;      % [ft²]

% blown wing prop (for each, total of 12):
S_wetted_blown = 6909.16*InToFtsq;       % [ft²]
S_projected_blown = 144.95*InToFtsq;     % [ft²] 

% thermal prop:
S_wetted_prop = 1138.64*InToFtsq;        % [ft²]
S_projected_prop = 373.29*InToFtsq;      % [ft²] 

% main landing gear (for each, total of 2):
S_wetted_MLG = 18.7;         % [ft²] 
S_projected_MLG = 255.84*InToFtsq;       % [ft²] 

% nose landing gear: 
S_wetted_NLG = 1695.17*InToFtsq;         % [ft²]
S_projected_NLG = 267.08*InToFtsq;       % [ft²]

% In a vector
S_wetted  = [S_wetted_tailv S_wetted_tailh S_wetted_wing S_wetted_fus S_wetted_blown S_wetted_prop S_wetted_MLG S_wetted_NLG];
S_projected  = [S_projected_tailv S_projected_tailh Sref S_projected_fus S_projected_blown S_projected_prop S_projected_MLG S_projected_NLG];

%% Reynolds number
mu = 1.458*10^(-6)*(T)^(1.5)*(1/(T+110.4));   % [] Viscosity dynamic

% Characteristics length
C_th = 2.5935*ft2m;   % [m]
C_wing  = 4.5005*ft2m ;  % [m]
C_fus   = 9.455;   % [m]        / par mean width 
C_blown = 0.8;     % [m]
C_tv    = 2.5*ft2m;       % [m]
C_MLG   = 1.14;    % [m]        
C_NLG   = 1.14;    % [m]

C = [C_tv C_th C_wing C_fus C_blown  C_MLG C_NLG];
Re = rho*V*C/mu;
k = 5*10^(-6);
Re_cutoff = 38.21*(C/k).^(1.053);

%% Laminar et Turbulent
Cf_lam  = 1.328./sqrt(Re);
Cf_turb = 0.455./((log10(Re).^(2.58))*(1+0.144*Mach^2)^(0.65));
%% Skin friction : Rapport turbulent & laminar
% From table 12.4 Raymer
a_lam = 0.10;   % Plane surface : wings & tail
b_lam= 0.05;    % Rest of the aircraft

C_f = zeros(size(Cf_lam));
C_f(1:3) = a_lam*Cf_lam(1:3) + (1-a_lam)*Cf_turb(1:3);
C_f(4:end) = b_lam.*Cf_lam(4:end) + (1-b_lam).*Cf_turb(4:end);

%% Form and Interference factor :  Wing, Tail (h & t), Fuselage
% Wing
t_cw = 0.12;             % airfoil thickness-to-chord length ratio
x_c_meanw = 0.298;         % chordwise location of the airfoil maximum thickness point.
Sweep_meanw = 0;         % 
FF_w = (1+0.6*t_cw/x_c_meanw+ 100*t_cw^4)*(1.34*Mach^0.18*(cos(Sweep_meanw)^0.28));
IF_w = 1;

% Tails OK
t_cv = 0.12;          % airfoil thickness-to-chord length ratio
t_ch = 0.12;          % airfoil thickness-to-chord length ratio
x_c_meanv = 0.3;     % chordwise location of the airfoil maximum thickness point.
x_c_meanh = 0.3;
Sweep_meanv = 20*pi/180;
Sweep_meanh = 2.3*pi/180;
FF_th = (1+0.6*t_ch/x_c_meanh+ 100*t_ch^4)*(1.34*Mach^0.18*(cos(Sweep_meanh)^0.28));
FF_tv = (1+0.6*t_cv/x_c_meanv+ 100*t_cv^4)*(1.34*Mach^0.18*(cos(Sweep_meanv)^0.28));
IF_t = 1.04;

% Fuselage OK
A_maxfus = 2.1;
l_fus = 9.455;
f_fus = l_fus/sqrt(4/pi*A_maxfus);

FF_fus= (0.9 + 5/((f_fus)^1.5) + f_fus/400);
IF_fus = 1;

% Resume : Tail, Wing & Fuselage
IF = [IF_t IF_t IF_w IF_fus];
FF = [FF_tv FF_th FF_w FF_fus];


%% CD0 parasite drag : Wing, tails & fuselage
CD0 = S_wetted(1:4).*C_f(1:4).*IF.*FF/Sref;

%engines 
l_litle_engine = 0.7*m2ft ;
l_big_engine = 1.3*m2ft ;
S_wet_litle = (0.5766+0.1098)*m2ft^2 ;
S_wet_big = (1.851+0.5566)*m2ft^2 ;
Re_litle_engine = rho * V * l_litle_engine*ft2m / mu ;
Re_big_engine = rho * V * l_big_engine*ft2m / mu ;
lambda_litle_engine = (0.34*m2ft+(C_wing))/(0.239*m2ft);
lambda_big_engine = (l_big_engine)/(0.406*m2ft);
Cf_le = 0.05*1.328/sqrt(Re_litle_engine)+0.95*0.455./((log10(Re_litle_engine).^(2.58))*(1+0.144*Mach^2)^(0.65)) ;
Cf_be = 0.05*1.328/sqrt(Re_big_engine)+0.95*0.455./((log10(Re_big_engine).^(2.58))*(1+0.144*Mach^2)^(0.65)) ;
CD_skin_engines = 12*Cf_le*(1+2.2/lambda_litle_engine^1.5+3.8/lambda_litle_engine^3)*S_wet_litle/Sref+2*Cf_be*(1+2.2/lambda_big_engine^1.5+3.8/lambda_big_engine^3)*S_wet_big/Sref ;

CD0 = [CD0, CD_skin_engines];
Xaxis = categorical({'V tail', 'H tail', 'Wings', 'Fuselage', 'Nacelles'});
%Xaxis = reordercats(Xaxis,{'V tail', 'H tail','Nacelles' , 'Fuselage', 'Wings'});
bar(Xaxis,CD0);
CD0_sum = sum(CD0) ;


%% Variation of Dmin : Engines and landing gear
% CDmin added at the end see torrenbek
% Nacelle of thermal engine 
%{
A_maxth = pi*(0.9)^2;                              % about 1.8/2
A_maxelec = pi*(0.5)^2;

DCDmin_elec = 0.5*A_maxelec/S_projected_blown ;  % 10 electric propellers, one thermal
DCDmin_th = 0.25*A_maxth/S_projected_prop  ;     % one thermal


% Landing gear (724 graph)
wheelh = 0.3861;
hightfus =  1.06;
widthwheel= 0.17530;
NoseToLG = 1.3661;

THR = hightfus/wheelh;
curve = NoseToLG/wheelh;
CD_NLG = 0.7*wheelh/widthwheel/S_projected_NLG;
% Landing gear (table)
DCDmin_MLG = 0.017*Sref/S_projected_MLG*2;                 % Open/!\ there are two 
DCDmin_NLG = 0.017*Sref/S_projected_NLG;                   % Open
%}
% tires
d_main = 17.3*inch2m*m2ft;   % [ft]
w_main = 6.3*inch2m*m2ft;   % [ft] 
d_nose = 13.8*inch2m*m2ft;   % [ft]
w_nose = 4.9*inch2m*m2ft;     % [ft]
dw_nose = d_nose*w_nose;
dw_main = d_main*w_main;

% Fairing
H_fair = 0.329*m2ft;   % [ft]
W_fair = 0.235*m2ft;   % [ft]
HW_fair = H_fair*W_fair;


dCDs_nose  = dw_nose*0.242/Sref; % -> we use this
dCDs_main  = dw_main*0.458/Sref; % -> we use this : already multiplyed by 2
mm=1;

%% Miscellaneous
% Trim
% /!\ h can have different 
W_lb = W/lb2kg;                            % weight at condition           [lbs]
k = 1/(pi*AR*0.77);                          % lift-induced drag constant    [-]
q = 1/2*rho/515*(V/knot2m*1.69)^2;   % dynamic pressure              [lb/ft²]
CMGC = C_wing ;           % wing mean geometric chord     [ft]
l_HT = 8.6+2.5935/4-(2.235+C_wing/4);       % distance between the wing AC  [ft]
                                            % and the C/4 of the HT         
h = 3.11-2.235;                          % dis btw the wing LE and CG    [ft]
h_AC = CMGC/4;                              % dis btw the wing LE and AC    [ft]
C_mw = -0.3411*(3.11-(2.235+C_wing/4));

z = 0.592;
M_W = q*Sref*CMGC*C_mw;

dCd_trim = (k/((q*Sref)^2))  *  ((1+(h-h_AC)/l_HT))^2  *  (W - (M_W - Tr*z)/(h_AC-h+l_HT))^2 ...
            - (k*W^2/(q*Sref)^2);


% Screen shields
dCD_screen =  0.002*S_projected_fus/Sref;   % Curved windscreen with 
                                            % a round upper edge
                                            
% No compressiblity effect cause mach <0.6 see graph p730 TORENBOOK

% Washout
phi_tip = deg2rad(2);                       % 
phi_MGC = 0;                                % To update
dCd_washout = 0.00004*(phi_tip - phi_MGC);


%% Leakage and protuberance drag : antennas and fails of facturing
% Can be added in the misc part.            RAYMER p430

% For a propeller aircraft, these protuberances increase of 5-10% the drag.
% and thus the drag coefficient. This increase must be added at the end of
% the computation. 

C_leak = 0.1;                % To multiply the total drag of coefficient 

CD_misc = (dCd_washout + dCD_screen + dCd_trim + dCDs_main + dCDs_nose)*(1+C_leak);
% Wings washout, antennas, windshield, trim, LGs
%% Induce drag from lift : p.687 General aviation method
delta = 0.11;                                             % Fig 15-22
CDi_wing(m) = Cl^2/pi/AR*(1+delta);                        % P.689

% To update then ok
%% ALL CD
CRUD = 1.25;
CD(m) = (CD0_sum + CD_misc)*CRUD + CDi_wing(m);
D = max(rho*Sref*ft2m^2*(V)^2*CD(m)/2);
% if (abs(Cl-0.3783) < 0.5*10^(-4) || abs(Cl-0.3487) < 0.5*10^(-4) )
%     plot(CD(m),Cl,'o')
%    
% end

end 
end
end


% fig.CLCLCD = figure(1);
% hold on 
% plot(DragToLiftRatio,DragToLiftRatio./CD,'LineWidth',1.5)
% xlabel('$CL$  $[-]$','Interpreter','latex')
% ylabel('$CL/CD$  $[-]$','Interpreter','latex')
% grid minor
% fig2pdf(fig.CLCLCD, 'CLCLCD', 1, 1.5, '.Figures/')
% hold off

% fig.CDCL = figure(1);
% hold on 
% plot(CD,DragToLiftRatio,'LineWidth',1.5)
% xlabel('$CD$  $[-]$','Interpreter','latex')
% ylabel('$CL$  $[-]$','Interpreter','latex')
% grid minor
% fig2pdf(fig.CDCL, 'CDCL', 1, 1.5, '.Figures/')
% hold off

%% Landing and take-off configuration 
% Flaps of 30° of deflection p727
dCD_flaps_TO  = 0.02;
dCD_flaps_L  = 0.05;

% the altitude & velocity are != so we need to update rho
% We also use the propellers electrical so we have z!= 0 and 

% Lift induced drag 
[T,c,p,rho_TOL] = atmosisa(0);
V_TO = 18.5;                                      % TO end ground run
V_L  = 18.5;                                       % Landing speed

Cl_TO = 2.8;                                     % End ground run   
%Cl_TO = 7.51;                                       % end turn
Cl_L  = 8.35;
CDi_wing_TO = Cl_TO^2/pi/8*(1+delta);              %                        P.689
CDi_wing_L = Cl_L^2/pi/8*(1+delta);                %                        P.689

% Trim drag
q_TO = 1/2*rho_TOL/515*(V_TO/knot2m*1.69)^2;    % dynamic pressure       [lb/ft²]
q_L = 1/2*rho_TOL/515*(V_L/knot2m*1.69)^2;    % dynamic pressure       [lb/ft²]

z = 0.7585/ft2m;                                   % COG-Thrustelec         [ft]
Tr = 247;                                          % Thrust of elec         [lbf]
dCd_trim_TO = (k/((q_TO*Sref)^2))  *  ...
              ((1+(h-h_AC)/l_HT))^2 * ...
              (W - (M_W-Tr*z)/(h_AC-h+l_HT))^2 ...
               - (k*W^2/(q_TO*Sref)^2);            % The trim drag is lower close to the ground : validated

dCd_trim_L = (k/((q_L*Sref)^2))  *  ...
              ((1+(h-h_AC)/l_HT))^2 * ...
              (W - (M_W-Tr*z)/(h_AC-h+l_HT))^2 ...
               - (k*W^2/(q_TO*Sref)^2);            % The trim drag is lower close to the ground : validated
           
% Total & miscellaneous     
CD_misc_TO = (dCd_washout + dCD_screen + dCd_trim_TO + dCDs_main + dCDs_nose + dCD_flaps_TO)*(1+0.075);
CD_misc_L = (dCd_washout + dCD_screen + dCd_trim_L + dCDs_main + dCDs_nose + dCD_flaps_L)*(1+0.075);

CD_TO = (CD0_sum + CD_misc_TO)*CRUD + CDi_wing_TO ;
CD_L = (CD0_sum + CD_misc_L)*CRUD + CDi_wing_L ;

% D_TO = rho_TO*Sref*FtoM^2*(V_TO)^2*CD_TO/2;      % Might be useless since
                                                   % variation huge for landing and take-off speeds

stop = 1;