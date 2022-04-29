clear
clc

%% DATA
NACA = 23012;   % Airfoil
AR = 8;         % Aspect Ratio [-]
lambda = 0.55;  % Taper Ratio [-]
MAC = 4.63;     % Mean Aerodynamic Chord [ft]
S = 164.69;     % Area [ft²]
C_root = 5.81;  % Root Chord [ft]
C_tip = 3.20;   % Tip Chord [ft]
b = 18.14*2;      % Wingspan [ft] 
AoI = 4;        % Angle Of Incidence [°]
CD0 = 0.029;    % Zero-Lift Drag Coefficient [-]
e = 0.83;       % Oswald's efficiency factor [-]
cl = 0.355;     % Lift Coefficient @ Cruise [-]
cd = 0.0354;    % Drag Coefficient @ Cruise [-]
AoA = 0;        % Angle Of Attack @ cruise [°]
LAMBDA = 0;     % Sweep Angle [°]  

%% CL 3D approximation page 346 gundmunsson
% XFOIL DATA
polar1 = importdata('polar_NACA_23012.txt');
polar = polar1.data;

aoa_xfoil = polar(:,1);
cl_xfoil = polar(:,2);
aoa_vec = aoa_xfoil(1:30);
cl_vec = cl_xfoil(1:30);

% p = polyfit(aoa_vec, cl_vec,1);
  
% cl_alpha = p(1); % 2D Lift Curve Slope
cl_alpha = (cl_vec(3)-cl_vec(2))/(aoa_vec(3)-aoa_vec(2));
aoa_0 = aoa_vec==0;
cl_0 = cl_vec(aoa_0);     % Lift @ 0 AoA (from XFOIL)

cruise_speed_kn = 170;                   % [knots]
cruise_speed = cruise_speed_kn*1.68781;  % [ft/s]
speed_sound = 1115.49;                   % [ft/s]
Mach_number = cruise_speed/speed_sound; 
beta = sqrt(1-Mach_number^2);     % Mach number parameter
kapa = cl_alpha*(180/pi)/(2*pi);  % Ratio 2D lift curve slope to 2pi

% CL_alpha = AR * cl_alpha /  ...
%            (2 + sqrt(AR^2*beta^2/kapa^2 * (1 + (tan(LAMBDA)/beta)^2))+4); 
                             % 3D Lift Curve Slope
CL_alpha = 0.995*cl_alpha/((1+(2*.55)/(8*(1+0.55)))+(cl_alpha/(pi*8)));                           
% CL_alpha = CL_alpha_per_rad / (180/pi);                               
alpha_0 = -cl_0/cl_alpha;    % 0 Lift angle [°]
CL_0 = -alpha_0 * CL_alpha;  % 3D Lift @ 0 AoA

% Plot the CL in fct of alpha curve 
m = CL_alpha;
y1 = m*-5 + CL_0;
y2 = m*10 + CL_0;
hold on
line([-5 10],[y1 y2],'LineWidth',1.5)
line([2 2],[-2 5],'LineWidth',1.5, 'LineStyle', '--')
plot(2,2*CL_alpha, 'ok','LineWidth', 1.5)
xlabel('Angle of Attack [$^\circ$]', 'Interpreter','latex', 'FontSize', 20);
ylabel('CL$_\alpha$ [-]', 'Interpreter','latex', 'FontSize', 20);
ylim([-0.5 2.5])
grid on
grid minor
hold off

fprintf('Approximation for 3D Lift Coefficient @ 0° AoI : %f\n', (0)*CL_alpha+CL_0);
fprintf('Approximation for 3D Lift Coefficient @ 1° AoI : %f\n', (1)*CL_alpha+CL_0);
fprintf('Approximation for 3D Lift Coefficient @ 2° AoI : %f\n', (2)*CL_alpha+CL_0);
fprintf('Approximation for 3D Lift Coefficient @ 3° AoI : %f\n', (3)*CL_alpha+CL_0);
fprintf('Approximation for 3D Lift Coefficient @ 4° AoI : %f\n', (4)*CL_alpha+CL_0);

%% CL,max 3D approximation
Re_approx = 6400 * MAC * cruise_speed;
cl_max_root = 1.8839;  % From XFOIL using Re_approx [-]
cl_max_tip = cl_max_root*lambda;
aoa_stall = 19.1;      % From XFOIL using Re_approx [°]

MGC = 2/3 * C_root * (1+ lambda + lambda^2) / (1 + lambda);
cl_max = cl_max_root + 2*MGC/b*(cl_max_tip-cl_max_root);
fprintf('Approximation for 2D cl_max @ cruise : %f\n', cl_max);

CL_max = 0.9 * cl_max;
fprintf('Approximation for 3D CL max @ cruise : %f\n', CL_max);

%%

S_ = 15.3; % surface
AR = 8; % aspect ratio
V = 88.0; % speed
[~,~,~,rho_0] = atmosisa(0); % air density at sea level


Di = [48.630993 10.680111 14.081778 58.907230 145.164777 272.799804 ...
      441.694808 651.669474 902.480797 1193.823238 1525.328927 ...
      1896.567909 2307.048436 2756.217304 3243.460241 3768.102346 ...
      4329.408584]; 

F_z = [-10758.642421 -4880.722894 1020.851058 6940.875528 12874.133487 ...
       18815.397365 24759.431634 30700.995393 36634.844959 42555.736463 ...
       48458.428451 54337.684496 60188.275813 66004.983890 71782.603122 ...
       77515.943464 83199.833093]; 

C_L = F_z/(0.5*rho_0*V^2*S_);
C_Di = Di/(0.5*rho_0*V^2*S_);
e = C_L.^2/pi/AR/C_Di;

C_L_alpha = (C_L(6)-C_L(2))/(6-2) * 180/pi;

%% comp vlm analytics

AoA_plot = -6:10;

CL_vlm_analytic = figure(1);
hold on
AoI = 4;
CL_analytic = CL_alpha*((-6:10)+AoI) + CL_0;
plot(AoA_plot,CL_analytic,'LineWidth',1.5)
plot(AoA_plot,C_L,'Color',[0.8500 0.3250 0.0980], 'LineWidth',1.5)
lgd2 = legend('Analytical','VLM');
lgd2.Interpreter='latex';
lgd2.Location = 'northwest';
xlabel('$\alpha$ $[^\circ]$','Interpreter','latex')
ylabel('$C_L$  $[-]$','Interpreter','latex')
grid minor
hold on
ylim([-0.5 2])
% fig2pdf(CL_vlm_analytic, 'vlm_analytic_comp', 1, 1.5, '.Figures/')

%% CDi VLM

polar = figure(2);
plot(C_Di,C_L, 'LineWidth',1.5)
% lgd = legend('Analytical','VLM');
% lgd.Interpreter='latex';
% lgd.Location = 'northwest';
xlabel('$C_{Di}$ $[-]$','Interpreter','latex')
ylabel('$C_L$  $[-]$','Interpreter','latex')
grid minor
% fig2pdf(polar, 'drag_polar', 1, 1.5, '.Figures/')


%% CDi VLM vs Analytics
CDi_vlm_analytic = figure(3);
hold on
CDi_analytic = CL_analytic.^2/(pi*AR*e);
plot(AoA_plot,CDi_analytic,'LineWidth',1.5)
plot(AoA_plot,C_Di,'Color',[0.8500 0.3250 0.0980], 'LineWidth',1.5)
lgd2 = legend('Analytical','VLM');
lgd2.Interpreter='latex';
lgd2.Location = 'northwest';
xlabel('$\alpha$ $[^\circ]$','Interpreter','latex') 
ylabel('$C_{Di}$  $[-]$','Interpreter','latex')
grid minor
fig2pdf(CDi_vlm_analytic, 'CDi_analytic_comp', 1, 1.5, '.Figures/')

%% CDi POLAR VLM vs Analytics
CDi_vlm_POLAR_analytic = figure(3);
hold on

CDi_analytic = (CL_analytic).^2/(pi*AR*e);
plot(CL_analytic,CDi_analytic,'LineWidth',1.5)
plot(C_L,C_Di,'Color',[0.8500 0.3250 0.0980], 'LineWidth',1.5)
lgd2 = legend('Analytical','VLM');
lgd2.Interpreter='latex';
lgd2.Location = 'northwest';
xlabel('$C_{L}$ $[-]$','Interpreter','latex')
ylabel('$C_{Di}$  $[-]$','Interpreter','latex')
grid minor
% fig2pdf(CDi_vlm_POLAR_analytic, 'CDi_POLAR_analytic_comp', 1, 1.5, '.Figures/')

%% 
lift_to_dragi_ratio = figure(5);
L_D = C_L./C_Di;
plot(C_L,L_D, 'LineWidth',1.5)
ylabel('$\frac{C_L}{C_{Di}}$ $[-]$','Interpreter','latex')
xlabel('$C_L$ $[-]$','Interpreter','latex')
grid minor
% fig2pdf(lift_to_dragi_ratio, 'lift_to_dragi_ratio', 1, 1.5, '.Figures/')

%%
Cfe = 0.0045;
S_front = ((1.77+1.77*.55)*12/100)*11.06;
Swet = 2*S_ + S_front;
CD0 = Cfe*Swet/S_;
CD0_Tom =  0.53;
C_D = CD0 + C_Di;

%%
L_to_D_ratio_CL = figure(6);
L_D = C_L./C_D;
plot(C_L,L_D, 'LineWidth',1.5)
ylabel('$\frac{C_L}{C_{Di}}$ $[-]$','Interpreter','latex')
xlabel('$C_L$ $[-]$','Interpreter','latex')
grid minor

%%
L_to_D_ratio_AoA = figure(7);
L_D = C_L./C_D;
plot(AoA_plot,L_D, 'LineWidth',1.5)
ylabel('$\frac{C_L}{C_D}$ $[-]$','Interpreter','latex')
xlabel('$\alpha$ $[^\circ]$','Interpreter','latex')
grid minor

%%
CD_AoA = figure(8);
plot(AoA_plot,C_D, 'LineWidth',1.5)
ylabel('$C_D$ $[-]$','Interpreter','latex')
xlabel('$\alpha$ $[^\circ]$','Interpreter','latex')
grid minor
% fig2pdf(CD_AoA, 'CD_AoA', 1, 1.5, '.Figures/')

%% FLAPS 2slot
% D_Clmax = 1.6*1.1+1.6*1.1*0.8;
% def_angle = 10;
% S_flap = (11.06/2)*0.6*(1.77*1.1)*.3; % [m]
% S_ref = 15.3/2;
% D_CLmax = 0.9*D_Clmax*(S_flap/(S_ref))*cosd(def_angle);
% 
% CLmax_Flaps = CL_max + D_CLmax;

%% FLAPS 1slot
D_Clmax = cl_max*1.25+cl_max*1.25*0.85;
def_angle = 10;
S_flap = (11.06/2)*0.6*(1.77*1.25)*.4; % [m]
S_ref = 15.3/2;
D_CLmax = 0.9*D_Clmax*(S_flap/(S_ref))*cosd(def_angle);

CLmax_Flaps = CL_max + D_CLmax;