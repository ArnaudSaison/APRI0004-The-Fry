drawwing
close all
nb_panels_chord = 10; % number of panels along x
nb_panels_span = 24; % number of panels along y

%%
mean_forces_Z_integral = zeros(2*nb_panels_span-1,nb_panels_chord-1);
mean_forces_Z = zeros(2*nb_panels_span,1);
for i=0:nb_panels_span*2-1
    for j=1:nb_panels_chord-1
        mean_forces_Z_integral(i+1,j) = (Wing_force_z(i*nb_panels_chord+j)+Wing_force_z(i*nb_panels_chord+j+1))/2;
    end
    mean_forces_Z(i+1) = sum(mean_forces_Z_integral(i+1,:));
end

figure
hold on
quiver3(wingcp(3:nb_panels_chord:end,1),wingcp(3:nb_panels_chord:end,2),wingcp(3:nb_panels_chord:end,3),zeros(48,1),zeros(48,1),mean_forces_Z)
view(90,0)

figure 
hold on
quiver3((1:24)',zeros(24,1),zeros(24,1),zeros(24,1),zeros(24,1),mean_forces_Z(1:24),0)
quiver3((-23:0)',zeros(24,1),zeros(24,1),zeros(24,1),zeros(24,1),flip(mean_forces_Z(25:48)),0)
view(0,0)

figure 
plot((1:24)',mean_forces_Z(1:24),'o')

%%
S = 15.3; % surface
MAC = 1.41; 
AR = 8; % aspect ratio
V = 88.0; % speed
[~,~,~,rho] = atmosisa(3.9624e+03);
c = linspace(1.77, 1.77*0.55,nb_panels_span*2);
C_L_span = mean_forces_Z./(0.5*rho*V^2*c);

figure 
% plot((1:24)',C_L_span(1:24),'o')
plot(C_L_span,'o')

%% CL-alpha wing

S = 15.3; % surface
AR = 8; % aspect ratio
V = 88.0; % speed
[~,~,~,rho_0] = atmosisa(0); % air density at sea level
alpha = [0:0.5:10,15,20,25,30,35,40,48,49,50,55]; % AoA [°]

% induced drag from VLM
% Di = [473.527581 577.387071 691.526252 815.911519 950.505525 ...
%       1095.267224 1250.151916 1415.111291 1590.093481 1775.043102 ...
%       1969.901291 2174.605742 2389.090726 2613.287109 2847.122363 ...
%       3090.520562 3343.402377 3605.685080 3877.282489 4158.105055 ...
%       4448.059814]; 

% force Z from VLM  
F_z = [24655.602757 27613.117523 30568.906858 33522.276365 36472.531475 ...
       39418.977563 42360.920073 45297.664673 48228.517402 51152.784843 ...
       54069.774292 56978.793938 59879.153029 62770.162053 65651.132904 ...
       68521.379049 71380.215668 74226.959816 77060.930476 79881.448812 ...
       82687.838218 109827.439207 134830.736863 157080.320338  ...
       176011.509807 191115.175954 201947.447562 209535.045574 ...
       209601.585247 209467.428518 205746.302775]; 

C_L = F_z/(0.5*rho_0*V^2*S);
% C_Di = Di/(0.5*rho_0*V^2*S);

% e = C_L.^2/pi/AR/C_Di; % Oswald efficiency factor
slope = (C_L(7)-C_L(4))/(alpha(7)-alpha(4));

fprintf('Lift coefficient = %.4f \n', C_L)

plot(alpha, C_L, '-o', 'LineWidth', 1.5)
xlabel('Angle of attack')
ylabel('Lift Coefficient')

%% CL-alpha HT
% S = 4.30; % surface
% AR = 4; % aspect ratio
% V = 88.0; % speed
% [~,~,~,rho_0] = atmosisa(0); % air density at sea level
% alpha_tail = 0:0.5:4; % AoA [°]
% 
% F_z_tail = [0 660.125591 1320.104342 1979.789427 2639.034038 3297.691413 ...
%             3955.614852 4612.657745 5268.673593];
% 
% C_L_tail = F_z_tail/(0.5*rho_0*V^2*S);
% 
% slope_tail = (C_L_tail(4)-C_L_tail(3))/(alpha_tail(4)-alpha_tail(3));
% 
% plot(alpha_tail, C_L_tail, '-o', 'LineWidth', 1.5)
% xlabel('Angle of attack')
% ylabel('Lift Coefficient')

%% CDi CL

S = 15.3; % surface
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

   
C_L = F_z/(0.5*rho_0*V^2*S);
C_Di = Di/(0.5*rho_0*V^2*S);
e = C_L.^2/pi/AR/C_Di;
figure
plot(C_Di,C_L)

figure
plot(-6:10,C_L,'-o')
