%% esthetics
lw=1;
%% parameters
parameters
g= 9.80665; %m/s
L=W*g; % m*kg/sec^2 (newton)

R=8.31;%gaz constant j/(K*mol)
MM=0.0289652; %kg/mol
%% temperature evolution
steps=15000;
T0=288.15;% K
l=0.0065; %K/m
altitude=linspace(0,15000,steps);%m
T=zeros(1,steps);
[~,altstrat]=min(abs(altitude-11000));
for i=1:altstrat
    T(i)=T0-altitude(i)*l;
end
T(find(T==0))=T(altstrat);

altitude_ft=altitude*3.28084; % conversion to feet


%% Pressure evolution
p=zeros(1,steps);
p_0=101325; %Pa
for i=1:steps
    p(i)=p_0*(1-(l*altitude(i)/T0))^(g*MM/(R*l));
end
p_US=p*0.0208854342; %pound-force/ft^2



%% density evolution
rho=zeros(1,steps);
for i=1:steps
 rho(i)=p(i)*MM/(R*T(i));
end
rho_USA=rho*0.0624279606; %lb/ft^3

%% desing altitude
rho_star=2*L/(S_w*Vcruise^(2)*CL_cruise);% density for the load factor n to be 1

[~,emplacement]=min(abs(rho-(ones(1,steps)*rho_star)));
design_alt=altitude_ft(emplacement)



%% plots
turbulencelow=12000*ones(1,length(T));
turbulencehigh=20000*ones(1,length(T));
stratosphere=ones(1,length(T))*36089.2;
design=ones(1,length(T))*design_alt;
[~,design_index]=min(abs(altitude_ft-design(1)));
cieling=ones(1,length(T))*14000;

% figure(1)
% plot(T,altitude_ft,'LineWidth',lw);
% hold on
% x=linspace(210,290,length(T));
% plot(x,turbulencelow,'--r','LineWidth',lw)
% plot(x,turbulencehigh,'--r','LineWidth',lw)
% plot(x,stratosphere,'--k','LineWidth',lw)
% plot(x,design,'--g','LineWidth',lw)
% plot(x,cieling,'--m','LineWidth',lw)
% t1=text(270,turbulencelow(1)+4000,'Turbulence zone');
% t2=text(270,design(1)+2000,'Design altitude');
% t3=text(270,stratosphere(1)+2000,'Stratosphere limit');
% t4=text(270,cieling(1)+2000,'Ceiling altitude');
% t1(1).Color='red';
% t2(1).Color='green';
% t3(1).Color='black';
% t4(1).Color='magenta';
% ylabel('altitude (ft)')
% xlabel('Temperature(K)')
% hold off
% 
% figure(2)
% plot(p_US,altitude_ft,'LineWidth',lw);
% hold on
% x=linspace(200,2200,length(T));
% plot(x,turbulencelow,'--r','LineWidth',lw)
% plot(x,turbulencehigh,'--r','LineWidth',lw)
% plot(x,stratosphere,'--k','LineWidth',lw)
% plot(x,design,'--g','LineWidth',lw)
% plot(x,cieling,'--m','LineWidth',lw)
% xpl=1600;
% t1=text(xpl,turbulencelow(1)+4000,'Turbulence zone');
% t2=text(xpl,design(1)+2000,'Design altitude');
% t3=text(xpl,stratosphere(1)+2000,'Stratosphere limit');
% t4=text(xpl,cieling(1)+2000,'Ceiling altitude');
% t1(1).Color='red';
% t2(1).Color='green';
% t3(1).Color='black';
% t4(1).Color='magenta';
% ylabel('altitude (ft)')
% xlabel('Pressure (pound-force/feet^2)')
% hold off
% 
% figure(3)
% plot(rho_USA,altitude_ft,'LineWidth',lw);
% hold on
% x=linspace(0.01,0.08,length(T));
% plot(x,turbulencelow,'--r','LineWidth',lw)
% plot(x,turbulencehigh,'--r','LineWidth',lw)
% plot(x,stratosphere,'--k','LineWidth',lw)
% plot(x,design,'--g','LineWidth',lw)
% plot(x,cieling,'--m','LineWidth',lw)
% xpl=0.060;
% t1=text(xpl,turbulencelow(1)+4000,'Turbulence zone');
% t2=text(xpl,design(1)+2000,'Design altitude');
% t3=text(xpl,stratosphere(1)+2000,'Stratotsphere limit');
% t4=text(xpl,cieling(1)+2000,'Ceiling altitude');
% t1(1).Color='red';
% t2(1).Color='green';
% t3(1).Color='black';
% t4(1).Color='magenta';
% ylabel('altitude (ft)')
% xlabel('air density (lb/ft^3)')
% hold off
% 
% 


%% desing cruise mach
sound_speed=331.3+0.606*T; %m/s
rho_star=2*L/(S_w*VC_max^(2)*CL_cruise);% density for the load factor n to be 1
[~,emplacementVC]=min(abs(rho-(ones(1,steps)*rho_star)));
M_C=VC_max/sound_speed(emplacementVC);
ConstantMC=M_C*sound_speed;
%% airspeed with height
VC=M_C*sound_speed;%m/s
dyn_pressure=rho(design_index)*VC(design_index)^2*(1/2);% has to be constant to have constant drag on descent
for i=flip(1:design_index)
 VC(i)=sqrt(dyn_pressure*2./rho(flip(i)));%m/s
end
%% maximum velocity
MD=1.25*M_C;
max_speed_MD=MD*sound_speed;
VD=max_speed_MD;
for i=flip(1:design_index)
 VD(i)=min([1.25*VC(i) max_speed_MD(i)]);%m/s
end
%% conversion to imperial

VD_f=VD*3.28084;%feet/s
VC_f=VC*3.28084;%feet/s
ConstantMC=ConstantMC*3.280084;
%% placard diagram
fplac=figure('Name','placard');
axis([200 450 0 20000])
hold on
grid on 
plot(VC_f,altitude_ft,'color',[0 0.5 0.5],'LineWidth',lw);
plot(VD_f,altitude_ft,'r','LineWidth',lw);
x=linspace(150,450,length(T));
plot(ConstantMC,altitude_ft,'--b','LineWidth',lw);
plot(x,stratosphere,'--','color',[0.5 0.5 0.5],'LineWidth',lw)
plot(x,design,'--k','LineWidth',lw)
plot(x,cieling,'--','color',[0.5 0.5 0.5],'LineWidth',lw)
xpl=205;
t2=text(xpl,design(1)-2000,'Design altitude');
t3=text(xpl,stratosphere(1)+2000,'Stratotspheric limit');
t4=text(xpl,cieling(1)+2000,'Service ceiling');
t2(1).Color='black';
t3(1).Color=[0.5 0.5 0.5];
t4(1).Color=[0.5 0.5 0.5];
t5=text(330,18000,'VC');
t6=text(420,18000,'VD');
t5(1).Color=[0 0.5 0.5];
t6(1).Color='red';
t7=text(333,10000,'Constant Mc');
t7(1).Color='Blue';
ylabel('altitude (ft)','Interpreter','LaTex')
xlabel('true airspeed (feet/s)','Interpreter','LaTex')
fig2pdf(fplac, 'placard', 1, 1.5)
hold off