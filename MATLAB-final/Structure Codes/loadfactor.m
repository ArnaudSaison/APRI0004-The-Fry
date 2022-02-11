configuration ='take_off';
%% values 
placard
computation_altitude=design_alt;%ft
[~,alt_index]=min(abs(altitude_ft-computation_altitude));



g=9.81;%m/s
rho0 = rho(1) ;%kg/m^3 
VD = VD(alt_index)*sqrt(rho(alt_index)/rho0) %m/s
VC = VC(alt_index)*sqrt(rho(alt_index)/rho0) %m/s
Ve=linspace(0,VD,1000);%m/s
Clmaxud = 1.32 ; % Cl if the wing is upside down
%stall speed
Vs=50;% takeoff
    

%% max & min
n_max = 3.8 * ones(1,length(Ve)) ;
n_min = -1.8 * ones(1,length(Ve))  ;
n_pos=n_max;
n_neg=n_min;

%% manouever flaps up
S_w
n_m_pos = rho0*Ve.^2*S_w*Clmax_flats_up/(2*W);
n_m_neg = -rho0*Ve.^2*S_w*Clmaxud/(2*W);
slope=(VD-VC)^-1;
oo=-slope*VD;
interp=slope*Ve+oo;
n_m_neg=max(n_m_neg,interp);

[~,emplVS1] = min(abs(n_m_pos-1)) ;
Vs1 = Ve(emplVS1) ;
[~,emplVA] = min(abs(n_m_pos-n_max)) ;
VA = Ve(emplVA) ;
n_pos=min(n_m_pos,n_pos);
n_neg=max(n_m_neg,n_neg);

%% manoeuver flaps down
out=0;
j=2;
Vs0=Vs1;
Vf=VD;
if configuration == 'take_off'
    Vf=1.6*Vs1;
    [~,index_Vf] = min(abs(Ve-Vf));
    V_short=Ve(1:index_Vf);
    n_fd = rho0*V_short.^2*S_w*Clmax_flats_down/(2*W);
elseif configuration == 'approach'
    Vf=1.8*Vs1;
    [~,index_Vf] = min(abs(Ve-Vf));
    V_short=Ve(1:index_Vf);
    n_fd = rho0*V_short.^2*S_w*Clmax_flats_down/(2*W);
elseif configuration == 'landing_'
    n_fd=zeros(1,length(Ve));
    Vs0ok=0;
    while out==0
       if Ve(j)>=Vf
           out=1;
       end
       n_fd(j) = rho0*Ve(j).^2*S_w*Clmax_flats_down/(2*W); 
       if n_fd(j)>=1 && Vs0ok==0
           Vs0=Ve(j);
           Vf=1.8*Vs0;
           Vs0ok=1;
       end
       j=j+1;
    end
    V_short=Ve(1:j-1);
    n_fd=n_fd(1:j-1);
else
    print('input a correct configuration: landing_ , approach or take_off')
end
[~,emplV0] = min(abs(n_fd-1)) ;
V0 = Ve(emplV0) ;


%% gusts
load('cl23012.txt');
Cla=0.34;
c=1.369;%m
mu=(2*W)/(rho0*Cla*c*g*S_w);
F=(0.88*mu)/(5.3+mu);
%% gust 2
%Ue
UeVC=interp1([0 15000],[56 44],altitude_ft);
UeVD=interp1([0 15000],[44 22],altitude_ft);
UeVC=UeVC(design_index);
UeVD=UeVD(design_index);
%VB
VDval=1+(rho0*VD*S_w*F*Cla*UeVD)/(2*W);
VBdroite=1+(rho0*Ve*S_w*F*Cla*UeVC)/(2*W);

[~,vbind]=min(abs(VBdroite-n_m_pos));
VB=Ve(vbind);
[~,vcempl]=min(abs(Ve-VC));
gustplus=[VBdroite(1:vcempl) interp1([VC VD ],[VBdroite(vcempl) VDval],Ve(vcempl+1:end))];
gustmoins=2-gustplus;

[~,Vinit_ind]=min(abs(n_m_pos-gustmoins));
gust_envpos=min(n_m_pos,gustplus);
gust_envpos=gust_envpos(Vinit_ind:end);
gust_envneg=gustmoins(Vinit_ind:end);
nlimit=max(gustplus)
nultimate=1.5*nlimit
%% plots 
lw=2;
yspeeds=-10:12;
xspeeds=ones(1,length(yspeeds));
% maneuver
figure(5)
axis([0 100 -3 6])
grid on
hold on 
plot([flip(Ve) Ve],[flip(n_m_neg) n_m_pos],'blue','Linewidth',lw);
plot(V_short,n_fd,'green','Linewidth',lw);
plot([Ve flip(Ve)],[n_pos flip(n_neg)],'red','Linewidth',lw)

xlabel('Equivalent airspeed [ft/s]')
ylabel('Load factor [g]')

plot(xspeeds*Vs1,yspeeds,'k--')
plot(xspeeds*VC,yspeeds,'k--')
plot(xspeeds*VD,yspeeds,'k--')
plot(xspeeds*VA,yspeeds,'k--')
set(gca,'xtick',[Vs1 VA VC VD],'xticklabel',["Vs1" "VA" "VC" "VD"])
legend( 'Maneuver flaps up', 'Maneuver flaps down','Maneuver envelope')
hold off
% gust
figure(6)
plot( Ve,n_m_pos,'blue','Linewidth',lw);
grid on
hold on
plot([Ve flip(Ve)],[gustmoins flip(gustplus) ],'green','Linewidth',lw)
plot([Ve(Vinit_ind:end) flip(Ve(Vinit_ind:end))],[gust_envpos flip(gust_envneg)],'red','Linewidth',lw)
xlabel('equivalent airspeed [ft/s]')
ylabel('load factor [g]')
plot(1:100,nultimate*ones(1,100),'c--')
t1=text(VD/2,nultimate+0.5,'Ultimate load factor');
t1.Color='Cyan';
plot(xspeeds*VB,yspeeds,'k--')
plot(xspeeds*VD,yspeeds,'k--')
plot(xspeeds*VC,yspeeds,'k--')
set(gca,'xtick',[VB VC VD],'xticklabel',["VB" "VC" "VD"])
legend('Maneuver flaps down','Gust load factor','Gust envelope')
axis([0 100 -6 11])
hold off

% %% dernier plot
% envplus=max(n_pos(Vinit_ind:end),gust_envpos);
% envmoins=min(n_neg(Vinit_ind:end),gust_envneg);
% figure(7)
% plot([Ve(Vinit_ind:end) flip(Ve(Vinit_ind:end))],[envmoins flip(envplus)],'red','Linewidth',lw)
% hold on



