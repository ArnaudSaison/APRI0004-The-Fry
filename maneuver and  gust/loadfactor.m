configuration ='take_off';
%% values 
placard
computation_altitude=design_alt;%ft
[~,alt_index]=min(abs(altitude_ft-computation_altitude));



g=9.81;%m/s
rho0 = rho(1) ;%kg/m^3 
VD = VD(alt_index)*sqrt(rho(alt_index)/rho0); %m/s
VC = VC(alt_index)*sqrt(rho(alt_index)/rho0); %m/s 
Ve=linspace(0,VD,1000);%m/s
Clmaxud = 1.32*0.9 ; % Cl if the wing is upside down
%stall speed
Vs=50;% takeoff
Ve_ft=Ve*3.28084;    

%% max & min
n_max = 3.8 * ones(1,length(Ve)) ;
n_min = -1.8 * ones(1,length(Ve))  ;
n_pos=n_max;
n_neg=n_min;
W=W*9.81;
%% manouever flaps up
%W=9.81*W;
n_m_pos = rho0*Ve.^2*S_w*Clmax_flats_up/(2*W);
n_m_neg = -rho0*Ve.^2*S_w*Clmaxud/(2*W);
slope=(VD-VC)^-1;
oo=-slope*VD;
interp=slope*Ve+oo;
n_m_neg=max(n_m_neg,interp);

[~,VS1_ind]=min(abs(n_m_pos-1));
Vs1=Ve_ft(VS1_ind);

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
load('polar.txt');
c=1.369;%m
%Ue
UeVC=interp1([0 15000],[56 44],altitude_ft);
UeVD=interp1([0 15000],[28 22],altitude_ft);
UeVC=UeVC(design_index);
UeVD=UeVD(design_index);

W_lbf=W*0.224809;
S_w_ft=S_w*10.7639;
rho_US=rho_star*0.00194032;%slugs/ft^3
c_ft=c*3.28084;
VC_ft=VC*3.28084;
VD_ft=VD*3.28084;

newalpha=AOA+UeVC/VC_ft;
%spline(polar(:,1),polar(:,2),newalpha)
CLa=2*pi/((1+2*pi/(pi*AR)));

mu_g=2*(W_lbf/S_w_ft)/(rho_US*c_ft*32.174*CLa);
K=0.88*mu_g/(5.3+mu_g);

n_VC_plus=1+K*UeVC*VC_ft*CLa/(498*(W_lbf/S_w_ft));
n_VC_moins=1-K*UeVC*VC_ft*CLa/(498*(W_lbf/S_w_ft));
n_VD_plus=1+K*UeVD*VD_ft*CLa/(498*(W_lbf/S_w_ft));
n_VD_moins=1-K*UeVD*VD_ft*CLa/(498*(W_lbf/S_w_ft));

gustplus=interp1([0 VC_ft VD_ft],[1 n_VC_plus n_VD_plus],Ve_ft) ;
gustmoins=interp1([0 VC_ft VD_ft],[1 n_VC_moins n_VD_moins],Ve_ft) ;

%% full envelope
fullplus_1=min(n_m_pos,gustplus);
fullplus=max(fullplus_1,n_pos);

%% speeds & n
[~,VA_emp]=min(abs(n_m_pos-n_max));
VA_ft=Ve_ft(VA_emp);
n_VA=fullplus(VA_emp);
nnVA=gustmoins(VA_emp);

[~,VB_emp]=min(abs(n_m_pos-gustplus));
VB_ft=Ve_ft(VB_emp);
n_VB=fullplus(VB_emp);
nnVB=gustmoins(VB_emp);

n_VC=n_VC_plus;
% plots 
[~,ind_inter]=min(abs(fullplus-gustmoins));
Ve_plot=Ve_ft(ind_inter:end);

nn=-3:6;
yy=ones(length(nn),1);
fvn=figure('Name','Vn');
axis([0 115*3.28084 -4 8])
grid on
hold on 
plot(V_short*3.28084,n_fd,'green','Linewidth',lw);
xlabel('Equivalent airspeed [ft/s]')
ylabel('Load factor [-]')
plot(Ve_ft,n_m_pos)

plot([Ve_plot flip(Ve_plot)],[fullplus(ind_inter:end) flip(gustmoins(ind_inter:end))],'red','Linewidth',lw)
plot(Ve_ft,interp1([0 VC_ft],[1 n_VC_plus],Ve_ft) ,'--','Color',[0.5 0.5 0.5],'Linewidth',1)
plot(Ve_ft,interp1([0 VC_ft],[1 n_VC_moins],Ve_ft) ,'--','Color',[0.5 0.5 0.5],'Linewidth',1)
plot(Ve_ft,interp1([0 VD_ft],[1 n_VD_plus],Ve_ft) ,'--','Color',[0.5 0.5 0.5],'Linewidth',1)
plot(Ve_ft,interp1([0 VD_ft],[1 n_VD_moins],Ve_ft) ,'--','Color',[0.5 0.5 0.5],'Linewidth',1)
%plot(Ve_ft,gustmoins)
plot(yy*VA_ft,nn ,'--k','Linewidth',1)
plot(yy*VB_ft,nn ,'--k','Linewidth',1)
plot(yy*VC_ft,nn ,'--k','Linewidth',1)
plot(yy*VD_ft,nn ,'--k','Linewidth',1)
plot(yy*Vs1,nn ,'--k','Linewidth',1)
xticks([0 50 100 150 200 250 300 350])
tva=text(VA_ft-20,1,'V_A');
tvb=text(VB_ft-15,1,'V_B');
tvc=text(VC_ft+10,1,'V_C');
tvd=text(VD_ft-20,1,'V_D');
tvd=text(Vs1+10,1,'V_{s1}');

plot([VA_ft VA_ft],[n_VA nnVA],'o','Color',[0 0.8 0])
plot([VB_ft VB_ft],[n_VB nnVB],'o','Color',[0 0.8 0])
plot(VC_ft,n_VC ,'o','Color',[0 0.8 0])
plot([VD_ft VD_ft],[3.8 gustmoins(end)],'o','Color',[0 0.8 0])
tvd=text(VA_ft-12,n_VA,'A','Color',[0 0.8 0]);
tvd=text(VB_ft-12,n_VB+0.1,'B','Color',[0 0.8 0]);
tvd=text(VC_ft+5,n_VC+0.1,'C','Color',[0 0.8 0]);
tvd=text(VD_ft+5,3.8,'D1','Color',[0 0.8 0]);
tvd=text(VD_ft+5,gustmoins(end),'D2','Color',[0 0.8 0]);
tvd=text(VB_ft-12,nnVB,'E','Color',[0 0.8 0]);
tvd=text(VA_ft-12,nnVA,'F','Color',[0 0.8 0]);
fig2pdf(fvn, 'Vn', 1, 2)
hold off

