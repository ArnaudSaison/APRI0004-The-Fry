% clear
% clc

%%
Re = 7.68e6;
naca = importdata('NACA23012.txt');
polar1 = importdata('polar_NACA_23012.txt');
polar = polar1.data;

aoa = polar(:,1);
cl = polar(:,2);
cd = polar(:,3);

%%
fig.airfoil = figure(1);
plot(naca(:,1),naca(:,2),'LineWidth',1.5)
axis equal
xlabel('$x/c$  $[-]$','Interpreter','latex')
ylabel('$y$  $[-]$','Interpreter','latex')
grid minor
% fig2pdf(fig.airfoil, 'naca_23012', 1, 1.5, '.Figures/')

%%
fig.lift_curve = figure(2);
plot(aoa,cl,'LineWidth',1.5)
hold on
plot(aoa(15),cl(15),'o','LineWidth',1.5)
xlabel('$\alpha$  $[^\circ]$','Interpreter','latex')
ylabel('$C_l$  $[-]$','Interpreter','latex')
grid minor
% fig2pdf(fig.lift_curve, 'lift_curve', 1, 1.5, '.Figures/')

%%
fig.drag_polar = figure(3);
plot(cd,cl,'LineWidth',1.5)
xlabel('$C_d$  $[-]$','Interpreter','latex')
ylabel('$C_l$  $[-]$','Interpreter','latex')
% hold on
% plot(cd(15),cl(15),'o','LineWidth',1.5)
grid minor
xlim([0 0.04])
% fig2pdf(fig.drag_polar, 'drag_polar', 1, 1.5, '.Figures/')

