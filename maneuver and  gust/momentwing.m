parameters
placard
nz=9.1;
b=10.95;%m
weight_light_engine=70*nz; %N
weight_heavy_engine=120*nz; %N
linear_weight_wing=nz*12*cmean; %N/m
l=(b/2);%m
x=linspace(0,l,1000);%m
L=0.5*CL_cruise*rho(design_index)*Vcruise^2*cmean;%N/m
space=b/10;
Rv=L*l-4*weight_light_engine-weight_heavy_engine;
M=zeros(1,1000);
T=zeros(1,1000);
for i=1:1000
    M(i)= M(i)+L*x(i)^2/2;
    M(i)=M(i)-weight_heavy_engine*x(i)-linear_weight_wing*x(i)^2/2;
    T(i)=T(i)-weight_heavy_engine-linear_weight_wing*x(i)/2;
    T(i)=T(i)+L*x(i)/2;
    if x(i)>=space
     M(i)=M(i)-weight_light_engine*x(i);
     T(i)=T(i)-weight_light_engine;
    end
    if x(i)>=2*space
     M(i)=M(i)-weight_light_engine*x(i);
     T(i)=T(i)-weight_light_engine;
    end
    if x(i)>=3*space
     M(i)=M(i)-weight_light_engine*x(i);
     T(i)=T(i)-weight_light_engine;
    end
    if x(i)>=4*space
     M(i)=M(i)-weight_light_engine*x(i);
     T(i)=T(i)-weight_light_engine;
    end
    
end
min(M)
min(T)
x=x*3.28084;
figure(10)
plot(x,zeros(1,1000),'black','Linewidth',2)
hold on
plot(x,M,'red','Linewidth',2)
plot(x,T,'green','Linewidth',2)
grid on
legend('wing','bending moment','shear force','Interpreter','LaTex')
xlabel('spanwise distance [ft]','Interpreter','LaTex')