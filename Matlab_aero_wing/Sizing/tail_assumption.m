%% Diamond DA42

S_ref = 161.459; %trapezoidal wing surface ft^2
AR_wing = 8; %aspect Ratio
b_ref=sqrt(AR_wing*S_ref); %wingspan
Lambda=0;%wing sweep at 25% MAC
lambda=0.35;%taper ratio
c_ref=b_ref/AR_wing;%mean chord

l_HT = 20; %ft
l_VT = l_HT;

V_HT = 0.8;
V_VT = 0.07;
    
lambda_HT = 0.481;
lambda_VT = 0.511;

S_HT = V_HT * S_ref * c_ref / l_HT;
b_HT = sqrt(AR_wing * S_HT);

S_VT = V_VT * S_ref * c_ref / l_VT;
b_VT = sqrt(AR_wing * S_VT);