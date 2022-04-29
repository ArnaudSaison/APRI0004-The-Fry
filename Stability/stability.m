function [ Kn , hn, xn] = stability(CG)  

 width_fus = 1.5 ; %largeur fuselage 
 length_fus = 9.455 ; %longueur du fuselage 
 
 taper_w = 0.55 ; %taper wing
 taper_t = 0.55 ; %taper horizontal tail
 root_chord_w = 1.77 ;
 root_chord_t = 1.02 ;
 MAC_w = root_chord_w * 2/3 * (( 1 + taper_w + taper_w^2 )/( 1 + taper_w ))  %mean aerodynamical chord of wing 
 MAC_t = root_chord_t * 2/3 * (( 1 + taper_t + taper_t^2 )/( 1 + taper_t ))  %mean aerodynamical chord of tail 
 
 S = 17.955  %wing surface 
 b = 12.3 ; %span wing 
 AR = b^2/S ; %aspect ratio 
 a = 0.0812*(180/pi) ; %cl slope of wing wrt aoa 
 a1 = 0.0739*(180/pi) ; %cl slope of tail 
 St = 3; %horizontal tail surface

 mb_2 = 3.9 ; %paramètre compliqué à expliquer voir slide 45 design course 
 LE_tail = 8.63 ; %x position of leading edge of horizontal tail 
 LE_wing = 2.235 ; %x position of leading edge of wing 
 M_fus = (LE_wing + MAC_w/4)/length_fus; 
 propeller_position = LE_wing ;

M=[0.1 0.2 0.3 0.4 0.5 0.6 0.7];
k=[0.115 0.172 0.344 0.487 0.688 0.888 1.146];
k_fus = interpn(M,k,M_fus);
    hc = CG - LE_wing ;
    h = hc / MAC_w ; 
    h0c = MAC_w / 4 ;
    h0 = 0.25 ;
    lT = LE_tail - CG + 0.25*MAC_t;
    lt = LE_tail + 0.25*MAC_t - (LE_wing + 0.25*MAC_w);

    m = 2*mb_2/b ;

    d_epsilon_d_alpha = 1.75*a/((pi*AR*(2*taper_w*lt/b)^0.25)*(1+abs(m)));

    d_Cmfus_d_CLw=(k_fus*width_fus^2*length_fus)/(S*MAC_w*a);

    d_CLT_d_CLw = (a1/a)*(1-d_epsilon_d_alpha);


    hn = h0 + d_CLT_d_CLw * (lT/MAC_w)*(St/S) - d_Cmfus_d_CLw;
    %stick free correction
   
    hn = hn*0.965 ;
    %power on correction
    if propeller_position < CG
        hn = hn*(1-0.02*(CG-propeller_position)/MAC_w) ;
    end
    xn = hn*MAC_w + LE_wing ;
    Kn = (hn - h);
    
 
 end