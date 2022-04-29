function pitches = HLP_pitch(a, b, c, nb)
% pitches
% 
% 

x = [0, 0.5, 1];
y = [a, b, c];

pitches = spline(x, y, linspace(0,1,nb));

end