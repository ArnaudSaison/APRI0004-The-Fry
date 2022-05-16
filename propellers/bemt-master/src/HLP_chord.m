function chords = HLP_chord(max_chord, nb)
% HLP_chord generates the chord for all nb elements using a spline
%

% x = [0, 0.2, 0.5, 0.9, 1];
% y = [2.2, 3,   2.6,   0.75,  0]/100;

x = ([0, 265, 779, 1038, 1296, 1558, 1819, 2077, 2337, 2856, 3117, 3376, 3503, 3635, 3670, 3706])/3706;
y = (380-[77, 45, 9, 4, 11, 26, 37, 53, 76, 130, 170, 216, 248, 305, 320, 380-1])/380*max_chord;


% x = [x, 2-flip(x(1:end-1))];
% y = [y, -flip(y(1:end-1))];

chords = interp1(x, y, linspace(0,1,nb), 'makima');

end