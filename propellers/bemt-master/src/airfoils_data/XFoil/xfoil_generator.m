clear all; close all; clc;


airfoil = 'mh114';
reynolds_list = 100e3:100e3:2000e3;
alpha_0 = -5;
alpha_fin = 30;
alpha_delta = 0.2;

for reynolds = reynolds_list
    cmd = sprintf(['load ', airfoil,'.dat', '\n', ...
                  'pane', '\n', ...
                  'oper', '\n', ...
                  'iter 500', '\n', ...
                  'v ', num2str(reynolds), '\n', ...
                  'pacc', '\n', ...
                  airfoil, '_', num2str(reynolds),'.txt', '\n', ...
                  ' ', '\n', ...
                  'aseq ', num2str(alpha_0),' ', num2str(alpha_fin), ' ', num2str(alpha_delta), '\n', ...
                  ' ', '\n', ...
                  'quit', '\n']);
    disp(cmd)
    clipboard('copy',cmd)
	input('Press enter for next value (or end the script if nothing left)\n')
end

