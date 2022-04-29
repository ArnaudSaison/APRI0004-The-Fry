function [Results, Elem, Param, Blades] = optiCollective(Param, Blades, target_T, a, b)
% optiCollective finds the collective pitch for a given thrust
%   Param, Blades:     	see BEMT code
%   target_T:    	[N] target thrust
%   a:             	[°] collective lower bound (must be determined beforehand through testing)
%   b:             	[°] collective upper bound (must be determined beforehand through testing)
%

% Iteration parameters
max_iter = 100; % [#] maximum nb of iterations
tol = 1e-5; % [N] tolerance on final 
i = 1;

while max_iter >= i && (b - a) > tol
    
    try
        % using middle value
        Blades.COLL_PITCHdeg = (a + b) / 2;
        [Results, ~, Param, Blades] = bemt(Param, Blades);

        % choosing new boundaries
        if target_T < Results.T
            b = Blades.COLL_PITCHdeg;

        else % greater than searched value
            a = Blades.COLL_PITCHdeg;

        end
        
    catch
        disp('error on boundaries')
        b = b - 1;
    end
    
    % iteration increment
    i = i+1;
end

warning on;
[Results, Elem, Param, Blades] = bemt(Param, Blades);
warning off;

if max_iter >= i
    disp(['Optimal collective pitch found after ', num2str(i), ' iterations']);
else
    warning('No optimal collective pitch found')
end
clear i;

end