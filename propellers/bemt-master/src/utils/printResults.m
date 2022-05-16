% Print results in console
if Param.Simul.PRINT_CONSOLE == 1
    disp('--------------------------------------')
    disp('BEMT code results')
    disp('--------------------------------------')
    fprintf('Total Thrust: %0.2f N\n', Results.T)
    fprintf('Total Torque: %0.2f Nm\n', Results.Q)
    fprintf('Total Power: %0.2f kW (%0.1f hp)\n', Results.P/1000, Results.P/1000/0.746)
    disp(' ');
    fprintf('Advance ratio: %0.2f \n', Results.adv_ratio)
    fprintf('Propulsive efficiency: %0.2f %%\n', Results.eff_prop*100)
    fprintf('Rotor CT: %0.6f\n', Results.CT)
    fprintf('Rotor CQ: %0.6f\n', Results.CQ)
    disp('--------------------------------------')
    disp(' ')
end

% Plot the results
if Param.Simul.SHOW_GRAPHS == 1
    plotResultsCustom(Elem);
end
if Param.Simul.SHOW_3DVIEW == 1
    plotRotor(Blades, Elem);
end

% Saves the results
if Param.Simul.SAVE_RESULTS == 1
    % Create dir if required then save all solver outputs
    if ~exist(Param.Simul.SAVE_PATH, 'dir')
        mkdir(Param.Simul.SAVE_PATH);
    end
    save([Param.Simul.SAVE_PATH, Param.Simul.SAVE_FILENAME, '.mat'], 'Results', 'Elem', 'Param', 'Blades');
end
