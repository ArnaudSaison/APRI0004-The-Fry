%% Conceptual Design of a STOL Aircraft 
% (AIAA 2022 Aircraft Deisgn competition)
%
% Code by: Adrien BOURGUIGNON, Arnaud SAISON, Augustin SANDRONT, Maxence 
% CASAGRANDE, Maxime DUMONT, Veysi ASLANCI, Robin VESTRAETE, Tom DETHIER
% 
% Academic year: 2021-2022
% University: Université de Liège - Faculté des Sciences Appliquées
% Master in Aerospace Engineering
% Course: Aerospace Design Project
% 
% All paremeters set in 'parameters.m'
%

% METHODOLOGY FOR GRAPHING THE PAYLOAD-RANGE DIAGRAM
% 1) start with max total mass: run 'parameters.m'
% 2) impose new range
% 3) find amount of fuel required: run 'performance.m', impose battery 
%    mass, and run 'batteries.m'
% 4) check if 
%       a) maximum power of ICE is exceeded
%       b) maximum volume of fuel is exceeded
%    then
%       a) reduce total mass by increment
%       b) reduce payload mass by same amount
%       c) check if payload mass is still > 0
%       d) recalculate tot mass dependent parameters
%       e) go back to (3) after saving exit code
%    elseif abs(mass left after fuel) > margin
%       a) substract mass from payload
%       b) check if payload mass is still > 0
%       c) go back to (3) after saving exit code
%    else (enough power, enough volume and mass approx tot mass)
%       a) save solution parameters
%       b) go back to (2) because new total mass should not exceed last
%          total mass, same for payload mass -> will speed up iterative
%          process


clc; clear; close all; %warning off;
debug = 0;

%% Loops for ranges and iterative process

% 1) start with max total mass: run 'parameters.m'
par = parameters();

pr.nb_ranges = 300;
pr.max_ICE_power = 175* 746; % [W]
pr.max_fuel = 75; % [l]
pr.margin = 0.5; % [kg]
pr.mass_decrement2 = 0.5; % [kg]
pr.mass_decrement = pr.mass_decrement2;
pr.min_payload = 0; % [kg]

pr.ranges = linspace(par.range, par.range * 3, pr.nb_ranges);
pr.exit_codes = zeros(1, pr.nb_ranges); % 0 = tot mass | 1 = power | 2 = volume | 3 = impossible (no payload)
pr.payload_masses = zeros(1, pr.nb_ranges);
pr.to_powers = zeros(1, pr.nb_ranges);
pr.percentage_electic = zeros(1, pr.nb_ranges);
pr.percentage_to = zeros(1, pr.nb_ranges);
pr.percentage_climb = zeros(1, pr.nb_ranges);

for range_i = 1:1:pr.nb_ranges
    % 2) impose new range
    par.range = pr.ranges(range_i);
    
    % using the slope allows for faster convergence
    if range_i > 2
        pr.mass_decrement = -(pr.payload_masses(range_i-1) - pr.payload_masses(range_i-2)) * 0.8;
        
        if pr.mass_decrement < pr.mass_decrement2
            pr.mass_decrement = pr.mass_decrement2;
        end
    end
    
    while 1
        % 3) find amount of fuel required: run 'performance.m', impose battery 
        %    mass, and run 'batteries.m'
        res = performance(par);

        par.mass_battery = 175;
        res.en.mass_battery = par.mass_battery;

        res = batteries(res, par);

        % 4) check if 
        %       a) maximum power of ICE is exceeded
        %       b) maximum volume of fuel is exceeded
        if (res.en.ICE_actual_power > pr.max_ICE_power || ...
            res.en.volume_fuel * 1000 > pr.max_fuel)
        
            % DEBUGGING
            if debug == 1; disp(['   decrement = ', num2str(pr.mass_decrement)]); end
            
            %then
            % a) reduce total mass by increment
            par.M_to = par.M_to - pr.mass_decrement;

            % b) reduce payload mass by same amount
            par.mass_payload = par.mass_payload - pr.mass_decrement;
            
            % c) check if payload mass is still > 0
            if par.mass_payload < pr.min_payload
                pr.exit_codes(range_i) = 3; % = impossible (no payload)
                
                % values restored in order for next iteration to use them
                par.M_to = par.M_to + pr.mass_decrement;
                par.mass_payload = par.mass_payload + pr.mass_decrement;
                
                % DEBUGGING
                if debug == 1; disp(['   mass left = ', num2str(res.en.mass_left_after_fuel)]); end
                
                break; % DID NOT CONVERGED STOP
                
                % all parameters are = 0 (default value is kept)
                % next range is still computed for sanity check
                % if two code 3 in a row, range increase stops
                
            end
            
            % the increased decrement can only be used for the first
            % iteration, after that it needs to be reset
            pr.mass_decrement = pr.mass_decrement2;

            % d) recalculate tot mass dependent parameters
            par.total_max_mass = par.M_to;
            par.W_to = par.M_to * par.g;

            % e) go back to (3) after saving exit code
            if res.en.ICE_actual_power > pr.max_ICE_power
                pr.exit_codes(range_i) = 1; % = not enough power
                
            elseif res.en.volume_fuel * 1000 > pr.max_fuel
                pr.exit_codes(range_i) = 2; % = not enough volume
                
            end
            
        % elseif abs(mass left after fuel) > margin
        elseif (abs(res.en.mass_left_after_fuel) > pr.margin)
            
            % a) substract mass from payload
            par.mass_payload = par.mass_payload + res.en.mass_left_after_fuel;
            
            % b) go back to (3) after saving exit code
            % default exit code = 0
        
        % else (enough power, enough volume and mass approx tot mass)
        else
            
            % a) save solution parameters
            pr.payload_masses(range_i) = par.mass_payload;
            pr.to_powers(range_i) = res.to.P_eng;
            pr.percentage_electic(range_i) = res.en.battery_tot_energy_final / res.en.tot_energy * 100;
            pr.percentage_to(range_i) = res.en.tot_to_energy / res.en.tot_energy * 100;
            pr.percentage_climb(range_i) = res.en.tot_climb_energy / res.en.tot_energy * 100;
            
            % b) go back to (2) because new total mass should not exceed last
            %    total mass, same for payload mass -> will speed up iterative
            %    process
            
            % DEBUGGING
            if debug == 1; disp(['   mass left = ', num2str(res.en.mass_left_after_fuel)]); end
            
            break; % CONVERGED STOP
            
        end
        
    end
    
    % printing exit code in command window
    disp(['Exit code = ', num2str(pr.exit_codes(range_i)), ' for payload mass = ', num2str(pr.payload_masses(range_i))]);
    disp(' ')
    
    % if two code 3 (= minimum payload reached) in a row, range increase stops
    if (range_i >= 2 && ...
        pr.exit_codes(range_i) == 3 && ...
        pr.exit_codes(range_i-1) == 3)
        
        disp(['Impossible above range of ', num2str(pr.ranges(range_i)/1.852e3), ' [nmi]']);
        break;
        
    end
    
end

% cleanup of the results
pr.i.lim_pow = find(pr.exit_codes==1);
pr.i.lim_vol = find(pr.exit_codes==2);
pr.i.lim_imp = find(pr.exit_codes==3);

pr.ranges = [0 pr.ranges(1:1:pr.i.lim_imp)];
pr.payload_masses = [pr.payload_masses(1) pr.payload_masses(1:1:pr.i.lim_imp)];
pr.to_powers = [pr.to_powers(1) pr.to_powers(1:1:pr.i.lim_imp)];
pr.percentage_electic = [pr.percentage_electic(1) pr.percentage_electic(1:1:pr.i.lim_imp(1)-1) pr.percentage_electic(pr.i.lim_imp(1)-1)];
pr.percentage_to = [pr.percentage_to(1) pr.percentage_to(1:1:pr.i.lim_imp)];
pr.percentage_climb = [pr.percentage_climb(1) pr.percentage_climb(1:1:pr.i.lim_imp)];



%% Graph

close all;

graph.size = 0.75;
graph.AR = 1.5;
graph.style1 = '-r';
graph.style2 = '-b';
graph.folder = 'Figures/payload-range/';

% Figure 1
fig.payload_range_kg = figure('Name', 'Payload-range diagram (kg)');

plot(pr.ranges/1.852e3, pr.payload_masses, graph.style1);
xlabel('range [nmi]');
grid('on');
ylabel('payload mass [kg]');
ylim([0 450]);

make_fig(fig.payload_range_kg, 'Figures/payload_range_kg', graph.size, graph.size, graph.AR);


% Figure 2
fig.payload_range_lbs = figure('Name', 'Payload-range diagram (lbs)');

colororder({'r','b'})

yyaxis left
plot(pr.ranges/1.852e3, pr.payload_masses * 2.205, graph.style1);
ylabel('payload mass [lb]');
yticks(sort(...
            round(...
                  [0, ...
                   pr.payload_masses(1)*2.205, ...
                   pr.payload_masses(1)*2.205/2], ...
                  0 ...
                 ) ...
           ) ...
      );
ylim([0 950]);

% %, ...
%                    pr.payload_masses(pr.i.lim_pow(2))*2.205, ...
%                    pr.payload_masses(pr.i.lim_vol(2))*2.205

xlabel('range [nmi]');
xticks(round([0, 300, pr.ranges(pr.i.lim_pow(2))/1.852e3, pr.ranges(pr.i.lim_vol(2))/1.852e3, pr.ranges(pr.i.lim_imp(2))/1.852e3], 0));
xtickangle(90);
grid('on');

yyaxis right
plot(pr.ranges/1.852e3, pr.percentage_electic, graph.style2); hold on;
%plot(pr.ranges(1:end-1)/1.852e3, pr.percentage_to(1:end-1) + pr.percentage_climb(1:end-1), graph.style2);
ylabel('Battery energy / total energy [%]');
ylim([0 100]);
yticks(sort([linspace(0, 100, 3), round(pr.percentage_electic(1), 1), round(pr.percentage_electic(end), 1)]));

grid('on');

make_fig(fig.payload_range_lbs, [graph.folder 'payload_range_lbs'], graph.size, graph.size, graph.AR);











