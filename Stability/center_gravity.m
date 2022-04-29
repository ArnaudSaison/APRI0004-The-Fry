function [min_CG, max_CG] = center_gravity(empty_weight, CG_empty)
max_fuel_mass = 55;
tank_pos = 2.235+1.7/2 ;
luggage_pos = 4 ;
seat_1st_row = 1.5;
seat_2nd_row = 1.5 + 1.2 ;
luggage_weight = [0 13.61 2*13.61 3*13.61];
human_weight = 86.1825503;
payload=max_fuel_mass+4*human_weight+3*luggage_weight; 

mass_mat=zeros(2,6);
CG_mat=zeros(2,6);

for i=1:2
fuel_prct = i-1;
fuel_mass=max_fuel_mass*fuel_prct;

%% no passenger (with pilot)

CG_mat(i,1)= (empty_weight*CG_empty + human_weight*seat_1st_row + fuel_mass*tank_pos)/(empty_weight+fuel_mass + human_weight);
mass_mat(i,1)= empty_weight + fuel_mass + human_weight;

%% 1 passenger devant (with pilot)

CG_mat(i,2)= (empty_weight*CG_empty + 2*human_weight*seat_1st_row + luggage_weight(1)*luggage_pos + fuel_mass*tank_pos)/(empty_weight + 2*human_weight + luggage_weight(1) + fuel_mass);
mass_mat(i,2)= empty_weight+fuel_mass+2*human_weight;

%% 1 passenger derrière  (with pilot)

CG_mat(i,3)= (empty_weight*CG_empty + human_weight*seat_1st_row + human_weight*seat_2nd_row + luggage_weight(1)*luggage_pos + fuel_mass*tank_pos)/(empty_weight + 2*human_weight + luggage_weight(1) + fuel_mass);
mass_mat(i,3)= empty_weight+fuel_mass+2*human_weight;
%% 2 passengers 1 devant 1 derrière (with pilot)

CG_mat(i,4)= (empty_weight*CG_empty + 2*human_weight*seat_1st_row + human_weight*seat_2nd_row + luggage_weight(2)*luggage_pos + fuel_mass*tank_pos)/(empty_weight + 3*human_weight + luggage_weight(2) + fuel_mass);
mass_mat(i,4)= empty_weight+fuel_mass+3*human_weight;
%% 2 passengers derrière (with pilot)

CG_mat(i,5)= (empty_weight*CG_empty + human_weight*seat_1st_row + 2*human_weight*seat_2nd_row + luggage_weight(2)*luggage_pos + fuel_mass*tank_pos)/(empty_weight + 3*human_weight + luggage_weight(2) + fuel_mass);
mass_mat(i,5)= empty_weight+fuel_mass+3*human_weight;
%% 3 passengers (with pilot)

CG_mat(i,6)= (empty_weight*CG_empty + 2*human_weight*seat_1st_row + 2*human_weight*seat_2nd_row + luggage_weight(3)*luggage_pos + fuel_mass*tank_pos)/(empty_weight + 4*human_weight + luggage_weight(3) + fuel_mass)
mass_mat(i,6)= empty_weight+fuel_mass+4*human_weight;

end

max_CG = max(max(CG_mat));
min_CG = min(min(CG_mat));

CG_vect = [min_CG, max_CG];

end


