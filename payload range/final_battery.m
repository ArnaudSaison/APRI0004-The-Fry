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
% All outputs defined in 'ouput.m'
% 

clear; close all; clc;

%% 
%==========================================================================
% Initialization of the parameters (set by user in parameters.m)
%==========================================================================
par = parameters();

%==========================================================================
% Calculating the different phases of the flight
%==========================================================================
res = performance(par);

%==========================================================================
% Batteries based on chosen design parameters
%==========================================================================
% PARAMETER (only one that should be changed here)
par.mass_battery = 175;
res.en.mass_battery = par.mass_battery;

res = batteries(res, par);

%==========================================================================
% Printing the output
%==========================================================================
output(par, res);