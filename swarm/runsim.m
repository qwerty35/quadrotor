close all;
clear all;
clc;

addpath(genpath('./'));
addpath('/opt/ibm/ILOG/CPLEX_Studio128/cplex/matlab/x86-64_linux')

% Plan path 
disp('Planning ...');
map = load_map('maps/map_empty.txt', 0.5, 0.5, 0.5);
stop  = {[6.1, 6.1, 6.1], [1.1, 1.1, 1.1], [6.1, 6.1, 1.1], [1.1, 1.1, 6.1], [6.1, 1.1, 6.1], [1.1, 6.1, 1.1], [1.1, 6.1, 6.1], [6.1, 1.1, 1.1]};
start = {[1.1, 1.1, 1.1], [6.1, 6.1, 6.1], [1.1, 1.1, 6.1], [6.1, 6.1, 1.1], [1.1, 6.1, 1.1], [6.1, 1.1, 6.1], [6.1, 1.1, 1.1], [1.1, 6.1, 6.1]};

%% ECBS
[path, makespan] = ECBS(map, start, stop);
% plot_ecbs_path(map, path);

%% Run trajectory
trajectory = test_trajectory(start, stop, map, path, makespan, false); % with visualization
