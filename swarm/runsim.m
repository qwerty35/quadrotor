close all;
clear all;
clc;

addpath(genpath('./'));

% Plan path 
disp('Planning ...');
map = load_map('maps/map7.txt', 0.5, 0.5, 0.5);
start = {[-0.1, 3.1, 3.0], [6.1, 3.1, 3.1]};
stop  = {[6.1, 3.1, 3.0], [0.0, 3.1, 3.1]};

%% ECBS
[path, makespan] = ECBS(map, start, stop);
plot_path(map, path);

%% Run trajectory
trajectory = test_trajectory(start, stop, map, path, makespan, true); % with visualization
