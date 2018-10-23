close all;
clear all;
clc;

addpath(genpath('./'));

% % Plan path 1 (Hole)
% disp('Planning ...');
% map = load_map('maps/map1.txt', 0.5, 0.5, 0.25);
% start = {[0.0  -4.9 0.2]};
% stop  = {[6.0  18.0-1 5.0]};
%  stop  = {[6.0  18.0-6 3.0]};

% % Plan path 2 (simple X,Z zigzag)
% disp('Planning ...');
% map = load_map('maps/map2.txt', 0.5, 0.5, 0.25);
% start = {[0.0, 2, 5]};
% stop  = {[10, 2, 5]};

%% Plan path 2.5 (simple X,Z zigzag)
%disp('Planning ...');
%map = load_map('maps/map2.txt', 0.5, 0.5, 0.25);
%start = {[0.0, 2, 5]};
%stop  = {[10, 2, 3]};

% %% Plan path 3 (X,Z zigzag)
% disp('Planning ...');
% map = load_map('maps/map3.txt', 0.5, 0.5, 0.25);
% start = {[0.0, 2, 5.0]};
% stop  = {[20, 2, 5]};

% Plan path 4 (narrow X,Z zigzag)
disp('Planning ...');
map = load_map('maps/map4.txt', 0.5, 0.5, 0.1);
start = {[0.0, 2, 5.0]};
stop  = {[20, 4, 5]};

% % Plan path 5 (X,Y zigzag)
% disp('Planning ...');
% map = load_map('maps/map5.txt', 0.5, 0.5, 0.1);
% start = {[0.0, 5.0, 2]};
% stop  = {[20, 5, 2]};

% % Plan path 6 (One box)
% disp('Planning ...');
% map = load_map('maps/map_one_block.txt', 0.5, 0.5, 0.1);
% start = {[0.0, 0.0, 0.0], [10, 10, 3]};
% stop  = {[10, 10, 3], [0.0, 0.0, 0.0]};

% % Plan path 7 (long X,Z zigzag)
% disp('Planning ...');
% map = load_map('maps/map6.txt', 0.5, 0.5, 0.1);
% start = {[0, 1, 3]};
% stop  = {[39, 2, 5]};

% % Plan path 8 (random tree)
% disp('Planning ...');
% map = load_map('maps/map_tree.txt', 0.5, 0.5, 0.1);
% start = {[0.0  -4.9 0.2]};
% stop  = {[6.0  18.5 5.0]};

% % Plan path 9 (corridor)
% disp('Planning ...');
% map = load_map('maps/map_corridor.txt', 0.5, 0.5, 0.5);
% % start = {[0.0, 0.0, 0.0], [10, 10, 3]};
% % stop  = {[10, 10, 3], [0.0, 0.0, 0.0]};
% start = {[0.0, -3.0, 0.0], [13, 0.0, 0.0]};
% stop  = {[10, 10, 3], [0.0, 10.0, 3.0]};


nquad = length(start);
for qn = 1:nquad
   tic
   path{qn} = dijkstra(map, start{qn}, stop{qn}, true);
%    [path{qn},dist] = fmm(map, start{qn}, stop{qn}); map.dist = dist;
   toc
end
plot_path(map, path);

%% Additional init script
% init_script;

%% Run trajectory
trajectory = test_trajectory(start, stop, map, path, true); % with visualization
