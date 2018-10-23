function [ desired_state ] = trajectory_box_generator(t, qn, map, path)
% TRAJECTORY_BOX_GENERATOR: Turn a Dijkstra or A* path into a box trajectory
%
% NOTE: This function would be called with variable number of input
% arguments. In init_script, it will be called with arguments
% trajectory_generator([], [], map, path) and later, in test_trajectory,
% it will be called with only t and qn as arguments, so your code should
% be able to handle that. This can be done by checking the number of
% arguments to the function using the "nargin" variable, check the
% MATLAB documentation for more information.
%
% map: The map structure returned by your load_map function
% path: This is the path returned by your planner (dijkstra function)
%
% desired_state: Contains all the information that is passed to the
% controller, as in phase 2
%
% It is suggested to use "persistent" variables to store map and path
% during the initialization call of trajectory_generator, e.g.
% persistent map0 path0
% map0 = map;
% path0 = path;
persistent map0 path0 total_time X ts Npoly;
if numel(t) == 0 | numel(qn) == 0
    map0 = map;
    path0 = path;
    
    [box, dup, total_time, ts] = generate_box(map0,path0{1},0); % speed 2로 간다고 가정하고 time segment를 나눔
    map.box = box;
    map.dup = dup;
    
    tic
    [X,Npoly,~,~] = traj_opt_Gao(path0{1},map,ts); 
    toc

    % plot box, segment points
    nblocks = size(map.blocks,1);
    map.blocks = [map.blocks; box]; % box or dup
    segments = [X(Npoly+1:Npoly+1:end,:); path0{1}(end,:)];
    plot_path(map,path0{1},segments);
    map.blocks(nblocks+1:end, :) = [];
    return
end

if nargin < 4
    map = map0;
    path = path0;
end

p = path{qn};
if t >= total_time
    pos = p(end,:);
    vel = [0;0;0];
    acc = [0;0;0];
    jerk = [0;0;0];
else
     % Npolyth order path planning
     m = find(ts<=t);
     m = m(end);
     
     pos = timevector(t-ts(m),Npoly,0)*X((Npoly+1)*(m-1)+1:(Npoly+1)*m,:);
     vel = timevector(t-ts(m),Npoly,1)*X((Npoly+1)*(m-1)+1:(Npoly+1)*m,:);
     acc = timevector(t-ts(m),Npoly,2)*X((Npoly+1)*(m-1)+1:(Npoly+1)*m,:);
     jerk = timevector(t-ts(m),Npoly,3)*X((Npoly+1)*(m-1)+1:(Npoly+1)*m,:);
end

yaw = 0;
yawdot = 0;

desired_state.pos = pos(:);
desired_state.vel = vel(:);
desired_state.acc = acc(:);
desired_state.jerk = jerk(:);
desired_state.yaw = yaw;
desired_state.yawdot = yawdot;

