function [ desired_state ] = trajectory_swarm_generator(t, qi, map, path, makespan, params)
% TRAJECTORY_SWARM_GENERATOR: Turn an ECBS path into a box trajectory
% 
% map : struct, The map structure returned by your load_map function
% path: cell{qn}, This is the path returned by your planner (ECBS)
%
% desired_state: Contains all the information that is passed to the
% controller, as in phase 2
%
% It is suggested to use "persistent" variables to store map and path
% during the initialization call of trajectory_generator, e.g.
% persistent map0 path0
% map0 = map;
% path0 = path;

persistent map0 path0 Coef0 total_time0 ts0 Npoly;
if numel(t) == 0 | numel(qi) == 0
    speed = 2.2;
    grid_size = 1;
    time_step = grid_size/speed; 
    
    map0 = map;
    path0 = path;
    total_time0 = time_step * (makespan + 2); % +2: start, end points
    
    disp('Generate swarm box ...');
    tic
    [map.box_cell, ts_cell, ts0, rel] = generate_swarm_box(map0,path0,total_time0,time_step,0);
    toc
    
    [Coef0,Npoly,total_cost] = traj_opt_Park_swarm(path0,map.box_cell,ts_cell, ts0, rel, 0.3); 
    
    disp('Total_cost is ...');
    disp(total_cost);
    
    % plot box, segment points 
    segments = [Coef0{1}(Npoly+1:Npoly+1:end,:); path0{1}(end,:)];
    plot_path(map,path0,segments);
    return
end

p = path0{qi};
X = Coef0{qi};
ts = ts0;
total_time = total_time0;

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

