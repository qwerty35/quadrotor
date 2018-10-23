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
persistent map0 path0 total_time0 Coef0 ts0 Npoly;
if numel(t) == 0 | numel(qn) == 0
    map0 = map;
    path0 = path;
    
    for qi = 1:size(path0,2)
        [map.box{qi}, map.dup{qi}, total_time0{qi}, ts0{qi}] = generate_box(map0,path0{qi},0); % speed 2로 간다고 가정하고 time segment를 나눔

        % Select trajectory optimizer
%         tic
%         [Coef0,Npoly,ts] = optimize_ts_heu(path0{1},map,ts,@traj_opt_Gao);
%         toc
%         tic
%         [Coef0,Npoly,ts] = optimize_ts(path0{1},map,ts,@traj_opt_Gao);
%         toc
        tic
        [Coef0{qi},Npoly,total_cost,control_point] = traj_opt_Gao(path0{qi},map.box{qi},ts0{qi}); 
        toc
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%input 통일 필요
    end
    
%     tic
%     [Coef0,Npoly,total_cost] = traj_opt_nonlinear(path0,map.box,ts0); 
%     toc

    % plot box, segment points
    segments = [Coef0{qi}(Npoly+1:Npoly+1:end,:); path0{qi}(end,:)]; %% 전체 qi cover하게 수정해야
    plot_path(map,path0,segments);
    return
end

% if nargin < 4
%     map = map0;
%     path = path0;
%     Coef = Coef0;
%     ts = ts0;
%     total_time = total_time0;
% end

p = path0{qn};
X = Coef0{qn};
ts = ts0{qn};
total_time = total_time0{qn};

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

