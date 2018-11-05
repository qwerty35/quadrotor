function [box_cell, ts_cell, ts_total, rel] = generate_swarm_box(map, path, total_time, time_step, margin)
% This code is written by
% Jungwon Park
% Seoul National University
% Intelligent Control Systems Laboratory
% wjstk335@gmail.com
%
% GENERATE_SWARM_BOX : generate box path that safe from collision
%
% input: map            struct, 
%        path           cell{ path_size x xyz }
%        total_time     scalar, 
%        time_step      scalar, 
%        margin         scalar, block collision margin
%
% output: box_cell      cell{ boxpath_size x box }
%         ts_cell       cell{ array: boxpath_size+1 }
%         ts_total      array, total time segment
%         rel           matrix[ rel_size x rel_info ]
%
% Representation
% box       : safe flight corridor, [x0,y0,z0,x1,y1,z1,colorR,G,B ]  
% box_array : [x0,y0,z0,x1,y1,z1]
% box_cell  : cell{ matrix : boxpath x box }
% rel_info  : [qi, qj, sector, t_trans]
% sector    : +1:+x, -1:-x, +2:+y, -2:-y, +3:+z, -3:-z

qn = size(path,1);
initial_size = [map.xy_res/2, map.xy_res/2, map.z_res/2];
box_color = [0 255 0];
box_cell = cell(qn,1);
ts_cell = cell(qn,1);
ts_total = [];
rel = [];

for qi = 1:qn
    path_max = size(path{qi},1);
    box = [];
    ts = [];
    
    % find maximal volume box by axis search
    if size(map.blocks, 1) < 1
        box = [map.boundary box_color];
    else
        for i = 2 : path_max-1
            last_box = [path{qi}(i,:)-initial_size path{qi}(i,:)+initial_size box_color];

            % axis search
            for axis = 1:6 % -x,-y,-z,+x,+y,+z
                new_box = last_box;
                while collision_box_check(map, new_box, margin) == 0
                    last_box = new_box;
                    if axis < 4
                        new_box(axis) = new_box(axis) - initial_size(axis) * 2;
                    else
                        new_box(axis) = new_box(axis) + initial_size(axis-3) * 2;
                    end
                end
            end
            box = [box; last_box];
        end
    end
    
    % delete duplicated block
    i = 1;
    while i < size(box,1)
        incl = isBoxinBox(box(i,:), box(i+1,:));
        if incl >= 0
            box(i+incl, :) = [];
        else
            i = i+1;
        end
    end
    
    % generate box time segment
    i = 2;
    j = 1;
    start = 0;
    
    while i < path_max
        if j < size(box,1)
            current = isPointinBox(path{qi}(i,:),box(j,:));
            next = isPointinBox(path{qi}(i,:),box(j+1,:));
            if start == 0 && current == 1 && next == 1 % Does path enter the duplicated region?
                start = i;
            elseif start > 0 && current == 0 && next == 1 % Does path exit the duplicated region?
                % update
                ts = [ts floor((i-1+start)/2)-1];
                j = j+1;
                start = 0;
                continue;
            end
            i = i+1;
        else
            break;
        end
    end
    
    ts_cell{qi} = [0 ts*time_step total_time];
    box_cell{qi} = box;
    ts_total = [ts_total ts];
end

% relative time segment
sector_range = [-3 -2 -1 1 2 3];

for qi = 1:qn
    for qj = qi+1:qn      
        [path_size,index] = sort([size(path{qi},1), size(path{qj},1)]);
        path_max = path_size(2);
        path_min = path_size(1);      
        sector_log = zeros(6,path_max);
        
        for iter = 1 : path_max
            % get rel_pose
            if iter <= path_min
                rel_pose = path{qj}(iter,:)-path{qi}(iter,:);
            elseif index(1) == 1
                rel_pose = path{qj}(iter,:)-path{qi}(end,:);
            else
                rel_pose = path{qj}(end,:)-path{qi}(iter,:);
            end
            
            % log sector information
            for i = 1:size(sector_range,2)
                sector = sector_range(i);
                if rel_pose(abs(sector))*sign(sector) > 0
                    if iter == 1
                        sector_log(i,iter) = 1;
                    else
                        sector_log(i,iter) = sector_log(i,iter-1)+1;
                    end
                end
            end
        end
        
        % find minimum jump sector path (heuristic greedy search)
        iter = path_max;
        [c_next,s_next] = max(sector_log(:,iter));
        rel = [[qi qj sector_range(s_next) total_time]; rel];
        iter = iter - c_next + 1;
        
        while iter > 1
            [c_curr,s_curr] = max(sector_log(:,iter));
            
            % if there is no intersection then allow to jump sector  
            % except jumping through quadrotor (e.g. +x -> -x jumping is not allowed)
            if c_curr <= 1
                iter = iter - 1;
                
                s_opp = size(sector_range,2) + 1 - s_next;
                [c_curr,s_curr] = max(sector_log(:,iter));
                                
                if c_curr <= 0
                    error("Invalid Path, missing link");
                elseif s_curr == s_opp
                    error("TODO: Debug duplicated max value...");
                end
                count = 2;
            else
                count = 1;
                while  sector_log(s_curr,iter+count) > 0
                    count = count + 1;
                end
            end
            
            rel = [[qi qj sector_range(s_curr) floor(iter-1+count/2)*time_step]; rel];
            ts_total = [ts_total floor(iter-1+count/2)];
            
            s_next = s_curr;
            c_next = c_curr;
            iter = iter - c_next + 1;
        end
        
    end
end

% generate total time segment
ts_total = unique(ts_total); % remove duplicated elements and sort
ts_total = [0 ts_total*time_step total_time];

end

function small = isBoxinBox(box_array0, box_array1)
    small = -1;
    
    if box_array0(1:3) <= box_array1(1:3) & box_array1(4:6) <= box_array0(4:6)
        small = 1;
    elseif box_array1(1:3) <= box_array0(1:3) & box_array0(4:6) <= box_array1(4:6)
        small = 0;
    end
end

function incl = isPointinBox(point, box_array)
    incl = 0;
    
    if box_array(1:3) <= point & point <= box_array(4:6)
        incl = 1;
    end
end

