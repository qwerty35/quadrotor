function [path,dist] = fmm(map, start, goal)
% FMM_STAR
%
% 
% This code is modified by
% Jungwon Park
% Seoul National University
% wjstk335@gmail.com

path = [];
start_xyz = start;
goal_xyz = goal;

start = points_to_idx(map, start);
goal = points_to_idx(map, goal);

% init distance map
dist = inf(size(map.occ_map)); % init as inf distance
time_field = 1./velocity_field(map,2);

% init state
State = double(map.occ_map); % Accepted: 1, Far: 0, considered: -1
Considered = [];
Considered_cost = [];
recompute = [];

% init back tracking map
prev_cost = zeros(size(map.occ_map,1), size(map.occ_map,2), size(map.occ_map,3), 3);

% 6-connected
dijk = [1,0,0;-1,0,0; 0,1,0; 0,-1,0; 0,0,1; 0,0,-1];

%
% FMM algorithm
%
% Initialization
dist(start(1), start(2), start(3)) = 0; % set start dist as 0
State(start(1), start(2), start(3)) = 1;

for d = 1:size(dijk,1)
    nijk = start + dijk(d,:);
    if nijk > 0 & [map.nx, map.ny, map.nz] >= nijk ...
        & State(nijk(1),nijk(2),nijk(3)) == 0
        Considered = [Considered; nijk];
        Considered_cost = [Considered_cost; inf];
        State(nijk(1),nijk(2),nijk(3)) = -1;
    end
    recompute = (1:size(Considered,1))';
end

while (~isempty(Considered))
    % Eikonal
    for index = 1:size(recompute,1)
        cijk = Considered(recompute(index),:);
        n = 0;
        sum = 0;
        square_sum = 0;
        u_cand_min = inf;
        
        for d = 1:size(dijk,1)
            nijk = cijk + dijk(d,:);
            if nijk > 0 & [map.nx, map.ny, map.nz] >= nijk ...
                & State(nijk(1),nijk(2),nijk(3)) == 1
                n = n + 1;
                sum = sum + dist(nijk(1),nijk(2),nijk(3));
                square_sum = square_sum + dist(nijk(1),nijk(2),nijk(3))^2;
                if dist(nijk(1),nijk(2),nijk(3)) < u_cand_min
                    u_cand_min = dist(nijk(1),nijk(2),nijk(3));
                end
            end
        end
        temp = sum^2 - n*(square_sum-time_field(cijk(1),cijk(2),cijk(3))^2);
        if temp > 0
            Considered_cost(recompute(index)) = (sum+sqrt(temp))/n;
        else
            Considered_cost(recompute(index)) = time_field(cijk(1),cijk(2),cijk(3)) + u_cand_min;
        end
    end
    
    [~,minIndex] = min(Considered_cost);
    aijk = Considered(minIndex,:);
    State(aijk(1),aijk(2),aijk(3)) = 1;
    dist(aijk(1),aijk(2),aijk(3)) = Considered_cost(minIndex);
    
    Considered(minIndex,:) = [];
    Considered_cost(minIndex,:) = [];
    recompute = [];
    
    for d = 1:size(dijk,1)
        nijk = aijk + dijk(d,:);
        if nijk > 0 & [map.nx, map.ny, map.nz] >= nijk ...
            & State(nijk(1),nijk(2),nijk(3)) == 0
            Considered = [Considered; nijk];
            Considered_cost = [Considered_cost; inf];
            recompute = [recompute; size(Considered,1)];
            State(nijk(1),nijk(2),nijk(3)) = -1;
        end
    end
end

if dist(goal(1), goal(2), goal(3)) == inf
    path = [];
else
    cur_p = goal;
    path = [goal(1), goal(2), goal(3)];
    while any(cur_p ~= start)
        prev_cost = inf;
        for d = 1:size(dijk,1)
             pijk = cur_p + dijk(d,:);
             if all(pijk > 0) & all(int32([map.nx, map.ny, map.nz]) >= int32(pijk))
                 dist_cand = dist(pijk(1),pijk(2),pijk(3));
                 if dist_cand < prev_cost
                    prev_p = pijk;
                    prev_cost = dist_cand;
                 end
             end
        end
        cur_p = prev_p;
        path = [cur_p;path];
    end
    path = idx_to_points(map, path);
end
if size(path, 1) > 0
    path = [start_xyz; path; goal_xyz];
else
    path = zeros(0, 3);
end

end

