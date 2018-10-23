function [vfield] = velocity_field(map, v_max)

% init distance map
dist = (~map.occ_map).*inf; % init as inf distance
dist(isnan(dist)) = 0;

% init visited map
% visited = map.occ_map;
unvisited = ones(map.nx, map.ny, map.nz); % unvisited: 1, visited: inf

% 6-connected
dijk = [1,0,0;-1,0,0; 0,1,0; 0,-1,0; 0,0,1; 0,0,-1];
% dijkstra algorithm
while (~all(unvisited(:)==inf))
    % find the min dist pos from unvisited node
    dist_unvisited = dist.*unvisited;
    [~, id] = min(dist_unvisited(:));
    % mark visited
    unvisited(id) = inf;
    % check 6-connected neighbors
    [i,j,k] = ind2sub([map.nx,map.ny,map.nz],id);
    for d = 1:size(dijk,1)
        nijk = bsxfun(@plus, int32([i,j,k]), int32(dijk(d,:)));
        % if neighbor (in the map) and (unvisited) and (unoccupied)
        if all(nijk > 0) & all(int32([map.nx, map.ny, map.nz]) >= int32(nijk)) ...
            & unvisited(nijk(1),nijk(2),nijk(3)) == 1 ...
            & map.occ_map(nijk(1),nijk(2),nijk(3)) ~= 1
            alt = dist(id) + sqrt(sum(dijk(d,:).^2));
            % Got rid of function call for optimization
            nid = (nijk(3)-1)*map.nx*map.ny + (nijk(2)-1)*map.nx + nijk(1);
            if alt < dist(nid)
                dist(nid) = alt;
            end
        end
    end
end

vfield = v_max * (tanh(dist-exp(1))+1)/2;

end
