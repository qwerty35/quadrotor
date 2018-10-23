function [collide, delta, normal, r_col] = collision_block_check(s, map, params)
% COLLISION_BLOCK_CHECK 
% simple collision check function
% CAUTION : it does NOT detect edge collision! plane collision only
%
% input:    s       17 x 1, state vector = [x, y, z, xd, yd, zd, qw, qx, qy, qz, p, q, r, w1, w2, w3, w4]
%           map     struct, map information
%
% output:   collide number of collided block
%           delta   dynamic array, local penetration at colide block
%           normal  dynamic cell, normal vector of collision
%           r_col   dynamic cell, vector that com -> collision point

collide = 0; 
delta = [];
normal = {};
r_col = {};

p = s(1:3)';        % current quad center of mass
q = s(7:10)';  
ez = quatrotate(q, [0 0 1]); % z basis vector of the body coordinate system

blocks = map.blocks;
n = {[-1 0 0] [0 -1 0] [1 0 0] [0 1 0]}; % block normal vector


% Get 4 colision point candidates(x -x y -y) pc
% Disc cage 
if(strcmp(params.cage_type, 'disc'))
    % Get 4 colision point candidates(x -x y -y) pc
    pc = zeros(4,3);
    r = zeros(4,3);
    for i = 1:4
        rc = 0;
        v = n{i} - dot(ez,n{i})*ez;
        if norm(v) ~= 0
            rc = -v/norm(v)*params.R; 
        end
        pc(i, :) = p+rc;
        r(i, :) = rc;
    end
% Sphere cage
elseif(strcmp(params.cage_type, 'sphere'))
    pc = zeros(4,3);
    r = zeros(4,3);
    for i = 1:4
        rc = -n{i}*params.R;
        pc(i, :) = p+rc;
        r(i, :) = rc;
    end
else
    error('cage_type error');
end

% Collision check
for i = 1 : size(blocks, 1)
    inc = -1;
    for j = 1 : 4
        if(blocks(i,1) < pc(j,1) && pc(j,1) < blocks(i,4) && ...
           blocks(i,2) < pc(j,2) && pc(j,2) < blocks(i,5) && ...
           blocks(i,3) < pc(j,3) && pc(j,3) < blocks(i,6))

           if(j<3)
               temp = -dot(pc(j,:),n{j})-blocks(i,j);
           else
               temp = -dot(pc(j,:),n{j})+blocks(i,j+1);
           end
           
           delta = [delta temp];
           normal = [normal n(j)];
           r_col = [r_col r(j,:)];
           inc = inc + 1;
         end
    end

    if(inc > 0)
        error('quad in the block!');
    elseif(inc == 0)
        collide = collide + 1;
    end
end
    


