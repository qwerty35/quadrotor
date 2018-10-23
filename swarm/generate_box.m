function [box, dup, total_time, ts] = generate_box(map, path, margin)
% This code is based on below paper
% Online Generation of Collision-Free Trajectories for Quadrotor Flight in
% Unknown Cluttered Environment
% Jing Chen, Tianbo Liu, and Shaojie Shen
%
% This code is written by
% Jungwon Park
% Seoul National University
% wjstk335@gmail.com
%
% GENERATE_BOX : generate box path that safe from collision
%
% input: map
%        path
%        margin         block collision margin
%
% output: box           matrix #path * 9 ([x0,y0,z0,x1,y1,z1,colorR,G,B ]) safe region 
%         dup           matrix #path-1 * 9 ([x0,y0,z0,x1,y1,z1,colorR,G,B]) duplicated
%                       volume between boxes
%
% Representation
% box_cell  :    {center_pts, boxsize(1,3)}, 
% box_array :    [x0,y0,z0,x1,y1,z1]
% window_array : [x0,y0,z0,x1,y1,z1]

initial_size = [map.xy_res/2, map.xy_res/2, map.z_res/2];
box = [];
dup = [];
ts = [];

% find maximal volume box by axis search
if size(map.blocks, 1) < 1
    box = [map.boundary [0 255 0]];
else
for i = 2 : size(path,1)-1
    last_box = box_cell2array({path(i,:), initial_size});
        
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
    
    box = [box; [last_box [0 255 0]]];
%     if isfield(map,'dist')
%         p = points_to_idx(map,path(i,:));
%         ts = [ts map.dist(p(1),p(2),p(3))];
%     end
end
end

% % find maximal volume box
% for i = 2 : size(path,1)-1
%     last_box = {path(i,:), initial_size};
%     new_box = {last_box{1}, last_box{2}+initial_size*2};
%     while ~collision_box_check(map, new_box, margin)  % if current box cannot grow up then go to next block
%         last_box = new_box;
%         new_box{2} = new_box{2} + initial_size * 2;
%     end
%     box = [box; [box_cell2array(last_box) [0 255 0]]];
% end
% 
% % concatenate sequential block
% i = 1;
% box_temp = zeros(1,9);
% while i < size(box,1)
%     if isSequentialBox(box(i,:), box(i+1,:))
%         box(i,1:3) = min(box(i,1:3), box(i+1,1:3));
%         box(i,4:6) = max(box(i,4:6), box(i+1,4:6));
%         box_temp = box(i+1,:);
%         box(i+1, :) = [];
%     else
%         if isSequentialBox(box_temp, box(i+1,:))
%             box(i+1,1:3) = min(box_temp(1:3), box(i+1,1:3));
%             box(i+1,4:6) = max(box_temp(4:6), box(i+1,4:6));
%         end
%         i = i+1;
%     end
% end

% delete duplicated block
i = 1;
while i < size(box,1)
    incl = isBoxinBox(box(i,:), box(i+1,:)); 
    if incl >= 0
        box(i+incl, :) = [];
%         if isfield(map,'dist')
%             ts(i+incl) = [];
%         end
    else
        i = i+1;
    end
end

% generate dup 
for i = 1 : size(box,1)-1
    dup = [dup; [max(box(i,1:3),box(i+1,1:3)) min(box(i,4:6),box(i+1,4:6)) [0 0 255]]];
end

% total_time
speed = 2.2;
path_len = sum(sqrt(sum((path(2:end, :) - path(1:end-1,:)).^2,2)));
total_time = path_len/speed;

% generate ts
% if isfield(map,'dist')
%     p = points_to_idx(map,path(end,:));
%     ts = [ts map.dist(p(1),p(2),p(3))];
%     ts = ts/ts(end);
%     ts = ts*total_time;
% else
    if size(map.blocks, 1) < 1
        path_seg_len = 1;
    else
    path_seg_len = zeros(1,size(box,1));
    path_seg_len(1) = norm(path(1,:)-(dup(1,1:3)+dup(1,4:6))/2);
    for i = 2:size(dup,1)
        path_seg_len(i) = norm((dup(i,1:3)+dup(i,4:6))/2 - (dup(i-1,1:3)+dup(i-1,4:6))/2);
        if path_seg_len(i) == 0
           path_seg_len(i) = norm((box(i,1:3)+box(i,4:6))/2 - (dup(i-1,1:3)+dup(i-1,4:6))/2);
        end
    
    path_seg_len(end) = norm(path(end,:)-(dup(end,1:3)+dup(end,4:6))/2);
    end
    end
    ts = cumsum(path_seg_len);
    ts = ts/ts(end);
    ts = [0 ts];
    ts = ts*total_time;
end

% end

function small = isBoxinBox(box_array0, box_array1)
    small = -1;
    
    if box_array0(1:3) <= box_array1(1:3) & box_array1(4:6) <= box_array0(4:6)
        small = 1;
    elseif box_array1(1:3) <= box_array0(1:3) & box_array0(4:6) <= box_array1(4:6)
        small = 0;
    end
end

function seq = isSequentialBox(box_array0, box_array1)
    seq = 0;
    
     v0 = box_array0(1:3) - box_array1(1:3);
     v1 = box_array0(4:6) - box_array1(4:6);
     
     if cross(v0, v1) == 0 & (sum(v0(:)==0) == 2 || sum(v1(:)==0) == 2)
         seq = 1;
     end
end
