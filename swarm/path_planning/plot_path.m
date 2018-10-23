function plot_path2(map, path, segments)
% PLOT_PATH Visualize a path through an environment
%   PLOT_PATH(map, path) creates a figure showing a path through the
%   environment.  path is an N-by-3 matrix where each row corresponds to the
%   (x, y, z) coordinates of one point along the path.

figure(1);
% path = points_to_idx(map, path);


% [x,y,z] = meshgrid(1:size(map.map3d,1),1:size(map.map3d,2),1:size(map.map3d,3));
% x = x(:); y = y(:); z = z(:);
% 
% xyz = [];
% for i = 1:numel(x)
%     if map.map3d(x(i),y(i),z(i),1)~= 255 | map.map3d(x(i),y(i),z(i),2)~= 255 | map.map3d(x(i),y(i),z(i),2)~= 255
%         xyz = [xyz; x(i),y(i),z(i)];
%     end
% end
% rgb = zeros(size(xyz,1),3);
% for i = 1:size(xyz,1)
%     rgb(i,:) = map.map3d(xyz(i,1),xyz(i,2),xyz(i,3),:);
% end
% xyz = idx_to_points(map, xyz);
% if size(xyz,1) > 0
% pcshow(xyz, rgb, 'MarkerSize', 20);
% end
% hold on;
% if size(path,1) > 0
% pcshow(path, [0,0,0],'MarkerSize', 20);
% end
% hold off;


hold on;

% additional block to show box or dup
nblocks = size(map.blocks,1);
if isfield(map,'box')
    addblock = map.box; % box or dup
    for qi = 1:size(addblock,2)
        map.blocks = [map.blocks; addblock{qi}]; 
    end
end

for i = 1:size(map.blocks,1)
    block = map.blocks(i, :);
    
    x = [ones(4,1) * block(1); ones(4,1) * block(4)];
    y = [ones(2,1) * block(5); ones(2,1) * block(2); ones(2,1) * block(5); ones(2,1) * block(2)];
    z = [block(3);block(6);block(3);block(6);block(3);block(6);block(3);block(6)];


    vert = [x, y, z];
    fac = [1 2 6 5;2 3 7 6;3 4 8 7;4 1 5 8;1 2 3 4;5 6 7 8];
    c = block(7:9)/255;
    if c == [0 1 0] | c == [0 0 1]% generated box
        patch('Vertices',vert,'Faces',fac,...
              'FaceVertexCData',hsv(6),'FaceColor',c,'FaceAlpha',0.1);
    else
        patch('Vertices',vert,'Faces',fac,...
              'FaceVertexCData',hsv(6),'FaceColor',c);
    end
    
    x = [ones(4,1) * block(1); ones(4,1) * block(4)];
    y = [block(2);block(5);block(2);block(5);block(2);block(5);block(2);block(5)];
    z = [ones(2,1) * block(3); ones(2,1) * block(6); ones(2,1) * block(3); ones(2,1) * block(6)];

    vert = [x, y, z];
    fac = [1 2 6 5;2 3 7 6;3 4 8 7;4 1 5 8;1 2 3 4;5 6 7 8];
    c = block(7:9)/255;
    if c == [0 1 0] | c == [0 0 1]% generated box
        patch('Vertices',vert,'Faces',fac,...
              'FaceVertexCData',hsv(6),'FaceColor',c,'FaceAlpha',0.1);
    else
        patch('Vertices',vert,'Faces',fac,...
              'FaceVertexCData',hsv(6),'FaceColor',c);
    end
end
map.blocks(nblocks+1:end, :) = [];

if size(path,1) > 0
    for qi = 1:size(path,2)
        pcshow(path{qi}, [0,0,0],'MarkerSize', 0.1);
    end
end
if exist('segments','var')
    pcshow(segments, [0,0,1],'MarkerSize', 20);
end
% axis([map.boundary(1)-1, map.boundary(4)-1, map.boundary(2)-1,map.boundary(5)+1,map.boundary(3)+1,map.boundary(6)+1])
hold off;
view(0, 90); % XY
% view(0,0); % XZ
end