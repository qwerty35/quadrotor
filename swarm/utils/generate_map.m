function [output] = generate_map(n, min_box_size, boundary)

% input n            : number of trees
%       min_box_size : minimum box size
%       boundary     : array [x_min y_min z_min x_max y_max z_max]

fileID = fopen('maps/map_tree.txt','w');
fprintf(fileID,'# random generated tree map\n\n');
fprintf(fileID,'boundary ');
fprintf(fileID,'%.2f ', boundary);
fprintf(fileID,'\n\n');

for i = 1:n
    size = min_box_size + rand;
    x = boundary(1) + rand*(boundary(4)-boundary(1));
    y = boundary(2) + rand*(boundary(5)-boundary(2));
    
    fprintf(fileID,'block  ');
    fprintf(fileID,'%.1f  ', [x-size/2 y-size/2 boundary(3) x+size/2 y+size/2 boundary(6) 255 0]); 
    fprintf(fileID,'%.1f', 0); % 255 0 0 이어야하는데 뒷부분 스페이스 때문에 0이 더 붙어버림.
    fprintf(fileID,'\n');
end

output = 1;


end

