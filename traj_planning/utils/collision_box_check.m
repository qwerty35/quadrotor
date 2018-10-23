function c = collision_box_check(map, box, margin)
    if size(box,2) == 2
        box_array = box_cell2array(box);
    else
        box_array = box;
    end
        
    if size(map.blocks, 1) < 1
        c = 0;
        return;
    end
    
    if any([box_array(1:3) < map.boundary(1:3) map.boundary(4:6) < box_array(4:6)])
        c = 2; % out of boundary
        return;
    end
    
    for i = 1 : size(map.blocks, 1)
        block = map.blocks(i,1:6) + [-1 -1 -1 1 1 1] * margin;
        c = 1;
        for j = 1:3
            if box_array(j+3) <= block(j)  || block(j+3) <= box_array(j)
               c = 0;
               break;
            end
        end
        if c == 1
            break;
        end
    end
end