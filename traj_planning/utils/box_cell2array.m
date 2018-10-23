function box_array = box_cell2array(box_cell)
    center = box_cell{1};
    boxSize = box_cell{2};
    box_array = [center-boxSize center+boxSize];
end