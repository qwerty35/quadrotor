function [path,makespan] = ECBS(map,start,stop)
    fileID = fopen('maps/output_8agents.yaml');
    textscan(fileID,'%s',2,'Delimiter','\n');
    line = textscan(fileID,'%s %f',1);
    makespan = line{2};
    textscan(fileID,'%s',5,'Delimiter','\n');
    line = textscan(fileID,'%s');
    parser = line{1};
    
    path = cell(sum(startsWith(parser, 'agent')),1);
    qi = 0;
    iter = 1;
    while iter <= size(parser,1)
        p = parser{iter};
        if startsWith(p, 'agent')
            qi = qi + 1;
            iter = iter + 1;
            continue;
        end
        
        path{qi} = [path{qi}; [str2double(parser{iter+2}) str2double(parser{iter+4}) str2double(parser{iter+6})]];
        iter = iter+9;
    end
    
    for qi = 1:size(path,1)
        path{qi} = [start{qi}; path{qi}; stop{qi}];
    end
end

