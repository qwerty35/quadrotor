function plot_ecbs_path(map,path)
    
    color = winter(size(path,1));
    set(gca, 'XLim', map.boundary([1,4]), 'YLim', map.boundary([2,5]), 'ZLim', map.boundary([3,6]))
    view(43,24);
    hold on;
    grid on;
    
    maxlength = 0;
    for qi = 1:size(path,1)
        if size(path{qi},1) > maxlength
            maxlength = size(path{qi},1);
        end
    end
    
    dt = 0.1;
    
    for i = 1/dt:maxlength/dt
        for qi = 1:size(path,1)
            
            i2 = ceil(i*dt);
            i1 = floor(i*dt);
            if i2 > size(path{qi},1)
                i2 = size(path{qi},1);
            elseif i1 > size(path{qi},1)
                i1 = size(path{qi},1);
            end
            
            if i*dt == floor(i*dt)
                a = 0;
                b = 1;
            else
                a = i*dt-floor(i*dt);
                b = ceil(i*dt)-i*dt;
            end
            
            
            h(qi) = scatter3(path{qi}(i2,1)*a + path{qi}(i1,1)*b, ...
                             path{qi}(i2,2)*a + path{qi}(i1,2)*b, ...
                             path{qi}(i2,3)*a + path{qi}(i1,3)*b, ...
                             'filled','MarkerFaceColor', color(qi,:));
        end
        drawnow
        pause(dt/3);
        delete(h);
    end
end

