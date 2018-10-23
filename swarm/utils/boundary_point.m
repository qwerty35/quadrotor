function [enter_time,exit_time] = boundary_point(X,Npoly,dup,ts,index)
    % search index th boundary point   
    margin = 1e-4;
    
    if index == 1
        error("Do not consider start point!");
    end
    
    % enter point
    valid_enter = 0;
    enter_time = inf;
    
    if is_inner_point(X((Npoly+1)*(index-1),:),dup(index-1,:),margin) % 저번게 같은 박스에 있을때
        enter_time = ts(index-1);
    else
        p_enter = X((Npoly+1)*(index-2)+1:(Npoly+1)*(index-1),:); 
        for k = 1:3
            for b = 0:1
                boundary = zeros(Npoly+1,1);
                boundary(end) = dup(index-1,k+3*b);
                enter_cand = roots(p_enter(:,k) - boundary);
                for j = 1:size(enter_cand,1)
                    if isreal(enter_cand(j)) && -margin < enter_cand(j) && enter_cand(j) < ts(index)-ts(index-1)+margin ...
                        && is_inner_point([polyval(p_enter(:,1),enter_cand(j)),polyval(p_enter(:,2),enter_cand(j)),polyval(p_enter(:,3),enter_cand(j))],dup(index-1,:), margin)
                        valid_enter = valid_enter + 1;
                        if enter_time > ts(index-1) + enter_cand(j)
                            enter_time = ts(index-1) + enter_cand(j);
                        end
                    end
                end
            end
        end
        if valid_enter < 1
            error("no valid point!");
        end
    end
    
    % exit point
    valid_exit = 0;
    exit_time = 0;
    if index+1 < size(ts,2) && is_inner_point(X((Npoly+1)*(index+1),:),dup(index-1,:), margin)
        exit_time = ts(index+1);
    else
        p_exit = X((Npoly+1)*(index-1)+1:(Npoly+1)*index,:);
        for k = 1:3
            for b = 0:1
                boundary = zeros(Npoly+1,1);
                boundary(end) = dup(index-1,k+3*b);
                exit_cand = roots(p_exit(:,k) - boundary);
                for j = 1:size(exit_cand,1)
                    if isreal(exit_cand(j)) && -margin < exit_cand(j) && exit_cand(j) < ts(index+1)-ts(index)+margin ...
                        && is_inner_point([polyval(p_exit(:,1),exit_cand(j)),polyval(p_exit(:,2),exit_cand(j)),polyval(p_exit(:,3),exit_cand(j))],dup(index-1,:), margin)
                        valid_exit = valid_exit + 1; 
                        if exit_time < ts(index) + exit_cand(j)                                
                            exit_time = ts(index) + exit_cand(j);
                        end
                    end
                end
            end
        end
        if valid_exit < 1
            exit_time = ts(index+1);        
        end
    end

end

