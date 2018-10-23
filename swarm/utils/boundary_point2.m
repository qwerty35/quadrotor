function [result] = boundary_point2(X,Npoly,dup,ts,index)
    % search index th boundary point   
    margin = 1e-4;
    r = (Npoly+1)/2;
    
    if index <= 1 || index >= size(ts,2)
        error("Do not consider start or end point!");
    end
    
    boundary_time_xyz = zeros(3,2);
    
    % Build C
    C = zeros(Npoly+1);
    for i = 1:r
        C(i,:) = timevector(0,Npoly,r-i);
        C(i+r,:) = timevector(ts(index+1)-ts(index-1), Npoly, r-i);
    end
    
    for k = 1:3
        % get d_enter, d_exit
        d_enter = zeros(r,1);
        d_exit = zeros(r,1);
        for i = 1:r
            d_enter(i) = X((Npoly+1)*(index-1)-r+i,k) * factorial(r-i);
        end
        if index+1 < size(ts,2)
            for i = 1:r
                d_exit(i) = X((Npoly+1)*(index+1)-r+i,k) * factorial(r-i);
            end
        else
            temp = X((Npoly+1)*(index-1)+1:(Npoly+1)*index,k);
            for i = 1:r
                d_exit(i) = polyval(temp,ts(index+1)-ts(index));
                temp = polyder(temp);
            end
        end
        
        a = C\[d_enter;d_exit];
        
        for b = 0:1
            valid = 0;
            boundary = zeros(Npoly+1,1);
            boundary(end) = dup(index-1,k+3*b);

            % get boundary_time
            cand = roots(a - boundary);
            for j = 1:size(cand,1)
                if isreal(cand(j)) && -margin < cand(j) && cand(j) < ts(index+1)-ts(index-1)+margin
                    valid = valid + 1;
                    boundary_time_xyz(k,b+1) = ts(index-1)+cand(j);
                end
            end
            if valid < 1
                boundary_time_xyz(k,b+1) = -1;
            elseif valid > 1
                error('fun');            
            end
        end
        
        if boundary_time_xyz(k,:) ~= -1
            continue
        elseif boundary_time_xyz(k,:) == -1
            boundary_time_xyz(k,:) = [ts(index-1) ts(index+1)];
            continue
        end
        if dup(index-1,k) < d_enter(r) && d_enter(r) < dup(index-1,k+3)
            boundary_time_xyz(k,boundary_time_xyz(k,:) == -1) = ts(index-1);
        elseif dup(index-1,k) < d_exit(r) && d_exit(r) < dup(index-1,k+3)
            boundary_time_xyz(k,boundary_time_xyz(k,:) == -1) = ts(index+1);
        else
            error('funfun')
        end
    end
    result = sort(boundary_time_xyz,2);
end

