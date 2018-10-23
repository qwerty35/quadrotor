function [X,Npoly,ts] = optimize_ts_grad(path, map, ts, traj_opt_handle)
%
% optimize time segment by using backtracking line search
%
% This code is written by
% Jungwon Park
% Seoul National University
% wjstk335@gmail.com
%

    h = 1e-5;
    end_condition = 1e-7;
    T = ts(2:end) - ts(1:end-1);
    m = size(T,2);
    iter = 1;
    count = 0;
    
    F = inf;
    F_update = inf;
    deltaF = 1;
    alpha = 0.1;
    beta = 0.3;
    while deltaF > end_condition 
        iter = iter + 1;
        
        % get gradient of F
        dF = zeros(1,m);
        ts = [0 cumsum(T)];
        [X,Npoly,F] = traj_opt_handle(path,map,ts);
%         if F < 40.2 
%             F
%             break
%         end
        
        g = (ones(m,m)-eye(m,m))*(-1/(m-1)) + eye(m,m);
        for i = 1:m
            ts_i = [0 cumsum(T+h*g(i,:))];
            [~,~,F_i] = traj_opt_handle(path,map,ts_i);
            if F_i < inf
                dF = dF + (F_i-F)/h * g(i,:);
            end
        end
        count = count + m+1;
        
        gamma = 1;
        % Before line search, find feasible gamma
        while true          
            ts_update = [0 cumsum(T-gamma*dF)];
            [~,index] = sort(ts_update);
            
            if any(index ~= (1:m+1))
                gamma = beta*gamma;
            else
                break;
            end
        end
        % Backtracking line search ->
        % linesearch·Î?%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        while true
            accept = abs(T-gamma*dF) > 0.05;
            ts_update = [0 cumsum(T-gamma*dF.*accept)];
            [~,~,F_update] = traj_opt_handle(path,map,ts_update);
            count = count + 1;
            if F_update > F - alpha*gamma*(norm(dF)^2)
                gamma = beta*gamma;
            else
                break;
            end
        end
        T = T - gamma * dF;
        ts = ts_update;
        deltaF = F - F_update;
    end
    
    iter
    count
    ts
    dF
    gamma
    F_update
end