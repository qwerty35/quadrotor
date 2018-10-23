function [X,Npoly,ts] = optimize_ts_heu(path, map, ts, traj_opt_handle)
%
% optimize time segment by using backtracking line search
%
% This code is written by
% Jungwon Park
% Seoul National University
% wjstk335@gmail.com
%
    end_condition = 0.0001;
    iter = 0;
    F_prev = inf;
    margin = 1e-4;
    r = 3; % min jerk
    point_type = [];
    alpha = 0.7;
    
    syms T t real
    [~,dWinv,~] = getWinv(r);
    
    while true      
        [X,Npoly,F] = traj_opt_handle(path,map,ts);
        dF = F-F_prev
        F
        if abs(dF) < end_condition %|| dF > 0 % point_type 3,4 없을 때 까지?
            break
        end
        segment_points = [X(Npoly+1:Npoly+1:end,:); path(end,:)];
        
        dts = zeros(size(ts));
        time_temp = [];
        point_type = [];
        for i = 2:size(ts,2)-1
            [enter_time, exit_time] = boundary_point(X,Npoly,map.dup,ts,i);
            time_temp = [time_temp; [enter_time exit_time]];
            
%             % inner point?
%             if is_inner_point(segment_points(i,:), map.dup(i-1,:),0)
%                 point_type = [point_type 1];
%                 continue
%             end
%             % tangential point?
%             elseif  abs(enter_time - ts(i)) < margin && abs(exit_time - ts(i)) < margin
%                 % dK ddK로 local optima 
%                 point_type = [point_type 2];
%                 
%                 p = 0;
%                 for k = 1:3
%                     if i < size(ts,2)-1
%                         c = [X((Npoly+1)*(i-1),k); X((Npoly+1)*(i-1)-1,k); 0.5*X((Npoly+1)*(i-1)-2,k); 
%                              X((Npoly+1)*i,k);
%                              X((Npoly+1)*(i+1),k); X((Npoly+1)*(i+1)-1,k); 0.5*X((Npoly+1)*(i+1)-2,k);];
%                     else
%                         c = [X((Npoly+1)*(i-1),k); X((Npoly+1)*(i-1)-1,k); 0.5*X((Npoly+1)*(i-1)-2,k); 
%                              X((Npoly+1)*i,k);
%                              path(end,k);0;0;]; % zero to zero 
%                     end
%                     p = p + sym2poly(collect(c'*subs(dWinv*t^6*(T-t)^6,T,ts(i+1)-ts(i-1))*c,t));
%                     
%                     
%                 end
%                 convex = 0;
%                 t_cand = roots(p);
%                 for j = 1:size(t_cand,1)
%                     if isreal(t_cand(j)) && ts(i-1) < t_cand(j) && t_cand(j) < ts(i+1)
%                         convex = convex + 1;
%                         if convex > 1
%                             error("non convex!")
%                         end
%                         dts(i) = (t_cand(j)-ts(i));
%                     end
%                 end
%                 continue
%             else
%                 % boundary point
%                 point_type = [point_type 3];
                dts(i) = dts(i)+(exit_time-2*ts(i)+enter_time)*alpha;
%             end
        end
        
        % update
        F_prev = F;
        ts = ts + dts;
        
%         % avoid numerical error
%         for i = 1:size(segment_points,1)-2
%             if ts(i+1) - ts(i) < 0.05
%                 if ts(i+2) - ts(i+1) > 0.1
%                     ts(i+1) = ts(i)+0.05;
%                 elseif ts(i) - ts(i-1) > 0.1
%                     ts(i) = ts(i+1) - 0.05;
%                 else
%                     error("fun")
%                 end
%             end
%         end
        
        iter = iter + 1;
    end
    iter
    ts
    F
end