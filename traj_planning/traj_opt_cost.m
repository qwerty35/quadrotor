function [cost] = traj_opt_cost(x,Q_total,ts,M,N,outdim,qn,margin)
%     heuristic?
%     rx = reshape(x,[M*(N+1),outdim,qn]);
% %     delta = 0.1./(x(1:M*(N+1)*outdim)-x(M*(N+1)*outdim+1:M*(N+1)*outdim*2)+1);
% 
%     cost = x'*Q_total*x+delta'*delta;

    dt = 0.1;    
    rx = reshape(x,[M*(N+1),outdim,qn]);
    collision_cost = 0;
    lambda_cost = 10000;
    
    ts1 = ts{1};
    for t = 0:dt:ts1(end)
        for qi = 1:qn
            ts_i = ts{qi};
            m = find(ts_i<=t);
            m = m(end);
            pi = rx((m-1)*(N+1)+1:m*(N+1),:,qi);
            
            t_b = (t-ts_i(m))/(ts_i(m+1)-ts_i(m));
            pos_i = zeros(1,outdim);
            for i = 0:N
                pos_i = pos_i + pi(i+1,:)*nchoosek(N,i)*(t_b)^i*(1-t_b)^(N-i);
            end
            for qj = qi+1:qn
                ts_j = ts{qj};
                m = find(ts_j<=t);
                m = m(end);
                pj = rx((m-1)*(N+1)+1:m*(N+1),:,qj);
                
                t_b = (t-ts_j(m))/(ts_j(m+1)-ts_j(m));
                pos_j = zeros(1,outdim);
                for j = 0:N
                    pos_j = pos_j + pj(j+1,:)*nchoosek(N,j)*(t_b)^j*(1-t_b)^(N-j);
                end
                
                dpos = pos_i - pos_j;
                dpos(3) = dpos(3)/300; % scale z-axis by 3 to make it ellipsoid
                dpos = norm(dpos);
                
                if dpos <= margin
                    collision_cost = collision_cost + lambda_cost*dt*(dpos-margin)^2;
                end
            end
        end
    end
    
    cost = x'*Q_total*x + collision_cost;
end

