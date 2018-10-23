function [p, N, cost] = traj_opt_nonlinear(path, box, ts)
% This code is based on below paper.
% Online Safe Trajectory Generation For Quadrotors Using Fast Marching
% Method and Bernstein Basis Polynomial
% Fei Gao, William Wu, Yi Lin and Shaojie Shen
%
% This code is written by
% Jungwon Park
% Seoul National University
% wjstk335@gmail.com
%
% input:  path                      {waypoints x outdim} x # of quadrotor
%
%         box                       matrix, # of path * 9([x0,y0,z0,x1,y1,z1,colorR,G,B ]) safe flight corridor 
%
%         ts                        {time segments} x # of quadrotor
%
% output: p                         {M*(N+1) x outdim} x # of quadrotor
%                                   Optimized coefficient vectors of time vector
%                                   [t^N; t^N-1; ... t^2; t; 1]
%                                   x_m(t) = p((N+1)*(m-1)+1:(N+1)*k,1)' * [t^N; t^N-1; ... t^2; t; 1]
%                                   outdim indicates number of output (x,y,z)
%
%         N                         order of trajectory planning

    %% Initialization    
    N = 5;             % order of trajectory planning
    M = size(box{1},1);   % number of segments 
    qn = size(box,2);
    
    count = 0;
    
    outdim = 3;        % only for x,y,z
    margin = 0.001;    % enforcing box constraint margin
    collision_radius = 0.4; % enforcing collision radius constraint
    
    Coef = zeros(1,N); % 1xN coefficient matrix of cost function J(T)
                       % CAUTION! you cannot touch C(0)
    
    % Set C by optimization purpose                 
    r = 3;
    Coef(r) = 1;
    % Coef(4) = 1;     % minimum snap trajectory
    % Coef(3) = 1;     % minimum jerk trajectory
    % Coef(2) = 1;     % minimum accelation trajectory

    gao_rescale = 0;   % multipling scale factor on traj
                       % 1 : true
                       % 0 : false

    
    %% Build cost function Q
    % Q : (N+1)M x (N+1)M    
    
    for n = 1:N     % derivative number
        if Coef(n) == 0
            continue;
        end
            
        % derivative b(t)
        db = zeros(N+1-n);
        b = zeros(N+1);
        for i = 0 : N
            ti = zeros(1,i+1);
            ti(1) = 1;
                
            tn = 1;
            for j = 1 : N-i
                tn = conv(tn,[-1 1]);
            end
            bi = nchoosek(N,i)*conv(ti,tn);
            b(i+1,:) = bi;
                
            for k = 1:n
                bi = polyder(bi);
            end
            db(i+1,:) = bi;
        end
            
        % integration
        Qp = zeros(N+1);
        for i = 1:N+1
            for j = 1:N+1
                db_temp = conv(db(i,:), db(j,:));
                db_temp = polyint(db_temp);
                Qp(i,j) = diff(polyval(db_temp,[0 1]));
            end
        end
        
        Q_total = [];
        for qi = 1:qn
            Q = zeros(M*(N+1),M*(N+1));
            ts_i = ts{qi};
            for m = 1:M     % segment number  
                Q((N+1)*(m-1)+1 : (N+1)*m, (N+1)*(m-1)+1 : (N+1)*m) = ...
                    Q((N+1)*(m-1)+1 : (N+1)*m, (N+1)*(m-1)+1 : (N+1)*m) + Coef(n) * Qp * (ts_i(m+1)-ts_i(m))^(-2*n+1 + gao_rescale*2);
            end
            Q = blkdiag(Q,Q,Q);
            Q_total = blkdiag(Q_total,Q);
        end
    end
    
    
    
    %% Build mapping matrix Aeq
    % 2*init_const + (M-1)*cont_const x (N+1)M
    init_const = r;
    cont_const = r;
    
    A_0 = zeros(N+1);
    A_T = zeros(N+1);
    
    tn = 1;
    for i = 0 : N
        A_0(1+i, 1:1+i) = tn;
        A_T(1+i, end-i:end) = tn;
        tn = conv(tn,[-1 1]);
    end
    
    Aeq_total = [];
    for qi = 1:qn
        ts_i = ts{qi};
        A_waypoints = zeros((2*init_const), (N+1)*M);
        A_cont = zeros((M-1)*cont_const, (N+1)*M);
        
        % Build A_waypoints
        nn = 1;
        for i = 1:init_const
            A_waypoints(i, 1:N+1) = (ts_i(2)-ts_i(1))^(gao_rescale-(i-1)) * nn * A_0(i, :);
            A_waypoints(init_const+i, (N+1)*(M-1)+1:(N+1)*M) = (ts_i(end)-ts_i(end-1))^(gao_rescale-(i-1)) * nn * A_T(i, :);
            nn = nn * (N-i+1);
        end

        % Build A_cont
        for m = 1 : M-1
            nn = 1;
            for j = 1:cont_const   
                A_cont(cont_const*(m-1)+j, (N+1)*(m-1)+1 : (N+1)*m) = (ts_i(m+1)-ts_i(m))^(gao_rescale-(j-1)) * nn * A_T(j, :);
                A_cont(cont_const*(m-1)+j, (N+1)*m+1 : (N+1)*(m+1)) = -(ts_i(m+2)-ts_i(m+1))^(gao_rescale-(j-1)) * nn * A_0(j, :);
                nn = nn * (N-j+1);
            end
        end

        Aeq = [A_waypoints; A_cont];
        Aeq = blkdiag(Aeq,Aeq,Aeq);
        Aeq_total = blkdiag(Aeq_total, Aeq);
    end
    
    %% Build mapping matrix Alq
    % Alq : 2(N+1)M*outdim x (N+1)M*outdim
    Alq = [eye((N+1)*M); -eye((N+1)*M)];
    Alq = blkdiag(Alq,Alq,Alq);
    
    Alq_total = [];
    for qi = 1:qn
        Alq_total = blkdiag(Alq_total, Alq);
    end
    
    %% Build deq 
    % deq :  (2*init_const+(M-1)*init_const_cont)*outdim x 1 
    deq_total = [];
    for qi = 1:qn
        path_i = path{qi};
        d_waypoints = zeros(2*init_const, outdim);
        d_cont = zeros((M-1)*cont_const, outdim);

        for k = 1 : outdim
            d_waypoints(1,k) = path_i(1,k);
            d_waypoints(init_const+1,k) = path_i(end,k); 
        end

        deq = [d_waypoints; d_cont];
        deq = deq(:);
        deq_total = [deq_total; deq];
    end
    
    %% Build dlq 
    % dlq :  2(N+1)M*outdim x 1
    dlq_total = [];
    for qi = 1:qn
        ts_i = ts{qi};
        box_i = box{qi};
        d_upper = zeros((N+1)*M, outdim);
        d_lower = zeros((N+1)*M, outdim);

        for k = 1 : outdim
            for m = 1 : M
                for n = 1 : N+1
                    d_upper((N+1)*(m-1)+n, k) = box_i(m,k+3)/((ts_i(m+1)-ts_i(m))^gao_rescale);
                    d_lower((N+1)*(m-1)+n, k) = -box_i(m,k)/((ts_i(m+1)-ts_i(m))^gao_rescale);
                end
            end
        end

        dlq = [d_upper; d_lower];
        dlq = dlq(:);
        dlq_total = [dlq_total; dlq];        
    end
    
    %% Solve QP
    % c: M*(N+1) x outdim
    c = zeros(M*(N+1)*outdim*qn,1);

    options = optimoptions('fmincon','MaxFunctionEvaluations',100000);
    [c, cost, exitflag] = fmincon(@(x)traj_opt_cost(x,Q_total,ts,M,N,outdim,qn,1),c,Alq_total,dlq_total,Aeq_total,deq_total,[],[],[],options);
    cost
%         if exitflag ~= 1
%             exitflag
%             error('quadprog error')
%         end
    c = reshape(c,[M*(N+1),outdim,qn]);
    
    %% translate c to p
    p = {};
    for qi = 1:qn
        ts_i = ts{qi};
        p_temp = zeros(M*(N+1), outdim);
        
        for k = 1:outdim
            for m = 1:M
                for i = 1:N+1
                    p_temp((N+1)*(m-1)+1:(N+1)*m, k) = p_temp((N+1)*(m-1)+1:(N+1)*m, k) + ...
                        ((ts_i(m+1)-ts_i(m))^gao_rescale)*c((N+1)*(m-1)+i,k,qi)*(b(i,:).*timevector(1/(ts_i(m+1)-ts_i(m)),N,0))';
                end
            end
        end
        
        p = [p,p_temp];
    end
    
    
   