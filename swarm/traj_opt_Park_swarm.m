function [p, N, total_cost] = traj_opt_Park_swarm(path, box_cell, ts_cell, ts_total, rel, margin)
% This code is written by
% Jungwon Park
% Seoul National University
% Intelligent Control Systems Laboratory
% wjstk335@gmail.com
%
% input:  path                      cell{ path_size x xyz }
%         box_cell                  cell{ boxpath_size x box }
%         ts_cell                   cell{ array: boxpath_size+1 }
%         ts_total                  array, total time segment
%         rel                       matrix[ rel_size x rel_info ]
%         margin                    inter collision margin
%
% output: p                         matrix[ M*(N+1) x outdim ]
%                                   Optimized coefficient vectors of time vector [t^N; t^N-1; ... t^2; t; 1]
%                                   i.e. x_m(t) = p((N+1)*(m-1)+1:(N+1)*k,1)' * [t^N; t^N-1; ... t^2; t; 1]
%                                   outdim indicates number of output (x,y,z)
%
%         N                         scalar, order of trajectory planning
%         total_cost                scalar,

    disp('Building Constraints matrices ...');
    tic
    %% Initialization    
    N = 5;             % order of trajectory planning
    M = size(ts_total,2)-1;   % number of total segments 
    qn = size(path,1);
    
    outdim = 3;        % only for x,y,z
%     box_margin = 0.001;    % enforcing box constraint margin
    
    Coef = zeros(1,N); % 1xN coefficient matrix of cost function J(T)
                       % CAUTION! you cannot touch C(0)
    
    % Set C by optimization purpose                 
    r = 3;
    Coef(r) = 1;
    % Coef(4) = 1;     % minimum snap trajectory
    % Coef(3) = 1;     % minimum jerk trajectory
    % Coef(2) = 1;     % minimum accelation trajectory

    gao_rescale = 0;   % multipling scale factor on traj
    dw = 1;            % downwash effect coefficient
                       
    %% Build cost function Q
    % Q : qn(N+1)M x qn(N+1)M
    
    % build cost matrix for one quad
    Q_ = zeros(M*(N+1),M*(N+1));
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
        
        for m = 1:M
            Qp_ = Qp * (ts_total(m+1)-ts_total(m))^(-2*n+1 + gao_rescale*2); 
                
            Q_((N+1)*(m-1)+1 : (N+1)*m, (N+1)*(m-1)+1 : (N+1)*m) = ...
                Q_((N+1)*(m-1)+1 : (N+1)*m, (N+1)*(m-1)+1 : (N+1)*m) + Coef(n) * Qp_;
        end
    end
    
    % build total cost matrix Q
    Q = Q_;
    for qi = 1:qn-1
        Q = blkdiag(Q,Q_);
    end        
    
    %% Build mapping matrix Aeq
    % qn(2*init_const + (M-1)*cont_const) x qn(N+1)M
    init_const = r;
    cont_const = r;
    
    A_waypoints = zeros((2*init_const), (N+1)*M);
    A_cont = zeros((M-1)*cont_const, (N+1)*M);
    
    A_0 = zeros(N+1);
    A_T = zeros(N+1);
    
    tn = 1;
    for i = 0 : N
        A_0(1+i, 1:1+i) = tn;
        A_T(1+i, end-i:end) = tn;
        tn = conv(tn,[-1 1]);
    end
    
    % Build A_waypoints
    nn = 1;
    for i = 1:init_const
        A_waypoints(i, 1:N+1) = (ts_total(2)-ts_total(1))^(gao_rescale-(i-1)) * nn * A_0(i, :);
        A_waypoints(init_const+i, (N+1)*(M-1)+1:(N+1)*M) = (ts_total(end)-ts_total(end-1))^(gao_rescale-(i-1)) * nn * A_T(i, :);
        nn = nn * (N-i+1);
    end
    
    % Build A_cont
    for m = 1 : M-1
        nn = 1;
        for j = 1:cont_const   
            A_cont(cont_const*(m-1)+j, (N+1)*(m-1)+1 : (N+1)*m) = (ts_total(m+1)-ts_total(m))^(gao_rescale-(j-1)) * nn * A_T(j, :);
            A_cont(cont_const*(m-1)+j, (N+1)*m+1 : (N+1)*(m+1)) = -(ts_total(m+2)-ts_total(m+1))^(gao_rescale-(j-1)) * nn * A_0(j, :);
            nn = nn * (N-j+1);
        end
    end
    
    Aeq_ = [A_waypoints; A_cont];
    
    % build total mapping matrix Aeq
    Aeq = Aeq_;
    for qi = 1:qn-1
        Aeq = blkdiag(Aeq,Aeq_);
    end
    
    %% Build mapping matrix Alq
    % Alq : qn2(N+1)M + qn(qn-1)(N+1)M x qn(N+1)M
    
    % build Alq_box
    Alq_box_ = [eye((N+1)*M); -eye((N+1)*M)];
    Alq_box = Alq_box_;
    for qi = 1:qn-1
        Alq_box = blkdiag(Alq_box,Alq_box_);
    end
    
    % build Alq_rel
    Alq_rel = []; % zeros(qn*(qn-1)*(N+1)*M, qn*(N+1)*M);
    iter = 0;
    for qi = 1:qn
        for qj = qi+1:qn
            Alq_rel((N+1)*M*2*iter+1 : (N+1)*M*(2*iter+1), (N+1)*M*(qi-1)+1:(N+1)*M*qi) = -eye((N+1)*M);
            Alq_rel((N+1)*M*2*iter+1 : (N+1)*M*(2*iter+1), (N+1)*M*(qj-1)+1:(N+1)*M*qj) = eye((N+1)*M);
            
            Alq_rel((N+1)*M*(2*iter+1)+1 : (N+1)*M*2*(iter+1), (N+1)*M*(qi-1)+1:(N+1)*M*qi) = eye((N+1)*M);
            Alq_rel((N+1)*M*(2*iter+1)+1 : (N+1)*M*2*(iter+1), (N+1)*M*(qj-1)+1:(N+1)*M*qj) = -eye((N+1)*M);
            
            iter = iter + 1;
        end
    end
    
    % build total mapping matrix Alq
    Alq = [Alq_box; Alq_rel];
    
    %% Build derivative vectors deq
    % deq :  qn(2*init_const+(M-1)*init_const_cont) x outdim
    deq = [];
    
    for qi = 1:qn
        d_waypoints = zeros(2*init_const, outdim);
        d_cont = zeros((M-1)*cont_const, outdim);
        for k = 1 : outdim
            d_waypoints(1,k) = path{qi}(1,k);
            d_waypoints(init_const+1,k) = path{qi}(end,k);
        end
        
        deq_ = [d_waypoints; d_cont];
        deq = [deq; deq_];
    end
    
    %% Build derivative vectors dlq 
    % dlq :  qn2(N+1)M + qn(qn-1)(N+1)M x outdim
    dlq_box = [];
    for qi = 1:qn
        d_upper = zeros((N+1)*M, outdim);
        d_lower = zeros((N+1)*M, outdim);

        for k = 1 : outdim
            bi = 1;
            for m = 1 : M
                for n = 1 : N+1
                    % find box number
                    while ts_cell{qi}(bi+1) < ts_total(m+1) %%?
                        bi = bi + 1;
                    end
                    d_upper((N+1)*(m-1)+n, k) = box_cell{qi}(bi,k+3)/((ts_total(m+1)-ts_total(m))^gao_rescale);
                    d_lower((N+1)*(m-1)+n, k) = -box_cell{qi}(bi,k)/((ts_total(m+1)-ts_total(m))^gao_rescale);
                end
            end
        end

        dlq_box_ = [d_upper; d_lower];
        dlq_box = [dlq_box; dlq_box_];
    end
    
    dlq_rel = [];
    for qi = 1:qn
        for qj = qi+1:qn
            d_upper = zeros((N+1)*M, outdim);
            d_lower = zeros((N+1)*M, outdim);
            
            % find box number
            for k = 1 : outdim
                for m = 1 : M
                    for n = 1 : N+1
                        % find box number
                        ri = 1;
                        while rel(ri,1) ~= qi || rel(ri,2) ~= qj || rel(ri,4) < ts_total(m+1)
                            ri = ri + 1;
                        end
                        
                        sector = rel(ri,3);
                        
                        d_upper((N+1)*(m-1)+n, k) = 10000000;
                        d_lower((N+1)*(m-1)+n, k) = 10000000;
                        if abs(sector) == k
%                             if sector < 0 && k ~= 3
%                                 d_upper((N+1)*(m-1)+n, k) = -margin/((ts_total(m+1)-ts_total(m))^gao_rescale);
%                             elseif sector < 0 && k == 3
% %                                 d_upper((N+1)*(m-1)+n, k) = -dw*margin/((ts_total(m+1)-ts_total(m))^gao_rescale);
%                             elseif sector > 0 && k ~= 3
%                                 d_lower((N+1)*(m-1)+n, k) = -margin/((ts_total(m+1)-ts_total(m))^gao_rescale);
%                             elseif sector > 0 && k == 3
% %                                 d_lower((N+1)*(m-1)+n, k) = -dw*margin/((ts_total(m+1)-ts_total(m))^gao_rescale);
%                             end
                            if sector < 0
                                d_upper((N+1)*(m-1)+n, k) = -2*margin/((ts_total(m+1)-ts_total(m))^gao_rescale);
                            elseif sector > 0
                                d_lower((N+1)*(m-1)+n, k) = -2*margin/((ts_total(m+1)-ts_total(m))^gao_rescale);
                            end
                        end
                    end
                end
            end
            
            dlq_rel_ = [d_upper; d_lower];
            dlq_rel = [dlq_rel; dlq_rel_];
        end
    end
    
    dlq = [dlq_box; dlq_rel];
    
    
%     save const_matrix
    
    toc
    %% Solve QP and get c
    % c: M*(N+1) x 1
    disp('Solve QP and get coefficients ...');
    tic
    c = zeros(qn*M*(N+1), outdim);
    total_cost = 0;
    
    for k = 1 : outdim
%       % quadprog
%         options = optimoptions('quadprog','MaxIterations',200,'Display','none');
%         options = optimoptions('quadprog','MaxIterations',200);
%         [c(:,k), cost, exitflag] = quadprog(Q,[],Alq,dlq(:,k),Aeq,deq(:,k),[],[],[],options);
        % cplex
        [c(:,k), cost, exitflag] = cplexqp(Q,zeros(1,qn*M*(N+1)),Alq,dlq(:,k),Aeq,deq(:,k));
        
        if cost < -1e-4
            error('cost < 0');
        end
%         if exitflag ~= 1
%             exitflag
%             error('quadprog error')
%         end
        total_cost = total_cost + cost;
    end   
    toc
    
    %% translate c to p
    p = cell(qn,1);
    for qi = 1:qn
        p{qi} = zeros(M*(N+1), outdim);

        for k = 1:outdim
            for m = 1:M
                for i = 1:N+1
                    p{qi}((N+1)*(m-1)+1:(N+1)*m, k) = p{qi}((N+1)*(m-1)+1:(N+1)*m, k) + ...
                        ((ts_total(m+1)-ts_total(m))^gao_rescale)*c((qi-1)*(N+1)*M+(N+1)*(m-1)+i,k)*(b(i,:).*timevector(1/(ts_total(m+1)-ts_total(m)),N,0))';
                end
            end
        end
    end
    
   