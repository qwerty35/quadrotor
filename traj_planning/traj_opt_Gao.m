function [p, N, total_cost,c] = traj_opt_Gao(path, box, ts)
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
% input:  path                      waypoints x outdim 
%
%         box                       matrix #path * 9 ([x0,y0,z0,x1,y1,z1,colorR,G,B ]) safe region 
%
%         ts                        time segments
%
% output: p                         M*(N+1) x outdim 
%                                   Optimized coefficient vectors of time vector
%                                   [t^N; t^N-1; ... t^2; t; 1]
%                                   x_m(t) = p((N+1)*(m-1)+1:(N+1)*k,1)' * [t^N; t^N-1; ... t^2; t; 1]
%                                   outdim indicates number of output (x,y,z)
%
%         N                         order of trajectory planning

    %% Initialization    
    N = 5;             % order of trajectory planning
    M_ori = size(box,1);   % number of segments 
    
    count = 0;
    redun = [];
    for m = 1:M_ori
        if ts(m+1)-ts(m) < 0.05
            count = count + 1;
            redun = [redun m];
        end
    end
    M = M_ori - count;
    box(redun,:) = [];
    ts(redun) = [];
    
    outdim = 3;        % only for x,y,z
    margin = 0.001;    % enforcing box constraint margin
    
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
    Q = zeros(M*(N+1),M*(N+1));
    
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
        
        for m = 1:M     % segment number  
            Qp_ = Qp * (ts(m+1)-ts(m))^(-2*n+1 + gao_rescale*2); 
                
            Q((N+1)*(m-1)+1 : (N+1)*m, (N+1)*(m-1)+1 : (N+1)*m) = ...
                Q((N+1)*(m-1)+1 : (N+1)*m, (N+1)*(m-1)+1 : (N+1)*m) + Coef(n) * Qp_;
        end
    end
    
    %% Build mapping matrix Aeq
    % 2*init_const + (M-1)*cont_const x (N+1)M
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
        A_waypoints(i, 1:N+1) = (ts(2)-ts(1))^(gao_rescale-(i-1)) * nn * A_0(i, :);
        A_waypoints(init_const+i, (N+1)*(M-1)+1:(N+1)*M) = (ts(end)-ts(end-1))^(gao_rescale-(i-1)) * nn * A_T(i, :);
        nn = nn * (N-i+1);
    end
    
    % Build A_cont
    for m = 1 : M-1
        nn = 1;
        for j = 1:cont_const   
            A_cont(cont_const*(m-1)+j, (N+1)*(m-1)+1 : (N+1)*m) = (ts(m+1)-ts(m))^(gao_rescale-(j-1)) * nn * A_T(j, :);
            A_cont(cont_const*(m-1)+j, (N+1)*m+1 : (N+1)*(m+1)) = -(ts(m+2)-ts(m+1))^(gao_rescale-(j-1)) * nn * A_0(j, :);
            nn = nn * (N-j+1);
        end
    end
    
    Aeq = [A_waypoints; A_cont];
    
    %% Build mapping matrix Alq
    % Alq : 2(N+1)M x (N+1)M
    Alq = [eye((N+1)*M); -eye((N+1)*M)];
        
    %% Build derivative vectors deq 
    % deq :  2*init_const+(M-1)*init_const_cont x outdim
    d_waypoints = zeros(2*init_const, outdim);
    d_cont = zeros((M-1)*cont_const, outdim);
    
    for k = 1 : outdim
        d_waypoints(1,k) = path(1,k);
        d_waypoints(init_const+1,k) = path(end,k); 
    end
    
    deq = [d_waypoints; d_cont];
    
    %% Build derivative vectors dlq 
    % dlq :  2(N+1)M x outdim
    d_upper = zeros((N+1)*M, outdim);
    d_lower = zeros((N+1)*M, outdim);
    
    for k = 1 : outdim
        for m = 1 : M
            for n = 1 : N+1
                d_upper((N+1)*(m-1)+n, k) = box(m,k+3)/((ts(m+1)-ts(m))^gao_rescale);
                d_lower((N+1)*(m-1)+n, k) = -box(m,k)/((ts(m+1)-ts(m))^gao_rescale);
            end
        end
    end
    
    dlq = [d_upper; d_lower];
    
    %% Solve QP and get c
    % c: M*(N+1) x 1
    c = zeros(M*(N+1), outdim);
    total_cost = 0;
    
    for k = 1 : outdim
%         options = optimoptions('quadprog','MaxIterations',200,'Display','none');
        options = optimoptions('quadprog','MaxIterations',200);
        [c(:,k), cost, exitflag] = quadprog(Q,[],Alq,dlq(:,k),Aeq,deq(:,k),[],[],[],options);
%         [c(:,k), cost, ~] = quadprog(Q,[],[],[],Aeq,deq(:,k),-d_lower(:,k),d_upper(:,k),[],options);
        if cost < -1e-4
            error('cost < 0');
        end
%         if exitflag ~= 1
%             exitflag
%             error('quadprog error')
%         end
        total_cost = total_cost + cost;
    end

    %% translate c to p
    p = zeros(M_ori*(N+1), outdim);

    for k = 1:outdim
        count = 0;
        for m = 1:M_ori
            if any(redun == m)
                count = count+1;
                continue
            end
            for i = 1:N+1
                p((N+1)*(m-1)+1:(N+1)*m, k) = p((N+1)*(m-1)+1:(N+1)*m, k) + ...
                    ((ts(m-count+1)-ts(m-count))^gao_rescale)*c((N+1)*(m-count-1)+i,k)*(b(i,:).*timevector(1/(ts(m-count+1)-ts(m-count)),N,0))';
            end
        end
    end
    
    
   