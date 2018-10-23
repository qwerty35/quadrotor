function [p, N, total_cost] = traj_opt_Chen(path, map, ts)
% This code is based on below paper.
% Online Generation of Collision-Free Trajectories for Quadrotor Flight in
% Unknown Cluttered Environment
% Jing Chen, Tianbo Liu, and Shaojie Shen
%
% This code is written by
% Jungwon Park
% Seoul National University
% wjstk335@gmail.com
%
% input:  path                      waypoints x outdim 
%
%         map                       map should have box, dup
%         box                       matrix #path * 9 ([x0,y0,z0,x1,y1,z1,colorR,G,B ]) safe region 
%         dup                       matrix #path-1 * 9 ([x0,y0,z0,x1,y1,z1,colorR,G,B]) duplicated
%                                   volume between boxes
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
    box = map.box;
    dup = map.dup;
    
    N = 5;            % order of trajectory planning
    M = size(box,1);   % number of segments 
    
    outdim = 3;        % trajectory only for x,y,z
    init_const = 3;    % initialize only pos, vel, acc
    margin = 0.001;    % enforcing box constraint margin
    
    Coef = zeros(1,N); % 1xN coefficient matrix of cost function J(T)
                       % CAUTION! you cannot touch C(0)
    
    % Set C by optimization purpose                 
    r = 3;
    Coef(r) = 1;
    % Coef(4) = 1;       % minimum snap trajectory
    % Coef(3) = 1;     % minimum jerk trajectory
    % Coef(2) = 1;     % minimum accelation trajectory


    

    %% Build cost function Q
    % Q : (N+1)M x (N+1)M
    Q = zeros(M*(N+1),M*(N+1));
    for m = 1:M     % segment number
        for n = 1:N     % derivative number
            if Coef(n) == 0
                continue;
            end
            
            % derivative P(t)
            q = ones(1,N+1);
            for idx = 1:n
                q = polyder(q);
            end
            Q_temp = q'*q;
            
            % integration
            for i = 1:size(Q_temp,1)
                for j = 1:size(Q_temp,2)
                    k = size(Q_temp,1)+size(Q_temp,2)-i-j+1;
                    Q_temp(i,j) = Q_temp(i,j)/k * (ts(m+1)-ts(m))^k;
                end
            end
            
            Q_temp(N+1,N+1) = 0; % rescale
            Q_temp = flip(flip(Q_temp,1),2); % lowest coeff go first
            
            Q((N+1)*(m-1)+1 : (N+1)*m, (N+1)*(m-1)+1 : (N+1)*m) = ...
                Q((N+1)*(m-1)+1 : (N+1)*m, (N+1)*(m-1)+1 : (N+1)*m) + Coef(n) * Q_temp;
        end
    end
%     Q = Q + eye(M*(N+1))*reg; % regularize term
%     cond(Q)
    
    %% Build mapping matrix Aeq
    % Aeq : 2*init_const + (M-1)*init_const x (N+1)M
    init_const_cont = r;
    if init_const > N+1
        error('init_const is TOO BIG!');
    end
    A_waypoints = zeros((2*init_const), (N+1)*M);
    A_cont = zeros((M-1)*init_const_cont, (N+1)*M);
    
    A_0 = zeros(N+1,N+1);
    for j = 1:size(A_0,2)
        for i = 1:size(A_0,1)
            if(i == j)
                A_0(i,j) = factorial(j-1)/factorial(j-i);
            else
                A_0(i,j) = 0;
            end
        end
    end
    
    for m = 1 : M
        A_T = zeros(N+1, N+1);
        for j = 1:size(A_T,2)
            for i = 1:size(A_T,1)
                if(i > j)
                    A_T(i,j) = 0;
                else
                    A_T(i,j) = factorial(j-1)/factorial(j-i) * (ts(m+1)-ts(m))^(j-i);
                end
            end
        end 
        
        % Build A_waypoints
        if m == 1
            A_waypoints(1:init_const, 1:N+1) = A_0(1:init_const, :);
        elseif m == M
            A_waypoints(init_const+1:2*init_const, (N+1)*(M-1)+1:(N+1)*M) = A_T(1:init_const, :);
        end
        
        % Build A_cont
        if m < M
            A_cont(init_const_cont*(m-1)+1 : init_const_cont*m, (N+1)*(m-1)+1 : (N+1)*m) = A_T(1:init_const_cont, :);
            A_cont(init_const_cont*(m-1)+1 : init_const_cont*m, (N+1)*m+1 : (N+1)*(m+1)) = -A_0(1:init_const_cont, :);  
        end
    end
    
    Aeq = [A_waypoints; A_cont];
    
    %% Build mapping matrix Alq
    % Alq : 2*M x (N+1)M
    Alq_partial = zeros(M, M*(N+1));
    
    A_0_1 = zeros(1,N+1);
    A_0_1(1) = 1;
    
    Alq_partial(1, 1:N+1) = zeros(1,N+1);
    for m = 2 : M        
        % Build A_waypoints
        Alq_partial(m, (N+1)*(m-1)+1:(N+1)*m) = A_0_1;
    end
    
    Alq = [Alq_partial; -Alq_partial];
    
    %% Build derivative vectors deq 
    % deq :  2*init_const+(M-1)*init_const x outdim
    d_waypoints = zeros(2*init_const, outdim);
    d_cont = zeros((M-1)*init_const_cont, outdim);
    
    for k = 1 : outdim
        d_waypoints(1,k) = path(1,k);
        d_waypoints(init_const+1,k) = path(end,k); 
    end
    
    deq = [d_waypoints; d_cont];
    
    %% Build derivative vectors dlq 
    % dlq :  2*M x outdim
    d_upper = zeros(M, outdim);
    d_lower = zeros(M, outdim);
    
    for k = 1 : outdim
        d_upper(1,k) = 1;
        d_lower(1,k) = 1;
        for m = 1 : M-1
            d_upper(m+1,k) = dup(m,k+3);
            d_lower(m+1,k) = -dup(m,k);
        end
    end
    
    dlq = [d_upper; d_lower];
    
    %% Solve QP and get p
    % p: M*(N+1) x 1
    p = zeros(M*(N+1), outdim);
    total_iter = 0; % number of QP execution
    total_cost = 0;
    
    for k = 1 : outdim 
        isfeasible = 0;
        Alq_k = Alq;
        dlq_k = dlq(:,k);
        while isfeasible == 0    
            options = optimoptions('quadprog','MaxIterations',200);
            [p(:,k), cost, exitflag] = quadprog(Q,[],Alq_k,dlq_k,Aeq,deq(:,k),[],[],[],options);
            
            if exitflag ~= 1
                cost = inf;
                break;
            end
            
            % Flip p in order [t^N; t^N-1; ... t^2; t; 1] for consistency
            for m = 1 : M
                p((N+1)*(m-1)+1:(N+1)*m,k) = flip(p((N+1)*(m-1)+1:(N+1)*m,k));
            end

            % Enforcing Box Constraints
            isfeasible = 1;
            for m = 1 : M
                p_m = p((N+1)*(m-1)+1:(N+1)*m,k);
                t_e = roots(polyder(p_m));
                real_index = abs(imag(t_e)) < 1e-13;
                t_e = t_e(real_index);
                for i = 1 : size(t_e, 1)
                    if 0 < t_e(i) && t_e(i) < ts(m+1)-ts(m) 
                        f_m = polyval(p_m, t_e(i));
                        if f_m < box(m,k)
                            Alq_partial = zeros(1,M*(N+1));
                            Alq_partial(1,(N+1)*(m-1)+1:(N+1)*m) = flip(timevector(t_e(i),N,0));

                            Alq_k = [Alq_k; -Alq_partial];
                            dlq_k = [dlq_k; -box(m,k)-margin];
                            
                            isfeasible = 0;
                            break;
                        elseif f_m > box(m,k+3)
                            Alq_partial = zeros(1,M*(N+1));
                            Alq_partial(1,(N+1)*(m-1)+1:(N+1)*m) = flip(timevector(t_e(i),N,0));
                            
                            Alq_k = [Alq_k; Alq_partial];
                            dlq_k = [dlq_k; box(m,k+3)-margin];
                            
                            isfeasible = 0;
                            break;
                        end
                    end    
                end
            end
            total_iter = total_iter + 1;
        end
        total_cost = total_cost + cost;
    end
    total_iter
    
    