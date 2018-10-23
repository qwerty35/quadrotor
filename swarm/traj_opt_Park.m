function [p, N, total_cost] = traj_opt_Park(path, map, ts)
%
% This code is written by
% Jungwon Park
% Seoul National University
% wjstk335@gmail.com
%

% CAUTION safe traj is not implemented!!!!

% input:  path                      waypoints x outdim 
%         map                       
%         ts                        time segments
%
% output: p                         M*(N+1) x outdim 
%                                   Optimized coefficient vectors of time vector
%                                   [t^N; t^N-1; ... t^2; t; 1]
%                                   x_m(t) = p((N+1)*(m-1)+1:(N+1)*k,1)' * [t^N; t^N-1; ... t^2; t; 1]
%                                   outdim indicates number of output (x,y,z)
%
%         N                         order of trajectory planning
%         total_cost                total QP cost

    %% Initialization
    N = 5; % order of trajectory planning, N should be odd, N > 3
    [waypoints,outdim] = size(path); % waypoints should be waypoints > 2 
    M = waypoints - 1; % number of segments
    
    % Set C by optimization purpose
    Coef = zeros(1,N); % 1xN coefficient matrix of cost function J(T)
    
    r = 3;
    Coef(r) = 1;
    % Coef(4) = 1;     % minimum snap trajectory
    % Coef(3) = 1;     % minimum jerk trajectory
    % Coef(2) = 1;     % minimum accelation trajectory
    
    % Set number of constraints
    init_const = r;    % initialize initial condition only pos, vel, acc
    
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
    
    %% Build mapping matrix A
    % A : (2*init_const+M-1)+cont_const*(M-1) x (N+1)M
%     if init_const > N+1
%         error('init_const is TOO BIG!');
%     end
    
%     A_waypoints = zeros((2*init_const+M-1), (N+1)*M);
%     A_cont = zeros(init_const*(M-1), (N+1)*M);
%     
%     A_0 = zeros(init_const,N+1);
%     for j = 1:size(A_0,2)
%         for i = 1:size(A_0,1)
%             if(i == j)
%                 A_0(i,j) = factorial(j-1)/factorial(j-i);
%             else
%                 A_0(i,j) = 0;
%             end
%         end
%     end
%     
%     A_0_1 = zeros(1,N+1);
%     A_0_1(1) = 1;
%     
%     A_waypoints(1:init_const, 1:N+1) = A_0(1:init_const, :);
%     for m = 1 : M
%         A_T = zeros(init_const, N+1);
%         for j = 1:size(A_T,2)
%             for i = 1:size(A_T,1)
%                 if(i > j)
%                     A_T(i,j) = 0;
%                 else
%                     A_T(i,j) = factorial(j-1)/factorial(j-i) * (ts(m+1)-ts(m))^(j-i);
%                 end
%             end
%         end 
%         
%         % Build A_waypoints
%         if m > 1
%             A_waypoints(init_const+m-1, (N+1)*(m-1)+1:(N+1)*m) = A_0_1;
%         end
%         if m == M
%             A_waypoints(init_const+M:2*init_const+M-1, (N+1)*(M-1)+1:(N+1)*M) = A_T;
%         end
%         
%         % Build A_cont
%         if m < M
%             A_cont(init_const*(m-1)+1 : init_const*m, (N+1)*(m-1)+1 : (N+1)*m) = A_T;
%             A_cont(init_const*(m-1)+1 : init_const*m, (N+1)*m+1 : (N+1)*(m+1)) = -A_0;  
%         end
%     end
%     
%     A = [A_waypoints; A_cont];

    A = zeros(M-1+init_const*(M+1), (N+1)*M);
    
    A_0 = zeros(init_const,N+1);
    for j = 1:size(A_0,2)
        for i = 1:size(A_0,1)
            if(i == j)
                A_0(i,j) = factorial(j-1)/factorial(j-i);
            else
                A_0(i,j) = 0;
            end
        end
    end
    
    A_0_1 = zeros(1,N+1);
    A_0_1(1) = 1;
    
    A(1:init_const, 1:N+1) = A_0(1:init_const, :);
    for m = 1 : M
        A_T = zeros(init_const, N+1);
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
        if m > 1
            A(init_const+m-1, (N+1)*(m-1)+1:(N+1)*m) = A_0_1;
        end
        if m == M
            A(init_const+M:2*init_const+M-1, (N+1)*(M-1)+1:(N+1)*M) = A_T;
        end
        
        % Build A_cont
        if m < M
            A(M+init_const*(m+1) : M-1+init_const*(m+2), (N+1)*(m-1)+1 : (N+1)*m) = A_T;
            A(M+init_const*(m+1) : M-1+init_const*(m+2), (N+1)*m+1 : (N+1)*(m+1)) = -A_0;  
        end
    end

    %% Build derivative vectors d 
    % d :  (2*init_const+M-1)+init_const_cont*(M-1) x outdim
%     d_waypoints = zeros(2*init_const+M-1, outdim);
%     d_cont = zeros(init_const*(M-1), outdim);
%     
%     for k = 1 : outdim
%         d_waypoints(1,k) = path(1,k);
%         for m = 2 : M+1
%             d_waypoints(init_const+m-1, k) =  path(m,k);
%         end
%     end
%     
%     d = [d_waypoints; d_cont];
    
    d = zeros(2*init_const+M-1+init_const*(M-1), outdim);    
    for k = 1 : outdim
        d(1,k) = path(1,k);
        for m = 2 : M+1
            d(init_const+m-1, k) =  path(m,k);
        end
    end
    
    %% Solve QP and get p
    % p: M*(N+1) x 1    
    p = zeros(M*(N+1), outdim);
    
    Q_ihat = inv(Q+A'*A);
    total_cost = 0;
    
    for k = 1 : outdim 
%         p(:,k) = Q_ihat*A'*inv(A*Q_ihat*A') * d(:,k);
        p(:,k) = Q_ihat * (A'*((A*Q_ihat*A') \ d(:,k)));
%         total_cost = total_cost + d(:,k)'*(inv(A*Q_ihat*A')-eye(size(A,1)))*d(:,k);
    end
        
    % Flip p in order [t^N; t^N-1; ... t^2; t; 1] for consistency
    for m = 1 : M
        p((N+1)*(m-1)+1:(N+1)*m,:) = flip(p((N+1)*(m-1)+1:(N+1)*m,:));
    end
end
