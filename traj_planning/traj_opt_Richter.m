function [p, N, cost] = traj_opt_Richter(path, map, ts)
% Polynomial Trajectory Planning
% This code is based on below paper.
% Polynomial Trajectory Planning for Aggressive Quadrotor Flight in Dense
% Indoor Environments
% Charles Richter, Adam Bry, and Nicholas Roy
%
% This code is written by
% Jungwon Park
% Seoul National University
% wjstk335@gmail.com
%
% input:  path                      waypoints x outdim 
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
    N = 5; % order of trajectory planning, N must be odd, N > 3
    [waypoints,outdim] = size(path); % waypoints should be waypoints > 2 
    M = waypoints - 1; % number of segments
    
    Coef = zeros(1,N); % 1xN coefficient matrix of cost function J(T)
                       % CAUTION! C(0) not considered
    
    % Set C by optimization purpose                 
%     Coef(4) = 1;       % minimum snap trajectory
    Coef(3) = 1;     % minimum jerk trajectory
    %Coef(2) = 1;     % minimum accelation trajectory
    
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
            Q_temp = flip(flip(Q_temp,1),2); % flip for lowest coeff go first
            
            Q((N+1)*(m-1)+1 : (N+1)*m, (N+1)*(m-1)+1 : (N+1)*m) = ...
                Q((N+1)*(m-1)+1 : (N+1)*m, (N+1)*(m-1)+1 : (N+1)*m) + Coef(n) * Q_temp;
        end
    end
    
    %% Build mapping matrix A
    % A : (N+1)M x (N+1)M    
    % The mapping matrix has the following structure:
    % if N = 5,
    % [ x 0 0 0 0 0 ]
    % [ 0 x 0 0 0 0 ]
    % [ 0 0 x 0 0 0 ]
    % [ x x x x x x ]
    % [ 0 x x x x x ]
    % [ 0 0 x x x x ]
    % from https://github.com/ethz-asl/mav_trajectory_generation
    
    A = zeros((N+1)*M, (N+1)*M);
    
    A_0 = zeros((N+1)/2,N+1);
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
        A_T = zeros((N+1)/2, N+1);
        for j = 1:size(A_T,2)
            for i = 1:size(A_T,1)
                if(i > j)
                    A_T(i,j) = 0;
                else
                    A_T(i,j) = factorial(j-1)/factorial(j-i) * (ts(m+1)-ts(m))^(j-i);
                end
            end
        end 
       
        A((N+1)*(m-1)+1 : (N+1)*m, (N+1)*(m-1)+1 : (N+1)*m) = [A_0; A_T];
    end
    
    % Get Ainv
    % Ainv : (N+1)M x (N+1)M
    Ainv = inv(A); % ethz-asl used schur-complement method to get A_inv but in my case I use simplest one
    
    %% Build reordering matrix C, C_transpose
    % C : (N+1)(M+1)/2 x (N+1)M
    % C_transpose : (N+1)M x (N+1)(M+1)/2
    C_transpose_1 = zeros((N+1)*M, (N+1)*(M+1)/2);
    C_transpose_2 = zeros((N+1)*M, (N+1)*(M+1)/2);
    
    % reorder fixed startpoint 
    C_transpose_1(1:(N+1)/2, 1:(N+1)/2) = eye((N+1)/2,(N+1)/2);
    iter = (N+1)/2+1;
    
    % reorder fixed mid point
    for m = 1 : M-1
        C_transpose_1((N+1)/2*(2*m-1)+1,iter) = 1;
        C_transpose_2((N+1)/2*(2*m)+1, iter) = 1;
        iter = iter + 1;
    end

    % reorder fixed endpoint
    for j = (N+1)/2*(2*M-1)+1 : (N+1)*M
        C_transpose_1(j,iter) = 1;
        iter = iter + 1;
    end
    
    % reorder free constraints    
    for m = 1 : M-1
        dm_T_start = (N+1)*(2*m-1)/2+1;
        dm_0_start = (N+1)*(2*m)/2+1;
        C_transpose_1( dm_T_start+1 : dm_T_start+(N-1)/2 , iter:iter+(N-1)/2-1) = eye((N-1)/2, (N-1)/2);
        C_transpose_2( dm_0_start+1 : dm_0_start+(N-1)/2 , iter:iter+(N-1)/2-1) = eye((N-1)/2, (N-1)/2);
        iter = iter + (N-1)/2;
    end
    
    C = C_transpose_1';
    C_transpose = C_transpose_1 + C_transpose_2;
    
    
    %% Build derivative vectors d 
    % d :  (N+1)M x outdim
    d = zeros((N+1)*M, outdim);
    for k = 1 : outdim 
        for m = 1 : M
            d((N+1)*(m-1)+1, k) =  path(m,k);
            d((N+1)/2*(2*m-1)+1, k) =  path(m+1,k);
        end
    end
    
    %% Build R and get d_opt
    % R : (N+1)M x (N+1)M
    % R = C*(A^-T)*Q*(A^-1)*(C^T); 
    H = Ainv' * Q * Ainv;
    R = C_transpose' * H * C_transpose;
    
    d_f_end = N+M;
    d_p_start = N+M+1;
    R_pp = R(d_p_start:end, d_p_start:end);
    R_fp = R(1:d_f_end, d_p_start:end);
    d_opt = zeros((N+1)*M, outdim);
    
    % d_p_opt = -R_pp\((R_fp')*d_f);
    % To invert R_pp, use QR decomposition
    for k = 1 : outdim 
        d_fp = C*d(:,k);
        d_f = d_fp(1:d_f_end);
        
%         uptriang = dsp.UpperTriangularSolver;
%         [QQ,u] = qr(R_pp);
%         b = (-QQ'*R_fp'*d_f);
%         d_p_opt = uptriang(u,b);
        d_p_opt = -inv(R_pp)*R_fp'*d_f;
        d_opt(:,k) = C_transpose * [d_f; d_p_opt];
    end
    
    %% Solve p
    % p: M*(N+1) x 1
    p = zeros(M*(N+1), outdim);
    for k = 1 : outdim 
        p(:,k) = A\d_opt(:,k);
    end
    
    cost = 0;
    for k = 1 : outdim 
        cost = cost + p(:,k)'*Q*p(:,k);
    end
        
    % Flip p in order [t^N; t^N-1; ... t^2; t; 1] for consistency
    for m = 1 : M
        p((N+1)*(m-1)+1:(N+1)*m,:) = flip(p((N+1)*(m-1)+1:(N+1)*m,:));
    end

end
