function [K,dK] = getK(ts)
%% Initialization
    M = size(ts,2)-1;
    N = 5;
    r= 3;
    % Set C by optimization purpose
    Coef = zeros(1,N); % 1xN coefficient matrix of cost function J(T)
    
    Coef(r) = 1;
    % Coef(4) = 1;     % minimum snap trajectory
    % Coef(3) = 1;     % minimum jerk trajectory
    % Coef(2) = 1;     % minimum accelation trajectory
    
    % Set number of constraints
    init_const = r;    % initialize initial condition only pos, vel, acc
    cont_const = r;    % Acceptable Jump?
    
       %% Build cost function Q
        % Q : (N+1)M x (N+1)M
        Q = zeros(M*(N+1),M*(N+1));
        dQ = zeros(M*(N+1),M*(N+1));
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
        if init_const > N+1
            error('init_const is TOO BIG!');
        end

        A_waypoints = zeros((2*init_const+M-1), (N+1)*M);
        A_cont = zeros(cont_const*(M-1), (N+1)*M);
        
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

        A_0_1 = zeros(1,N+1);
        A_0_1(1) = 1;

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
            dA_T = zeros(N+1, N+1);
            for j = 1:size(dA_T,2)
                for i = 1:size(dA_T,1)
                    if(i > j)
                        dA_T(i,j) = 0;
                    else
                        dA_T(i,j) = factorial(j-1)/factorial(j-i) * (ts(m+1)-ts(m))^(j-i);
                    end
                end
            end 

            % Build A_waypoints
            A_waypoints(1:init_const, 1:N+1) = A_0(1:init_const, :);
            if m > 1
                A_waypoints(init_const+m-1, (N+1)*(m-1)+1:(N+1)*m) = A_0_1;
            end
            if m == M
                A_waypoints(init_const+M:2*init_const+M-1, (N+1)*(M-1)+1:(N+1)*M) = A_T(1:init_const, :);
          end

            % Build A_cont
            if m < M
                A_cont(cont_const*(m-1)+1 : cont_const*m, (N+1)*(m-1)+1 : (N+1)*m) = A_T(1:cont_const, :);
                A_cont(cont_const*(m-1)+1 : cont_const*m, (N+1)*m+1 : (N+1)*(m+1)) = -A_0(1:cont_const, :);
            end
        end
        
        A = [A_cont;A_waypoints];
        
        K = [Q A'; A zeros(size(A,1))];
        
       dQ1 = [0 0 0 0 0 0;
              0 0 0 0 0 0;
              0 0 0 0 0 0;
              0 0 0 36 144*(ts(2)-ts(1)) 360*(ts(2)-ts(1))^2;
              0 0 0 144*(ts(2)-ts(1)) 576*(ts(2)-ts(1))^2 1440*(ts(2)-ts(1))^3;
              0 0 0 360*(ts(2)-ts(1))^2 1440*(ts(2)-ts(1))^3 3600*(ts(2)-ts(1))^4];
       dQ2 = [0 0 0 0 0 0;
              0 0 0 0 0 0;
              0 0 0 0 0 0;
              0 0 0 36 -144*(ts(3)-ts(2)) 360*(ts(3)-ts(2))^2;
              0 0 0 -144*(ts(3)-ts(2)) 576*(ts(3)-ts(2))^2 -1440*(ts(3)-ts(2))^3;
              0 0 0 360*(ts(3)-ts(2))^2 -1440*(ts(3)-ts(2))^3 3600*(ts(3)-ts(2))^4];
        dQ = [dQ1 zeros(6);
              zeros(6) dQ2]; 
          
        dA_cont = [0 1 2*(ts(2)-ts(1)) 3*(ts(2)-ts(1))^2 4*(ts(2)-ts(1))^3 5*(ts(2)-ts(1))^4;
              0 0 2 6*(ts(2)-ts(1)) 12*(ts(2)-ts(1))^2 20*(ts(2)-ts(1))^3;
              0 0 0 6 24*(ts(2)-ts(1)) 60*(ts(2)-ts(1))^2];
   
        dA_waypoint = [0 1 -2*(ts(3)-ts(2)) 3*(ts(3)-ts(2))^2 -4*(ts(3)-ts(2))^3 5*(ts(3)-ts(2))^4;
                       0 0 2 -6*(ts(3)-ts(2)) 12*(ts(3)-ts(2))^2 -20*(ts(3)-ts(2))^3;
                       0 0 0 6 -24*(ts(3)-ts(2)) 60*(ts(3)-ts(2))^2];
        dA = [dA_cont zeros(3,6);
              zeros(4,12);
              zeros(3,6) dA_waypoint];
        dK = [dQ dA'; dA zeros(size(dA,1))];
end

