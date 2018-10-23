function [Qp,Ap] = QpAp(a,b,c,d,r)
    %% Initialization
    M = 2; % number of segments
    N = 2*r-1;
    
    % Set C by optimization purpose
    Coef = zeros(1,N); % 1xN coefficient matrix of cost function J(T)

    Coef(r) = 1;
    % Coef(4) = 1;     % minimum snap trajectory
    % Coef(3) = 1;     % minimum jerk trajectory
    % Coef(2) = 1;     % minimum accelation trajectory
    
        for n = 1:N     % derivative number
            if Coef(n) == 0
                continue;
            end
            
            % derivative P(t)
            q = ones(1,N+1);
            for idx = 1:n
                q = polyder(q);
            end
            
            % integration
            Qp = zeros(N+1);
            for i = 1:(N+1)/2
                for j = 1:(N+1)/2
                    q_i = zeros(1,size(q,2));
                    q_j = zeros(1,size(q,2));
                    q_i(i) = q(i);
                    q_j(j) = q(j);
                    q_temp = polyint(conv(q_i,q_j));
                    
                    Qp(i,j) = diff(polyval(q_temp,[a, b]));
                end
            end
        end
        
       Qp = flip(flip(Qp,1),2);
        
       Ap = zeros(N+1);
       for i = 1:r
            Ac_temp = flip(timevector(c,N,i-1),2);
            Ad_temp = flip(timevector(d,N,i-1),2);
            Ap(i,:) = Ac_temp;
            Ap(i+r,:) = Ad_temp;
       end
       
        
end

