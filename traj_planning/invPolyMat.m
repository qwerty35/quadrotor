function [detA,A_pinv] = invPolyMat(A)
% This code is based on below paper.
% Computation of the Generalized Inverse of a Polynomial Matrix and
% Applications
% N.P.Karampetakis
%
% This code is written by
% Jungwon Park
% Seoul National University
% wjstk335@gmail.com
%
% input:  A                 cell{ A_0, A_1, ... , A_q }
%
% output: detA              array 1 x k_terminate+1
%         A_pinv            cell{ A_pinv0, A_pinv1, ... , A_pinv(2*k_terminate-1)*q }
%
%         pseudo inverse of A = (detA.*[s^0 s^1 s^2 ... s^q])^(-1) *
%             (A_pinv0*s^0 + A_pinv1*s^1 + ... + A_pinv0*s^((2*k_terminate-1)*q))

q = size(A,2)-1;
n = size(A{1},1);
m = size(A{1},2);

a_hat = zeros(n+1,2*n*q+1);
B = cell(n,2*(n-1)*q+1);

% Initial and Boundary conditions
B{0+1,0+1} = eye(n);
for i = 0:n-1
    for j  = 2*i*q+1:2*(n-1)*q
        B{i+1,j+1} = 0;
    end
end
a_hat(1,1) = 1;

% Recursive relations for a_hat, B
k_terminate = n;
prev_terminate_flag = 0;
for i = 0:n-1
    
    terminate_flag = 1;
    
    for j = 0:2*(i+1)*q
        temp = 0;
        for k = 0:j
            for l = 0:j-k
                if j-k-l > q || l > q || k > 2*(n-1)*q
                    continue;
                end
                temp = temp + A{j-k-l+1}*A{l+1}'*B{i+1,k+1};
            end
        end
        a_hat(i+1+1,j+1) = -trace(temp)/(i+1);
        
        if i < n-1
            B{i+1+1,j+1} = temp + a_hat(i+1+1,j+1)*eye(n);
        end
        
        if a_hat(i+1+1,j+1) ~= 0
            terminate_flag = 0;
            prev_terminate_flag = 0;
        end
    end
    
    % Terminate
    if terminate_flag == 1 && prev_terminate_flag == 0
        k_terminate = i;
        prev_terminate_flag = 1;
    end
end

% Output
detA = a_hat(k_terminate+1,:);
A_pinv = cell(1,(2*k_terminate-1)*q+1);
for j = 0:(2*k_terminate-1)*q
    A_pinv{j+1} = 0;
    for l = 0:j
        if j-l > q || l > 2*(k_terminate-1)*q
            continue;
        end
        A_pinv{j+1} = A_pinv{j+1} - A{j-l+1}'*B{k_terminate-1+1,l+1};
    end
end

end