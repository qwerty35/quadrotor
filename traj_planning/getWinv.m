function [Winv,dWinv,ddWinv] = getWinv(r)
%
% Get symbolic KKT matrix K for two segments problem
% ts = [0 t T]
%
if r ~= 3
    error("not yet...");
end
syms T t real
Q = [0 0 0 0 0 0 0 0 0 0 0 0;
     0 0 0 0 0 0 0 0 0 0 0 0;
     0 0 0 0 0 0 0 0 0 0 0 0;
     0 0 0 36*t 72*t^2 120*t^3 0 0 0 0 0 0;
     0 0 0 72*t^2 192*t^3 360*t^4 0 0 0 0 0 0;
     0 0 0 120*t^3 360*t^4 720*t^5 0 0 0 0 0 0;
     0 0 0 0 0 0 0 0 0 0 0 0;
     0 0 0 0 0 0 0 0 0 0 0 0;
     0 0 0 0 0 0 0 0 0 0 0 0;
     0 0 0 0 0 0 0 0 0 36*(T-t) 72*(T-t)^2 120*(T-t)^3;
     0 0 0 0 0 0 0 0 0 72*(T-t)^2 192*(T-t)^3 360*(T-t)^4;
     0 0 0 0 0 0 0 0 0 120*(T-t)^3 360*(T-t)^4 720*(T-t)^5;];
A = [1 t t^2 t^3   t^4    t^5   -1  0  0 0 0 0;
     0 1 2*t 3*t^2 4*t^3  5*t^4  0 -1  0 0 0 0;
     0 0 2   6*t   12*t^2 20*t^3 0  0 -2 0 0 0;
     1 0 0 0 0 0 0 0 0 0 0 0;
     0 1 0 0 0 0 0 0 0 0 0 0;
     0 0 2 0 0 0 0 0 0 0 0 0;
     0 0 0 0 0 0 1 0 0 0 0 0;
     0 0 0 0 0 0 1 (T-t) (T-t)^2 (T-t)^3 (T-t)^4 (T-t)^5;
     0 0 0 0 0 0 0 1 2*(T-t) 3*(T-t)^2 4*(T-t)^3 5*(T-t)^4;
     0 0 0 0 0 0 0 0 2 6*(T-t) 12*(T-t)^2 20*(T-t)^3;];

K = [Q A'; A zeros(size(A,1))];
Winv = inv([eye(size(Q)) zeros(size(A')); zeros(size(A)) eye(size(A,1))*(-1)] * K);

Winv = Winv(16:22,16:22);
dWinv = diff(Winv,t);
ddWinv = diff(dWinv,t);

end


