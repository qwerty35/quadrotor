close all 
clear
clc

[Winv,dWinv,ddWinv] = getWinv(3);

% syms T t
% ddWinvT1 = subs(ddWinv,T,1);
% result = [];
% for i = 0.01:0.01:0.99
%     asdf = double(subs(ddWinvT1,t,i));
%     result = [result [sort(eig(asdf)); i]];
% end


% syms T t p1 p2 v a positive
syms p1 p2 T t  v1 v2 a1 a2 alpha real
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

 
 
syms t2 real
Q2 = [0 0 0 0 0 0 0 0 0 0 0 0;
     0 0 0 0 0 0 0 0 0 0 0 0;
     0 0 0 0 0 0 0 0 0 0 0 0;
     0 0 0 36*t 72*t^2 120*t^3 0 0 0 0 0 0;
     0 0 0 72*t^2 192*t^3 360*t^4 0 0 0 0 0 0;
     0 0 0 120*t^3 360*t^4 720*t^5 0 0 0 0 0 0;
     0 0 0 0 0 0 0 0 0 0 0 0;
     0 0 0 0 0 0 0 0 0 0 0 0;
     0 0 0 0 0 0 0 0 0 0 0 0;
     0 0 0 0 0 0 0 0 0 36*t2 72*t2^2 120*t2^3;
     0 0 0 0 0 0 0 0 0 72*t2^2 192*t2^3 360*t2^4;
     0 0 0 0 0 0 0 0 0 120*t2^3 360*t2^4 720*t2^5;];
A2 = [1 t t^2 t^3   t^4    t^5   -1  0  0 0 0 0;
     0 1 2*t 3*t^2 4*t^3  5*t^4  0 -1  0 0 0 0;
     0 0 2   6*t   12*t^2 20*t^3 0  0 -2 0 0 0;
     1 0 0 0 0 0 0 0 0 0 0 0;
     0 1 0 0 0 0 0 0 0 0 0 0;
     0 0 2 0 0 0 0 0 0 0 0 0;
     0 0 0 0 0 0 1 0 0 0 0 0;
     0 0 0 0 0 0 1 t2 t2^2 t2^3 t2^4 t2^5;
     0 0 0 0 0 0 0 1 2*t2 3*t2^2 4*t2^3 5*t2^4;
     0 0 0 0 0 0 0 0 2 6*t2 12*t2^2 20*t2^3;];
 


Qp0t = [0 0 0 0 0 0;
      0 0 0 0 0 0;
      0 0 0 0 0 0;
      0 0 0 36*t 72*t^2 120*t^3;
      0 0 0 72*t^2 192*t^3 360*t^4;
      0 0 0 120*t^3 360*t^4 720*t^5;];
QptT = [0 0 0 0 0 0;
      0 0 0 0 0 0;
      0 0 0 0 0 0;
      0 0 0 36*T 72*T^2 120*T^3;
      0 0 0 72*T^2 192*T^3 360*T^4;
      0 0 0 120*T^3 360*T^4 720*T^5;] - Qp0t;
Ap = [1 0 0 0 0 0;
      0 1 0 0 0 0;
      0 0 2 0 0 0;
      1 t t^2 t^3 t^4 t^5;
      0 1 2*t 3*t^2 4*t^3 5*t^4;
      0 0 2 6*t 12*t^2 20*t^3];
Apt0 = [1 t t^2 t^3 t^4 t^5;
      0 1 2*t 3*t^2 4*t^3 5*t^4;
      0 0 2 6*t 12*t^2 20*t^3;
      1 0 0 0 0 0;
      0 1 0 0 0 0;
      0 0 2 0 0 0;];
AptT = [1 t t^2 t^3 t^4 t^5;
      0 1 2*t 3*t^2 4*t^3 5*t^4;
      0 0 2 6*t 12*t^2 20*t^3;
      1 T T^2 T^3 T^4 T^5;
      0 1 2*T 3*T^2 4*T^3 5*T^4;
      0 0 2 6*T 12*T^2 20*T^3;];
eyee = [[0 0;1 0;0 1;]; zeros(3,2)];



syms t3 real

Q3 = t;
Q3(1:6,1:6) = Qp0t;
Q3(7:12,7:12) = subs(Qp0t,t,t2);
Q3(13:18,13:18) = subs(Qp0t,t,t3);

A3 = [1 t t^2 t^3 t^4 t^5 -1 0  0 0 0 0 0 0 0 0 0 0;
     0 1 2*t 3*t^2 4*t^3 5*t^4 0 -1  0 0 0 0 0 0 0 0 0 0;
     0 0 2 6*t 12*t^2 20*t^3 0 0 -2 0 0 0 0 0 0 0 0 0;
     0 0 0 0 0 0 1 t2 t2^2 t2^3 t2^4 t2^5 -1 0  0 0 0 0;
     0 0 0 0 0 0 0 1 2*t2 3*t2^2 4*t2^3  5*t2^4 0 -1  0 0 0 0;
     0 0 0 0 0 0 0 0 2 6*t2 12*t2^2 20*t2^3 0 0 -2 0 0 0;
     1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
     0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
     0 0 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
     0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0;
     0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0;
     0 0 0 0 0 0 0 0 0 0 0 0 1 t3 t3^2 t3^3 t3^4 t3^5;
     0 0 0 0 0 0 0 0 0 0 0 0 0 1 2*t3 3*t3^2 4*t3^3 5*t3^4;
     0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 6*t3 12*t3^2 20*t3^3;];
 
K3 = [Q3 A3'; A3 zeros(14)];
K3inv = inv(K3);
W3 = K3inv(end-7:end,end-7:end);

C = [1 0 0 0 0 0;
     0 1 0 0 0 0;
     0 0 2 0 0 0;
     1 t t^2 t^3 t^4 t^5;
     0 1 2*t 3*t^2 4*t^3 5*t^4;
     0 0 2 6*t 12*t^2 20*t^3;];
 
asdf = (inv(C))' * Qp0t * inv(C); 

% simplify(QptT*inv(AptT))*eyee
% Qp0t*inv(Apt0)*eyee

% K = K{1}+K{2}*t + K{3}*t^2+ K{4}*t^3+K{5}*t^4+K{6}*t^5;
% K = [Q A'; A zeros(10)];
% % K2 = [Q2 A2'; A2 zeros(10)];
% 
% Kinv = inv([eye(12) zeros(12,10); zeros(10,12) eye(10)*(-1)] * K);
% K2inv = inv([eye(12) zeros(12,10); zeros(10,12) eye(10)*(-1)] * K2);
% W = K2inv(16:22,16:22);
% cW = collect(W*(t+t2)^3);
% M = diag([1 t^-1 t^-2 1 1 t2^-1 t2^-2]);
% M2 = diag([1 t^-1 t^-2 1 1 t2^-1 t2^-2]);
% collect(diag([1 t^-1 t^-2 1 1 t2^-1 t2^-2])*diag([1 t^-1 t^-2 1 1 t2^-1 t2^-2])*M*cW*M*diag([1 t^-1 t^-2 1 1 t2^-1 t2^-2])*diag([1 t^-1 t^-2 1 1 t2^-1 t2^-2]))
% cost = [];
% for i = 0.1:0.01:0.9
%     cost = [cost [[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 1 0 0] * subs(subs(Kinv,T,1),t,i) * [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 1 0 0]'; i]];
% end
% plot(cost(2,:), cost(1,:));
% hold on
% fplot([0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 1 0 0] * subs(Kinv,T,1) * [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 1 0 0]');


% 
% dK = diff(Kinv(16:22,16:22),t);
% % dKT1 = subs(dK,T,1);
% ddK = diff(dK,t);
% ddKT1 = subs(ddK,T,1);

% c = +-[p1;v1;a1;0;p2;-v2;a2]; %% 하나빼고 전부 + 계수 -> 계산시 무조건 convex
% c = [-p1;v1;a1;0;p2;v2;-a2];
% 
% asdf = ddK*t^7*(T-t)^7*T^3;


% % syms vv1 vv2 aa1 aa2
% asdf = subs(asdf,v1,v1/(t));
% asdf = subs(asdf,v2,v2/(t2));
% asdf = subs(asdf,a1,a1/(t^2));
% asdf = subs(asdf,a2,a2/(t2^2));


% result = [];
% c_max = [inf;inf;inf;inf;inf;inf;inf;];
% c_min = -c_max;
% for i = 0.1:0.01:0.9
%     [V,D] = eig(double(subs(ddKT1,t,i)));
%     c_cand = [];
%     for alpha1 = -23.2 : 0.1 : 23.2 
%         c_cand = [c_cand V*[alpha1;0;0;0;0;0;1]];
%     end
%     for j = 1:7
%         if max(c_cand(j,:)) < c_max(j)
%             c_max(j) = max(c_cand(j,:));
%         end
%         if min(c_cand(j,:)) > c_min(j)
%             c_min(j) = min(c_cand(j,:));
%         end
%     end
% %     alpha = inv(V)*c;
% %     result = [result [V(:,1); abs(diag(D)/D(1,1)); alpha.^2; i]];
% end
% for j = 1:size(result,1)-1
%     figure;
%     plot(result(end,:), result(j,:))
% end
