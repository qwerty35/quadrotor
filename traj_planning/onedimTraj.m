clear
clc
close all


M = 2;
r = 3;

start = 0.1;
interval = 0.01;
fin = 0.9;
j =1;

% result = [];
% for t = 0.1:0.01:0.9
%     [K,dK] = getK([0 t 1]);
%     Kinv = inv(K);
%     c = [zeros(15,1); 0; 2; 0; 1; 2; 2; 0];
%     result  = [result [c'*Kinv*dK*Kinv*c; t]];
% %     result  = [result [-c'*Kinv*c; t]];
% end

convex=[];
nonconvex=[];
while j < 2500

      c = [[0;0;0];
           [0; 0; 0;
            j;
            0; 0; 0;]; ] % convex?
%     c = [0;0;0;10*(rand(3,1)-0.5);0;10*(rand(3,1)-0.5);]

    cost = zeros(2,(fin-start)*(1/interval)+1);
    index = 1;
    for i = start:interval:fin
        ts = [0 i 1];
        W = getW(ts,5,3);
        cost(:,index) = [c'*W*c; i];
        index = index + 1;
    end
    dcost = cost(1,2:end) - cost(1,1:end-1);
    ddcost = dcost(1,2:end) - dcost(1,1:end-1);
%     if any(ddcost < 0)
        plot(cost(2,:),cost(1,:));
        hold on;
%         break;
%         nonconvex = [nonconvex c];
%     else
%         convex = [convex c];
%     end
%     
    j = j+1;
end