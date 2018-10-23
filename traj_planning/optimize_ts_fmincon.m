function [X,Npoly,ts] = optimize_ts_fmincon(path, map, ts_init, traj_opt_handle)
%
% optimize time segment by using fmincon
%
% This code is written by
% Jungwon Park
% Seoul National University
% wjstk335@gmail.com
%
    ts_size = size(ts_init,2);
    Alq = eye(ts_size-1,ts_size);
    for i = 1:ts_size-1
        Alq(i,i+1) = -1;
    end
    blq = -ones(ts_size-1,1)*0.1;
    Aeq = zeros(2,ts_size);
    Aeq(1,1) = 1;
    Aeq(2,ts_size) = 1;
    beq = [ts_init(1);ts_init(end)];
    lb = ts_init(1);
    ub = ts_init(end);
    [ts, ~,~] = fmincon(@(ts)traj_opt_wrapper(path,map,ts,traj_opt_handle),ts_init,Alq,blq,Aeq,beq,lb,ub);
    [X,Npoly,cost] = traj_opt_handle(path,map,ts);
    ts(1) = 0;
    ts
    cost
end