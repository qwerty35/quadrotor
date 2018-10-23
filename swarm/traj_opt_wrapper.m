function cost = traj_opt_wrapper(path, map, ts, traj_opt_handle)
    [~,~,cost] = traj_opt_handle(path,map,ts);
end
    
    
   