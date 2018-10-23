function sdot = quadEOM_CO(t, s, qi, controlhandle, trajhandle, params, map)
% revised from https://github.com/damanfb/LQR-Obstacles, quadEOM

OBSTACLE_STEPS = 45;
NUM_POINTS = 80;
dt = 1/10;

persistent expdtAB points s_dim u_dim qd
if numel(t) == 0 | numel(s) == 0
    F_norminal = params.mass * params.grav / 4;
    s_refm = [0;0;0;0;0;0;1;0;0;0;0;0;0;F_norminal*4;0;0];
    u = [F_norminal*4;0;0];
    s_dim = size(s_refm,1);
    u_dim = size(u,1);

    [A,B] = linearDiscretize(s_refm, u, params);
    AB = zeros(s_dim+u_dim+1);
    AB(1:s_dim,1:s_dim) = A;
    AB(1:s_dim, s_dim+1:s_dim+u_dim) = B;
    expdtAB = expm(dt*AB); 
    
    points = create_sphere(NUM_POINTS,1/2,3/2); 
    return
end
    tic
    
    C = [1 0 0, 0 0 0, 0 0 0 0, 0 0 0, 0 0 0;
         0 1 0, 0 0 0, 0 0 0 0, 0 0 0, 0 0 0;
         0 0 1, 0 0 0, 0 0 0 0, 0 0 0, 0 0 0;];
    bTr = [                    1,                      1,                     1,                      1;
                               0,      params.arm_length,                     0,     -params.arm_length;
              -params.arm_length,                      0,     params.arm_length,                      0;
           params.k_M/params.k_F, -params.k_M/params.k_F, params.k_M/params.k_F, -params.k_M/params.k_F];
    rTb = [0.25,                      0, -0.5/params.arm_length,  0.25*params.k_F/params.k_M;
           0.25,  0.5/params.arm_length,                      0, -0.25*params.k_F/params.k_M;
           0.25,                      0,  0.5/params.arm_length,  0.25*params.k_F/params.k_M;
           0.25, -0.5/params.arm_length,                      0, -0.25*params.k_F/params.k_M];
    
    
    orca_matrix = [];
    orca_constraint = [];
    
    % convert state to quad stuct for control
    qd{qi} = stateToQd(s);

    % Get desired_state
    desired_state = trajhandle(t, qi);

    % The desired_state is set in the trajectory generator
    qd{qi}.pos_des      = desired_state.pos;
    qd{qi}.vel_des      = desired_state.vel;
    qd{qi}.acc_des      = desired_state.acc;
    qd{qi}.yaw_des      = desired_state.yaw;
    qd{qi}.yawdot_des   = desired_state.yawdot;
    
    % reformat s
    curr_w = s(14:17);
    curr_thrusts = params.k_F*(curr_w.^2);
    si = [s(1:13); bTr(1:3,:)*curr_thrusts];
    qd{qi}.s_refm = si;
    
    qn = size(qd,2);
    for qj = 1:qn
        if qj==qi
            continue
        else
            ellipsoids = [];          
            if isfield(qd{qj},'s_refm')
                sj = qd{qj}.s_refm;
            else
                continue
            end
            
            expAB = eye(s_dim+u_dim+1);
            for k = 1:OBSTACLE_STEPS
                expAB = expAB * expdtAB;
                F = expAB(1:s_dim, 1:s_dim);
                G = expAB(1:s_dim, s_dim+1:s_dim+u_dim);
                J = inv(C*G);
                for l = 1:NUM_POINTS
                    reachable_point = J*(points(:,l) + C*F*(sj-si));
                    reachable_thrust = rTb(:,1:3)*reachable_point;
                    if ones(4,1) * (params.minF-params.maxF) <= reachable_thrust & reachable_thrust <= ones(4,1) * (params.maxF - params.minF)
                        ellipsoids = [ellipsoids; reachable_point'];
                    end
                end
            end
            if size(ellipsoids,1)>4
                curr_input = (sj(end-2:end)-si(end-2:end))';
                K = convhull(ellipsoids);
                [in, dist, normal] = inhull(curr_input, ellipsoids, K);
                if in
                else
                    S1.Vertices = ellipsoids; 
                    S1.Faces = K;
                    S1Obj = patch(S1);
                    S2.Vertices = curr_input;
                    S2Obj = patch(S2);
                    [dist,~,G,~]=GJK_dist(S1Obj,S2Obj);
                    normal = G'-curr_input;
                    normal = normal./sqrt(sum(normal.^2, 2));
                end
                orca_matrix =  [orca_matrix; -normal];
                orca_constraint = [orca_constraint; -dot(si(end-2:end)' - dist/2 * normal, normal)]; % +-?
            end
        end
    end
    
    % get pref control outputs
    [F_pref, M_pref, ~, ~] = controlhandle(qd, t, qi, params);
    pref_thrusts = rTb * [F_pref; M_pref];
    
    options = optimoptions('lsqlin','Display','off');
    [x,~,~,exitflag,~,~] = lsqlin(rTb(:,1:3), pref_thrusts, [orca_matrix; rTb(:,1:3)], [orca_constraint; ones(4,1) * params.maxF],[],[],[],[],[],options);
%     [x,~,~,exitflag,~,~] = lsqlin(eye(3), [F_pref; M_pref(1:2)], orca_matrix, orca_constraint,[],[],[],[],[],options);
    
    if exitflag ~= 1
        error("no sol");
    end
    
    F_new = x(1);
    M_new = [x(2:3); 0;]; 
    

    % compute derivative
    sdot = quadEOM_readonly(t, s, F_new, M_new, params, map);
    toc
end
