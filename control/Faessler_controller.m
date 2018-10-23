function[F_des, M_des, trpy, drpy] = Faessler_controller(qd, t, qn, params)
% Extracing parameters and references
% --------------------
% Dynamics of quadrotor suspended with load Constants
J = params.I;

g = [0; 0; -params.grav]; % convert to vector

r_ref = qd{qn}.pos_des;
v_ref = qd{qn}.vel_des;
a_ref = qd{qn}.acc_des;
psi_ref = 0;

% Extracing states
% ----------------
r = qd{qn}.pos;
v = qd{qn}.vel;
q = eul2quat(qd{qn}.euler','ZXY');%%%%%%%??
w = qd{qn}.omega;

% CONTROL
% -----------------
% Control parameter
% Position control parameter
k_p_xy = 1.5;
k_p_z = 15;
k_d_xy = 1.3;
k_d_z = 12;
        
kp = diag([k_p_xy,k_p_xy,k_p_z]);
kd = diag([k_d_xy,k_d_xy,k_d_z]);
    
% Attitude control parameter
k_p_rp = 5;
k_p_yaw = 1;
    
k_p_pq = 5;
k_p_r = 1;
    
k_p_att = diag([k_p_pq,k_p_pq,k_p_r]);
% -----------------
    
% Position Control
eB_z = quatrotate(q,[0 0 1])';
    
a_des = kp * (r_ref-r) + kd * (v_ref - v) + a_ref-g;
F_des = params.mass * dot(a_des, eB_z);
    
% Attitude Control
eB_z_des = a_des/norm(a_des);
    
temp = dot(eB_z,eB_z_des);
if abs(temp) > 1
    temp = sign(temp); 
end
    
alpha = acos(temp);
    
n_B = [0; 0; 0];
if(norm(cross(eB_z,eB_z_des)) ~= 0)
    n = cross(eB_z,eB_z_des)/norm(cross(eB_z,eB_z_des));
    n_B = quatrotate(quatconj(q),n')'; % body coordination of n
end
q_e_rp = [cos(alpha/2) sin(alpha/2)*n_B']; % error quaternion
    
w_p_des = sign(q_e_rp(1)) * 2 * k_p_rp * q_e_rp(2);
w_q_des = sign(q_e_rp(1)) * 2 * k_p_rp * q_e_rp(3);

eC_x = [cos(psi_ref) sin(psi_ref) 0]';
eC_y = [-sin(psi_ref) cos(psi_ref) 0]';
    
w_r_des = 0;
if(norm(cross(eC_y,eB_z_des)) ~= 0)
    eB_x_des = sign(eB_z_des(3))*cross(eC_y,eB_z_des)/norm(cross(eC_y,eB_z_des));
    eB_y_des = cross(eB_z_des,eB_x_des)/norm(cross(eB_z_des,eB_x_des));
             
    % with eB_x,y,z_des, reconstruct q_des 
    R_des = [eB_x_des eB_y_des eB_z_des];
    q_des = rotm2quat(R_des);
        
    q_e_y = quatmultiply(quatconj(quatmultiply(q,q_e_rp)),q_des);

    w_r_des = sign(q_e_y(1)) * 2 * k_p_yaw * q_e_y(4);
end
    
w_des = [w_p_des; w_q_des; w_r_des];
M_des = J*k_p_att*(w_des - w) + cross(w, J*w);

trpy = [0,0,0,0];
drpy = [0, 0,       0,         0];

end