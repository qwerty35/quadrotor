function[F_des, M_des, trpy, drpy] = Lee_controller(qd, t, qn, params)

% Extracing parameters and references
% --------------------
% Dynamics of quadrotor suspended with load Constants
m = params.mass;
J = params.I;

g = params.grav;

r_ref = qd{qn}.pos_des;
v_ref = qd{qn}.vel_des;
a_ref = qd{qn}.acc_des;

w_ref = [0; 0; 0];
M_ref = [0; 0; 0];

e1 = [1; 0; 0];
e2 = [0; 1; 0];
e3 = [0; 0; 1];

% Extracing states
% ----------------
r = qd{qn}.pos;
v = qd{qn}.vel;
R = RPYtoRot_ZXY2(qd{qn}.euler(1),qd{qn}.euler(2),qd{qn}.euler(3));
Omega = qd{qn}.omega;

Omegad = w_ref;
dOmegad = M_ref;

% CONTROL
% ------
% Position Control
eQ = r - r_ref;
deQ = v - v_ref;

k1 =  0.12*diag([4, 4 ,9.8055*1.2]);
k2 = 0.5*diag([4, 4, 10]);

A = -k1*eQ - k2*deQ + m*(a_ref + g*e3);
b3c = A/norm(A);
b3 = R(:,3);

F_des = vec_dot(A,b3);
    
% Attitude Control
b1d = e1;
b1c = -vec_cross(b3c,vec_cross(b3c,b1d));
b1c = b1c/norm(vec_cross(b3c,b1d));
Rc = [b1c vec_cross(b3c,b1c) b3c];
Rd = Rc;

if(norm(Rd'*Rd-eye(3)) > 1e-2)
    disp('Error in R') ; keyboard ;
end

kR = 15 ; kOm = 4;

err_R = 1/2 * vee_map(Rd'*R - R'*Rd) ;
err_Om = Omega - R'*Rd*Omegad ;
M_des = -kR*err_R - kOm*err_Om + vec_cross(Omega, J*Omega)...
    - J*(hat_map(Omega)*R'*Rd*Omegad - R'*Rd*dOmegad) ;

trpy = [0,0,0,0];
drpy = [0, 0,       0,         0];

end

