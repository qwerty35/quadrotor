function sdot = quadEOM_Obstacle(s, F, M, params)
% QUADEOM_READONLY Solve quadrotor equation of motion
%   quadEOM_readonly calculate the derivative of the state vector
%
% INPUTS:
% t      - 1 x 1, time
% s      - 17 x 1, state vector = [x, y, z, xd, yd, zd, qw, qx, qy, qz, p, q, r, F1, F2, F3, F4]
% F      - 1 x 1, thrust output from controller (only used in simulation)
% M      - 3 x 1, moments output from controller (only used in simulation)
% params - struct, output from crazyflie() and whatever parameters you want to pass in
% map    - struct, map information to check collision
%
% OUTPUTS:
% sdot   - 17 x 1, derivative of state vector s
%
% NOTE: You should not modify this function
% See Also: quadEOM_readonly, crazyflie

%************ EQUATIONS OF MOTION ************************
% Assign states
x = s(1);
y = s(2);
z = s(3);
xdot = s(4);
ydot = s(5);
zdot = s(6);
qW = s(7);
qX = s(8);
qY = s(9);
qZ = s(10);
p = s(11);
q = s(12);
r = s(13);
F1 = s(14);
F2 = s(15);
F3 = s(16);
F4 = s(17);

% Limit the force and moments due to actuator limits
A = [0.25,                      0, -0.5/params.arm_length,  0.25*params.k_F/params.k_M;
     0.25,  0.5/params.arm_length,                      0, -0.25*params.k_F/params.k_M;
     0.25,                      0,  0.5/params.arm_length,  0.25*params.k_F/params.k_M;
     0.25, -0.5/params.arm_length,                      0, -0.25*params.k_F/params.k_M];

prop_thrusts = A*[F;M];
% prop_thrusts_clamped = max(min(prop_thrusts, params.maxF/4), params.minF/4);
prop_thrusts_temp = min(prop_thrusts, params.maxF/4);
prop_thrusts_clamped = prop_thrusts_temp .* (prop_thrusts_temp > params.minF/4);
curr_thrusts = [F1; F2; F3; F4];

B = [                    1,                      1,                     1,                      1;
                         0,      params.arm_length,                     0,     -params.arm_length;
        -params.arm_length,                      0,     params.arm_length,                      0;
     params.k_M/params.k_F, -params.k_M/params.k_F, params.k_M/params.k_F, -params.k_M/params.k_F];
 
F = B(1,:)*curr_thrusts;
M = B(2:4,:)*curr_thrusts;

quat = [qW; qX; qY; qZ];
bRw = QuatToRot(quat);
wRb = bRw';


% Acceleration
accel = 1 / params.mass * (wRb * [0; 0; F] - [0; 0; params.mass * params.grav]);

% Angular velocity
K_quat = 2; %this enforces the magnitude 1 constraint for the quaternion
quaterror = 1 - (qW^2 + qX^2 + qY^2 + qZ^2);
qdot = -1/2*[0, -p, -q, -r;...
             p,  0, -r,  q;...
             q,  r,  0, -p;...
             r, -q,  p,  0] * quat + K_quat*quaterror * quat;

% Angular acceleration
omega = [p;q;r];
pqrdot   = params.invI * (M - cross(omega, params.I*omega));

% Motor angular velocity
Fdot = params.k_m* 2 * (prop_thrusts_clamped - curr_thrusts);

% Assemble sdot
sdot = zeros(13,1);
sdot(1)  = xdot;
sdot(2)  = ydot;
sdot(3)  = zdot;
sdot(4)  = accel(1);
sdot(5)  = accel(2);
sdot(6)  = accel(3);
sdot(7)  = qdot(1);
sdot(8)  = qdot(2);
sdot(9)  = qdot(3);
sdot(10) = qdot(4);
sdot(11) = pqrdot(1);
sdot(12) = pqrdot(2);
sdot(13) = pqrdot(3);
sdot(14) = Fdot(1);
sdot(15) = Fdot(2);
sdot(16) = Fdot(3);
sdot(17) = Fdot(4);

end
