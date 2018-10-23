function params = crazyflie()
% crazyflie: physical parameters for the Crazyflie 2.0
%
% 2016 Bernd Pfrommer
%
% This function creates a struct with the basic parameters for the
% Crazyflie 2.0 quad rotor (without camera, but with about 5 vicon
% markers)
%
% Model assumptions based on physical measurements:
%
% motor + mount + vicon marker = mass point of 3g
% arm length of mass point: 0.046m from center
% battery pack + main board are combined into cuboid (mass 18g) of
% dimensions:
%
%   width  = 0.03m
%   depth  = 0.03m
%   height = 0.012m
%

%%%%%cage무게도 고려 m, I모두

m = 0.030;  % weight (in kg) with 5 vicon markers (each is about 0.25g) 
g = 9.81;   % gravitational constant
I = [1.43e-5,   0,          0; % inertial tensor in m^2 kg
     0,         1.43e-5,    0;
     0,         0,          2.89e-5];
L_arm = 0.046; % arm length in m
L_leg = 0.017; % leg length in m
L_propeller = 0.022; % propeller length in m
H_propeller = 0.010; % propeller height in m

params.mass = m;
params.I    = I;
params.invI = inv(I);
params.grav = g;
params.arm_length = L_arm;
params.propeller = [L_leg; L_propeller; H_propeller];

params.maxangle = 40*pi/180; % you can specify the maximum commanded angle here
params.maxF     = 2.5*m*g;   % left these untouched from the nano plus
params.minF     = 0.05*m*g;  % left these untouched from the nano plus
params.k_F      = 6.11*1e-8; 
params.k_M      = 1.5*1e-9; 
params.k_m      = 20; 

% Collision parameters
% Collision model is based on below paper
% Dynamics of a Quadrotor Undergoing Impact with a wall
% Fiona et al
params.cage_type = 'sphere'; % cage type, 'sphere' or 'disc' 
params.R = 0.080; % cage radius (m)
params.e = 0.3;  % coeffcient of restitution
params.n = 0.54;  % collision type order
params.k_C = 40;  % constant stiffness coefficient (N/m)
params.u_C = 0.2; % coefficient of friction

% You can add any fields you want in params
% for example you can add your controller gains by
% params.k = 0, and they will be passed into controller.m

end
