syms x      y       z       theta       rt      mu 
syms x_dot y_dot  z_dot  theta_dot  rt_dot
syms x_ddot y_ddot  z_ddot  theta_ddot  rt_ddot

state_vec = [ x y z theta rt x_dot  y_dot   z_dot   theta_dot   rt_dot].';

x_ddot = @(x, y, z, theta, rt, x_dot, y_dot, z_dot, theta_dot, rt_dot) theta_dot^2*x + 2*theta_dot*y_dot - 2*theta_dot*rt_dot*y/rt + mu/rt^2 - mu*(rt+x)*((rt+x)^2+y^2+z^2)^(-3/2);
y_ddot = @(x, y, z, theta, rt, x_dot, y_dot, z_dot, theta_dot, rt_dot) theta_dot^2*y - 2*theta_dot*x_dot + 2*theta_dot*rt_dot*x/rt - mu*y*((rt+x)^2+y^2+z^2)^(-3/2);
z_ddot = @(x, y, z, theta, rt, x_dot, y_dot, z_dot, theta_dot, rt_dot) -mu*z*((rt+x)^2+y^2+z^2)^(-3/2);
theta_ddot = @(x, y, z, theta, rt, x_dot, y_dot, z_dot, theta_dot, rt_dot) -2*rt_dot*theta_dot/rt;
rt_ddot = @(x, y, z, theta, rt, x_dot, y_dot, z_dot, theta_dot, rt_dot) theta_dot^2*rt - mu/rt^2;

part_ddot_x = jacobian(x_ddot, [x, y, z, theta, rt, x_dot, y_dot, z_dot, theta_dot, rt_dot]);
part_ddot_y = jacobian(y_ddot, [x, y, z, theta, rt, x_dot, y_dot, z_dot, theta_dot, rt_dot]);
part_ddot_z = jacobian(z_ddot, [x, y, z, theta, rt, x_dot, y_dot, z_dot, theta_dot, rt_dot]);
part_ddot_theta = jacobian(theta_ddot, [x, y, z, theta, rt, x_dot, y_dot, z_dot, theta_dot, rt_dot]);
part_ddot_rt = jacobian(rt_ddot, [x, y, z, theta, rt, x_dot, y_dot, z_dot, theta_dot, rt_dot]);
%

H = [   eye(4,4)    zeros(4,1) zeros(4,3) zeros(4,2) 
        zeros(3,4)  zeros(3,1) eye(3,3)   zeros(3,2) ];
    
f = [ x_dot  y_dot   z_dot   theta_dot   rt_dot x_ddot y_ddot  z_ddot  theta_ddot  rt_ddot].';

%==========================================================================
%Constructing the Observability Matrix
tic;

%Zeroth-Order Lie Derivative Gradient
G0 = H;
toc

%First-Order Lie Derivative Gradient
L1 = H*f;
G1 = jacobian(L1, state_vec);
toc
%Second-Order Lie Derivative Gradient
L2 = G1*f;
G2 = jacobian(L2, state_vec);
toc

%Third-Order Lie Derivative Gradient
L3 = G2*f;
G3 = jacobian(L3, state_vec);
toc
%
%Fourth-Order Lie Derivative Gradient
L4 = G3*f;
G4 = jacobian(L4, state_vec);
toc

%Fifth-Order Lie Derivative Gradient
L5 = G4*f;
G5 = jacobian(L5, state_vec);
toc

%Sixth-Order Lie Derivative Gradient
L6 = G5*f;
G6 = jacobian(L6, state_vec);
toc

%Seventh-Order Lie Derivative Gradient
L7 = G6*f;
G7 = jacobian(L7, state_vec);
toc

%Eighth-Order Lie Derivative Gradient
L8 = G7*f;
G8 = jacobian(L8, state_vec);
toc

%{
%Nineth-Order Lie Derivative Gradient
L9 = G8*f;
G9 = jacobian(L9, state_vec);
%}

OM = [G0; G1; G2; G3; G4; G5; G6; G7; G8];% G9];
OM_rank = rank(OM)
time = toc

filename        = 'Observability8_results.mat';
save(filename)
       
