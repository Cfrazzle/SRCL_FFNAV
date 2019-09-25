function [Phi, F] = FFNAV_EKF_STM(state_pre, T)
% Description: This function calculates the discrete-time state transition
% matrix, given the current state and the time step
%
% Inputs:
%   state_pre - The previous state vector

% Outputs:
%   Phi - State Transition Matrix
%   F   - Jacobian of the Dynamic Model    
%
% Created by: Cory Fraser - August 31, 2017
% Last Edits: Cory Fraser - September 2, 2017
% =========================================================================
mu              = 398600.4418;          %Earth's Gravitational Constant (km^3/s^2)

% State Vector
 x          = state_pre(1);
 y          = state_pre(2);
 z          = state_pre(3);
 theta      = state_pre(4);
 rt         = state_pre(5);
 x_dot      = state_pre(6);
 y_dot      = state_pre(7);
 z_dot      = state_pre(8);
 theta_dot  = state_pre(9);
 rt_dot     = state_pre(10);
 
  rc        = sqrt((rt + x)^2 + y^2 + z^2);
  
%Partial derviations of the 5 non-linear equations to establish F

    %Derivatives of x-acceleration equation
    dxddot_dx           = (mu/rc^3)*(3*((rt+x)/rc)^2 - 1) + theta_dot^2;
    dxddot_dy           = 3*mu*((rt+x)/(rc^5))*y - ...
                          2*theta_dot*rt_dot/rt;
    dxddot_dz           = 3*mu*((rt+x)/rc^5)*z;
    dxddot_dtheta       = 0;
    dxddot_drt          = 2*theta_dot*rt_dot*y/rt^2 - ...
                          2*mu/rt^3 + (mu/rc^3)*(3*((rt+x)/rc)^2 - 1);
    dxddot_dx_dot       = 0; 
    dxddot_dy_dot       = 2*theta_dot;
    dxddot_dz_dot       = 0;
    dxddot_dtheta_dot   = 2*(theta_dot*x + y_dot - ...
                          (rt_dot/rt)*y);
    dxddot_drt_dot      = -2*theta_dot*y/rt;

    %Derivatives of y-acceleration equation
    dyddot_dx           = 2*theta_dot*rt_dot/rt + ...
                          3*mu*((rt+x)/rc^5)*y;
    dyddot_dy           = (mu/rc^3)*(3*(y/rc)^2-1) + theta_dot^2;
    dyddot_dz           = 3*mu*(y/rc^5)*z;
    dyddot_dtheta       = 0;
    dyddot_drt          = 3*mu*((rt+x)/rc^5)*y - ...
                          2*theta_dot*rt_dot*x/rt^2;
    dyddot_dx_dot       = -2*theta_dot; 
    dyddot_dy_dot       = 0;
    dyddot_dz_dot       = 0;
    dyddot_dtheta_dot   = 2*(theta_dot*y - x_dot + ...
                          rt_dot*x/rt);
    dyddot_drt_dot      = 2*theta_dot*x/rt;

    %Derivatives of z-acceleration equation
    dzddot_dx           = 3*mu*((rt+x)/rc^5)*z;
    dzddot_dy           = 3*mu*y*z/rc^5;
    dzddot_dz           = (mu/rc^3)*(3*(z/rc)^2-1);
    dzddot_dtheta       = 0;
    dzddot_drt          = 3*mu*((rt+x)/rc^5)*z;
    dzddot_dx_dot       = 0; 
    dzddot_dy_dot       = 0;
    dzddot_dz_dot       = 0;
    dzddot_dtheta_dot   = 0;
    dzddot_drt_dot      = 0;

    %Derivatives of theta-acceleration equation
    dtheta_ddot_dx          = 0;
    dtheta_ddot_dy          = 0;
    dtheta_ddot_dz          = 0;
    dtheta_ddot_dtheta      = 0;
    dtheta_ddot_drt         = 2*rt_dot*theta_dot/rt^2;
    dtheta_ddot_dx_dot      = 0; 
    dtheta_ddot_dy_dot      = 0;
    dtheta_ddot_dz_dot      = 0;
    dtheta_ddot_dtheta_dot  = -2*rt_dot/rt;
    dtheta_ddot_drt_dot     = -2*theta_dot/rt;
    
    %Derivatives of rt-acceleration equation
    drt_ddot_dx             = 0;
    drt_ddot_dy             = 0;
    drt_ddot_dz             = 0;
    drt_ddot_dtheta         = 0;
    drt_ddot_drt            = theta_dot^2 + 2*mu/rt^3;
    drt_ddot_dx_dot         = 0; 
    drt_ddot_dy_dot         = 0;
    drt_ddot_dz_dot         = 0;
    drt_ddot_dtheta_dot     = 2*theta_dot*rt;
    drt_ddot_drt_dot        = 0;

    %Defining the state matrix Ak
    F11 = zeros(5,5);
    F12 = eye(5,5);

    F21 = [dxddot_dx dxddot_dy dxddot_dz dxddot_dtheta dxddot_drt]; 
    F22 = [dxddot_dx_dot dxddot_dy_dot dxddot_dz_dot dxddot_dtheta_dot dxddot_drt_dot];

    F31 = [dyddot_dx dyddot_dy dyddot_dz dyddot_dtheta dyddot_drt];
    F32 = [dyddot_dx_dot dyddot_dy_dot dyddot_dz_dot dyddot_dtheta_dot dyddot_drt_dot];

    F41 = [dzddot_dx dzddot_dy dzddot_dz dzddot_dtheta dzddot_drt];
    F42 = [dzddot_dx_dot dzddot_dy_dot dzddot_dz_dot dzddot_dtheta_dot dzddot_drt_dot];

    F51 = [dtheta_ddot_dx dtheta_ddot_dy dtheta_ddot_dz dtheta_ddot_dtheta dtheta_ddot_drt];
    F52 = [dtheta_ddot_dx_dot dtheta_ddot_dy_dot dtheta_ddot_dz_dot ...
           dtheta_ddot_dtheta_dot dtheta_ddot_drt_dot];
    
    F61 = [drt_ddot_dx drt_ddot_dy drt_ddot_dz drt_ddot_dtheta drt_ddot_drt];
    F62 = [drt_ddot_dx_dot drt_ddot_dy_dot drt_ddot_dz_dot drt_ddot_dtheta_dot ...
           drt_ddot_drt_dot];

%Assembling the discrete state transition matrix Fk (10 x 10)
    F = [ F11 F12
            F21 F22
            F31 F32
            F41 F42
            F51 F52
            F61 F62 ];
   
    %Phi = expm(A_k*T); %Use the below Numerical Approximation;
    Phi = (F*T)*(((F*T)/2)*(((F*T)/3)*((F*T)/4 + eye(10,10))+eye(10,10))...
          +eye(10,10))+eye(10,10);     
end
% =========================================================================