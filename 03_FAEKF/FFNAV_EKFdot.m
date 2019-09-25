function [f] = FFNAV_EKFdot(state, mu)
% FFNAV EKF Dot ===========================================================
% Description: This function defines the nonlinear equations of relative
% motion for formation flying spacecraft. When passed a state vector, this
% function calculates the respective accelerations for the x, y, z, rt and
% theta functions, and returns them to the calling function in a vector.
%
% Inputs:
%   state - The current state vector
%   mu    - Earth's gravitational parameter
%
% Outputs:
%   f       - Diffential equation output accelerations
%
% Created by: Cory Fraser - FALL... 2015
% Last Edits: Cory Fraser - JUL 04, 2018
% Copyright(c) 2018 by Cory Fraser
% =========================================================================

%% State Vector
 x          = state(1);
 y          = state(2);
 z          = state(3);
 theta      = state(4);
 rt         = state(5);
 x_dot      = state(6);
 y_dot      = state(7);
 z_dot      = state(8);
 theta_dot  = state(9);
 rt_dot     = state(10);
 
  rc        = sqrt((rt + x)^2 + y^2 + z^2);
  
%% Derivatives of the state variables
 %x_dot     = x_dot_priori; 
 %y_dot     = y_dot_priori;
 %z_dot     = z_dot_priori;
 %theta_dot = theta_dot_priori;
 %rt_dot    = rt_dot_priori;
 x_ddot     = x*theta_dot^2 + 2*theta_dot*(y_dot-y*rt_dot/rt)+mu/rt^2-mu*(rt+x)/(rc^3);
 y_ddot     = y*theta_dot^2 - 2*theta_dot*(x_dot-x*rt_dot/rt)-mu*y/(rc^3);
 z_ddot     = -mu*z/rc^3;
 theta_ddot = -2*rt_dot*theta_dot/rt;
 rt_ddot    = rt*theta_dot^2-mu/rt^2;
   
%% Assembling the derivative of the state vector 
f = [ x_dot; y_dot; z_dot; theta_dot; rt_dot; x_ddot; y_ddot; z_ddot; theta_ddot; rt_ddot ];

end
% =========================================================================