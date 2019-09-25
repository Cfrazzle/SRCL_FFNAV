function [f] = FFNAV_EKFdot(state_k)
% Description: This function defines the nonlinear equations of relative
% motion for formation flying spacecraft. When passed a state vector, this
% function calculates the respective accelerations for the x, y, z, rt and
% theta functions, and returns them to the calling function in a vector.
%
% Inputs:
%   The priori state vector

% Outputs:
%   The calculate diffential equation outputs
%
% Created by: Cory Fraser - Fall 2015
% =========================================================================
mu              = 398600.4418;          %Earth's Gravitational Constant (km^3/s^2)

% State Vector
 x          = state_k(1);
 y          = state_k(2);
 z          = state_k(3);
 theta      = state_k(4);
 rt         = state_k(5);
 x_dot      = state_k(6);
 y_dot      = state_k(7);
 z_dot      = state_k(8);
 theta_dot  = state_k(9);
 rt_dot     = state_k(10);
 
  rc        = sqrt((rt + x)^2 + y^2 + z^2);
  
% Derivatives of the state variables
 %x_dot     = x_dot_priori; As above
 %y_dot     = y_dot_priori;
 %z_dot     = z_dot_priori;
 %theta_dot = theta_dot_priori;
 %rt_dot    = rt_dot_priori;
 x_ddot     = x*theta_dot^2 + 2*theta_dot*(y_dot-y*rt_dot/rt)+mu/rt^2-mu*(rt+x)/(rc^3);
 y_ddot     = y*theta_dot^2 - 2*theta_dot*(x_dot-x*rt_dot/rt)-mu*y/(rc^3);
 z_ddot     = -mu*z/rc^3;
 theta_ddot = -2*rt_dot*theta_dot/rt;
 rt_ddot    = rt*theta_dot^2-mu/rt^2;
    
    f = [ x_dot; y_dot; z_dot; theta_dot; rt_dot; x_ddot; y_ddot; z_ddot; theta_ddot; rt_ddot ];
end
% =========================================================================