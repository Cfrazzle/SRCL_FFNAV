function [f] = FFNAV_HCW_EKFdot(state_k, n)
% Description: This function defines the Hill-Clohessy-Wiltshire equations of relative
% motion for formation flying spacecraft. When passed a state vector, this
% function calculates the respective accelerations for the x, y, and z
% functions, and returns them to the calling function in a vector.
%
% Inputs:
%   The priori state vector

% Outputs:
%   The calculated diffential equation outputs
%
% Created by: Cory Fraser - Fall 2015
% Latest Edit: Cory Fraser - July 10, 2017
% =========================================================================
mu              = 398600.4418;          %Earth's Gravitational Constant (km^3/s^2)

% State Vector
 x          = state_k(1);
 y          = state_k(2);
 z          = state_k(3);
 x_dot      = state_k(4);
 y_dot      = state_k(5);
 z_dot      = state_k(6);
  
% Derivatives of the state variables
 %x_dot     = x_dot_priori; As above
 %y_dot     = y_dot_priori;
 %z_dot     = z_dot_priori;
 x_ddot     = 3*x*n^2 + 2*n*y_dot;
 y_ddot     = -2*n*x_dot;
 z_ddot     = -z*n^2;   
 f = [ x_dot; y_dot; z_dot; x_ddot; y_ddot; z_ddot ];
end
% =========================================================================