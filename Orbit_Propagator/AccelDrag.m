function [r_ddot] = AccelDrag(rho,R_Sat,V_Sat,T,Ad,m,Cd,w_Earth)
%% Atmospheric Drag Model =================================================
% =========================================================================
% Description: This function calculates the accelerations acting on a 
% spacecraft due to atmospheric drag.
%
% Inputs:
%   rho     - Density of the atmoshpere [kg/m^3]
%   R_Sat   - Position vector of the spacecraft in ECI frame [m]
%   V_Sat   - Velocity vector of the spacecraft in ECI frame [m/s]
%   T       - Transformation to true-of-date inertial system
%   Ad      - Area exposed to drag forces [m^2]
%   m       - Spacecraft mass [kg]
%   Cd      - Coefficient of Drag
%   w_Earth - Earth's angular velocity [rad/s]
%
% Outputs:
%   r_ddot - Acceleration of the spacecraft [m/s^2]
%
% Other Function Calls:
%   None
%
% Reference: Montenbruck & Gill, Chapter 3, 2011, pages 83-104
%  
% Created by: Cory Fraser - JUN 25, 2018
% Last Edit : Cory Fraser - JUN 25, 2018
% Copyright(c) 2018 by Cory Fraser
% =========================================================================

% Angular velocity vector of Earth 
omega_Earth = [0; 0; w_Earth];

%Transform position and velocity to true-of-date inertial system
%R_Sat = T * R_Sat;
%V_Sat = T * V_Sat;

% Relative velocity vector
v_rel = V_Sat - cross(omega_Earth, R_Sat);
e_v = v_rel/norm(v_rel); %Unit Vector

r_ddot = -(1/2)*rho*norm(v_rel)^2*(Cd*Ad/m)*e_v;

%Transform acceleration back to inertial system system
%r_ddot = T' * r_ddot_tod;