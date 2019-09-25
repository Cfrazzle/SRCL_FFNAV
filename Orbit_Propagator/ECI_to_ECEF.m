function [r_ECEF] = ECI_to_ECEF(r_ECI, t)
%% Orbital Mechanics - ECI to ECEF ========================================
% =========================================================================
% Description: This function coverts an Earth-Centered Inertia (ECI) 
% position vector into an Earth-Centered Earth-Fixed (ECEF) frame.
%
% Inputs:
%   
% Created by: Cory Fraser - January 29, 2018
% Last Edit : Cory Fraser - January 29, 2018
%% ========================================================================

t0 = 0;
omega_E = (2*pi + 2*pi/365.26)/(24*3600);

theta = omega_E*(t - t0);

R_ECEF_ECI =  [ cos(theta)   sin(theta)  0
                -sin(theta)   cos(theta)  0
                    0             0       1   ];
r_ECEF      = R_ECEF_ECI*r_ECI;

end

