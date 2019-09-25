%% GPS Simulation - Satellite Clock Offset Test ===========================
% =========================================================================
% Description: This 
%
% Inputs:
%   t_ref   - Satellite clock reference epoch [s]
%   t_GPS   - GPS system time at transmission [s]
%   alpha   - Clock correction coefficient matrix 
%              1) Clock bias        [s]
%              2) Frequency bias    [s/s]
%              3) Frequency drift   [s/s^2]
%
% Outputs:
%   R_ClockOffset - Range error due to ionospheric interference [m]
%   T_ClockOffset - Time error due to ionospheric interference  [s]

% Created by: Cory Fraser - Feb 01, 2018
% Last Edit : Cory Fraser - Feb 01, 2018
%% ========================================================================
clc; clear; close all;

%Correction Coefficient Set from Vigneron
af =    [ 3.76e-9   -1.60e-12   2.31e-17         %For PRN1
          7.28e-10  -1.76e-13   1.11e-17 ];      %For PRN2
      
t_GPS = [3000 5000];
t_ref = [75    74 ];

% Calculate the range error [m]
[R_ClockOffset, t_ClockOffset] = Error_GPS_Clock_Offset(af, t_GPS, t_ref)

% ========================================================================