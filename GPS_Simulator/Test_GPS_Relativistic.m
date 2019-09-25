%% GPS Simulation - Relativistic Clock Error Test =========================
% =========================================================================
% Description: This 
%
% Inputs:
%   ECC
%   Ek
%   SMA
%
% Outputs:
%   R_rel - Range error due to ionospheric interference [m]
%   T_rel - Time error due to ionospheric interference  [s]
%   figure - Delay variation throughout 24-hour local time
%
% Reference: Montenbruck & Gill, Page 212
%            Klobuchar, Page 329
%
% Created by: Cory Fraser - Feb 01, 2018
% Last Edit : Cory Fraser - Feb 01, 2018
%% ========================================================================
%clc; clear; close all;

ECC = 0.01;          % Approximate GPS orbit eccentricity
Ek  = 90*pi/180;     % Maximum point of relativistic effects
SMA = 26561.75e3;    % Approximate SMA of GPS Orbit [m]

c   = 2.99792458e8;  % Speed of light [m/s]
mu  = 3.986005e14;   % Earth's gravitational parameter [m^3/s^2] 

[R_rel, t_rel] = Error_GPS_Relativistic(ECC,Ek,SMA, mu, c)