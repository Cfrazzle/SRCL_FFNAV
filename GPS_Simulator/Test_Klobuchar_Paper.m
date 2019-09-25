%% GPS Simulation - Klobuchar Ionospheric Error Model =====================
% =========================================================================
% Description: This 
%
% Inputs:
%   r_Rec - Receiver position vector in ECEF frame      [m]
%   r_GPS - GPS satellite position matrix in ECEF frame [m]
%
% Outputs:
%   R_iono - Range error due to ionospheric interference [m]
%   T_iono - Time error due to ionospheric interference  [s]
%
% Reference: Montenbruck & Gill, Page 212
%            Klobuchar, Page 329
%
% Created by: Cory Fraser - January 29, 2018
% Last Edit : Cory Fraser - January 29, 2018
%% ========================================================================
clc; clear; close all
format long g
c        =  2.99792458e8;             % speed of light

%Coefficients 
alpha = [ 3.82e-8 
          1.49e-8
         -1.79e-7
         -0   ];

beta = [  1.43e5
         0
         -3.28e5
          1.13e5];

E = 20*pi/180;
%A = -100*pi/180;
A = 210*pi/180;

lat = 40*pi/180 /pi;
%lon = -116*pi/180 /pi;
lon = 260*pi/180 /pi;

t_GPS = 80720;
t_GPS = 593100;

% ========================================================================
% Klobuchar Model
% ========================================================================

% Compute the Earth-centered angle [semi-circles]
psi = (0.0137/(E/pi+0.11) - 0.022)';
%psi = (445/(E*180/pi + 20) - 4)'/180; %Degrees

% Determine the IPP sub-ionospheric latitude [semi-circles]
phi_L = lat + psi.*cos(A);
if phi_L > 0.416
    phi_L = 0.416;
elseif phi_L < -0.416
    phi_L = -0.416;
end

% Determine the sub-ionospheric longitude [semi-circles]
lambda_L = lon + psi.*sin(A)./cos(phi_L*pi);
if isnan(lambda_L)
    fprintf('\n NAN Detected \n')
end
%lambda_L = -0.6399;      %Conflicting "Correct" values from Klobuchar
%lambda_L = -124.795/180; = -0.6933 semi-circles...

% Calculate the IPP geomagnetic latitude [semi-circles]
phi_M =  phi_L + 0.064*cos((lambda_L - 1.617)*pi);
%phi_M = 0.2509; %"Correct" from Klobuchar

%Find the local time of day [s]              
t = 43200*lambda_L + t_GPS;
t = mod(t, 86400);
if t > 86400
    t = t - 86400;
elseif t < 0
    t = t + 86400;
end

% Compute the slant factor
F = 1.0 + 16.0*(0.53-E/pi).^3

% Period of the vTEC cosine approximation [s]
PER = beta(1) + beta(2).*phi_M + beta(3).*phi_M.^2 + beta(4).*phi_M.^3;
if PER < 72000
    PER = 72000;
end

% Phase of the vTEC cosine approximation [rad]
x = 2*pi*(t-50400)./PER;

% Amplitude of the vTEC cosine approximation [s]
AMP = alpha(1) + alpha(2).*phi_M + alpha(3).*phi_M.^2 + alpha(4).*phi_M.^3;
if AMP < 0
    AMP = 0;
end

% Calculate the signal delay
if abs(x) > 1.57
    T_iono = F*(5*10^-9);
else
    T_iono = F.*(5*10^-9 + AMP.*(1 - x.^2/2 + x.^4/24));
end

angles = [psi phi_L lambda_L phi_M]'*180

% =========================================================================
% Convert time delay to range error

R_iono = c*T_iono

% =========================================================================