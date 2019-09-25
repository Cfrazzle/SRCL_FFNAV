%% GPS Simulation - Klobuchar Ionospheric Model Test ======================
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
%   figure - Delay variation throughout 24-hour local time
%
% Reference: Montenbruck & Gill, Page 212
%            Klobuchar, Page 329
%
% Created by: Cory Fraser - Jan 29, 2018
% Last Edit : Cory Fraser - Feb 01, 2018
%% ========================================================================
%clc; clear; close all;

%Set of GPS Positions from Ephemerides (ECEF)
r_GPS = [   16126524.5052064, -15547762.6633447, 14383520.4858258;
            12603610.5415753,  12117327.4879380, 20031907.9625531;
            25942474.2728668, -4759620.41858934, 4338909.54401078;
            21058564.0351066,  16301904.2272398, 2284049.46872917;
            10072743.8616992,  21064317.2759496, 13011326.1051429;
            15934080.4069084, -4819743.53392942, 20534394.5855959];

r_Rec = [   3894192,  3894192, 0  ]; %ECEF
%r_Rec = [   0 0 -6400*1000  ]; %ECEF

%t_GPS = 57600;
t_GPS = -3*60*60;

%Ionosphere Cooefficients from GPS Message (Test Code)
%{
alpha = [ 1.21071934700012e-08 
         -7.45058059692383e-09
         -1.19209289550781e-07
          5.96046447753906e-08   ];

beta = [ 96256.0000000001
        -32768.0000000000
        -196608.000000000
         196608.000000000    ]; 
     %}

%Ionospheric coefficients from GPS - Yang et Al, 2014
%{
alpha = [ 0.6519e-8 
         -0.1490e-7
         -0.5960e-7
         -0.1192e-6   ];

beta = [  1.1525e5
         -1.7174e4
         -8.1250e5
          3.9824e6];
%}

%Ionospheric coefficients from Klobuchar
alpha = [ 3.82e-8 
          1.49e-8
         -1.79e-7
         -0   ];

beta = [  1.43e5
         0
         -3.28e5
          1.13e5];
%}

t_test = t_GPS:30:t_GPS+24*60*60;      
for i = 1:length(t_test)
    R_iono(:,i) = Error_Iono_Klobuchar(r_Rec, r_GPS, alpha, beta, t_test(i));
end

% =========================================================================
figure(1)
clf
plot((t_test-t_GPS)/60/60, R_iono(:,:))
axis([0 24 0 60]);
xlabel('t_G_P_S + time (hr)')
ylabel('Range Error [m]')
title('Klobuchar Ionospheric Delay Model');
set(gca,'xtick',[0 4 8 12 16 20 24],'ytick',[0 10 20 30 40 50 60] ,'Xgrid','off','Ygrid','on');
