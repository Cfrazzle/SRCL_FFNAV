%% FFNAV HCW EKF Parameters ===============================================
% =========================================================================
% Description: This script constructs the initial conditions for the 
% extended Kalman filtering algorith.
%
% Inputs:
%   a Priori State Estimates
%      x_k, y_k, z_k
%      x_dot_k, y_dot_k, z_dot_k
%
%   a Priori Covariance Estimate
%      P_k
%
% Measurements
%       x_m,  y_m, z_m, x_dot_m, y_dot_m, z_dot_m
%   
% Created by: Cory Fraser - July 2016
% Last Edit : Cory Fraser - July 10, 2017
%% ========================================================================
%General system constants
mu              = 398600.4418;          %Earth's Gravitational Constant (km^3/s^2)
R_Earth         = 6378;                 %Average Radius of the Earth
w_Earth         = 7.292115*10^-5;       %Mean angular velocity of the Earth (rad/sec)

%Definitions for RWOP - Initial COE for Target Spacecraft
theta0_target   = 0;                    %True Anomaly (rad)
a_target        = 7500;                 %Semi-major Axis (km)
e_target        = 0;                    %Eccentricity
i_target        = pi/6;                 %Inclination (rad)
RAAN_target     = 0;                    %Right-ascension of the Ascending Node (rad)
AoP_target      = 0;                    %Argument of Perigee (rad)
n               = sqrt(mu/a_target^3);  %Mean orbital motion of the target spacecraft

%Definitions for RWOP - Initial COE for Chaser Spacecraft
theta0_chaser   = theta0_target;        %True Anomaly (rad)
a_chaser        = a_target;             %Semi-major Axis (km)
e_chaser        = 0.000005*10^4;              %Eccentricity - delta e
i_chaser        = i_target;             %Inclination (rad)
RAAN_chaser     = RAAN_target;          %Right-ascension of the Ascending Node (rad)
AoP_chaser      = AoP_target;           %Argument of Perigee (rad)

%Additional parameters for initial conditions
p_target        = a_target*(1-e_target^2);                          %Semi-latus rectum  
h_target        = sqrt(mu*p_target);                                %Angular momentum
rt0_target      = p_target/(1+e_target*cos(theta0_target));         %Initial radius
theta_dot0_target = h_target/rt0_target^2;                          %Initial theta_dot 
rt_dot0_target  = sqrt(mu/p_target)*e_target*sin(theta0_target);    %Initial r_dot

%% ========================================================================
%%Initializing parameters for EKF propagation

% Known Initial Conditions for In-Plane Elliptical Formation
x0_chaser       = -0.075*0.5*10^4;               %Initial x-position of chaser(km)
y0_chaser       = 0;                   %Initial y-position of chaser(km)
z0_chaser       = 0;                    %Initial z-position of chaser(km)
theta0          = theta0_target;        %Initial true anomaly of target (rad)
rt0             = rt0_target;           %Initial radius of target (km)

x_dot0_chaser   = y0_chaser*n/2;        %Initial x-velocity of chaser(km/s)
y_dot0_chaser   = -2*n*x0_chaser;       %Initial y-velocity of chaser(km)
z_dot0_chaser   = 0;                    %Initial z-velocity of chaser(km)
theta_dot0      = theta_dot0_target;    %Initial true anomaly rate (rad/s)
rt_dot0         = rt_dot0_target;       %Initial change in radius (km/s)
                         
state0 =    [   x0_chaser            % x
                y0_chaser            % y
                z0_chaser            % z
                x_dot0_chaser        % x_dot
                y_dot0_chaser        % y_dot
                z_dot0_chaser ];     % z_dot


%Defining the initial state convariance matrix (6 x 6)            
cov0            =10^-3*eye(6);

%Arranging the initial covariance matrix into a vector for the EKF (1 x 36)
cov0            = [ cov0(1,:) cov0(2,:) cov0(3,:) cov0(4,:) cov0(5,:) cov0(6,:)]; 

%GPS accuracies - standard deviation values for Random Number Generator block
GPS_r_STD = (1.2*10^-3); %(km) - (-7) Equivalent of 1.2 m
GPS_v_STD = (0.03*10^-3); %(km/s) - (-8) Equivalent of 0.03 m/s
                
%Initial measurements - Assuming measured positions/velocities are correct
x0_m            = x0_chaser;   
y0_m            = y0_chaser;
z0_m            = z0_chaser;
x_dot0_m        = x_dot0_chaser;
y_dot0_m        = y_dot0_chaser;
z_dot0_m        = z_dot0_chaser;

%Arranging the initial measurement vector (6 x 1)
meas0 =     [   x0_m 
                y0_m 
                z0_m 
                x_dot0_m 
                y_dot0_m 
                z_dot0_m   ];

%Starting time for the simulation
t0 = 0;
%% ========================================================================
    sig2_r           = 2*10^-6; %Good - July 6   
    sig2_v           = 5*10^-9; % Better, velocity error std < 0.03
       
     Q0              = [ sig2_r          0 0 0 0 0 
                        0   sig2_r         0 0 0 0
                        0 0     sig2_r       0 0 0
                        0 0 0       sig2_v     0 0 
                        0 0 0 0         sig2_v   0 
                        0 0 0 0 0           sig2_v ];
                    
   sig2_rm           = 2*10^-4; %Good - July 7
   sig2_vm           = 5*10^-9; %Good - July 6
    
    R0              = [ sig2_rm                 0 0 0 0 0  
                        0   sig2_rm               0 0 0 0  
                        0 0     sig2_rm             0 0 0  
                        0 0 0       sig2_vm           0 0 
                        0 0 0 0         sig2_vm         0   
                        0 0 0 0 0           sig2_vm      ];                          
% =========================================================================
%Assembling initial vector to be passed to the EKF simulation (1x49) 
EKF0            = [ t0 state0' cov0 meas0'  ];
filename        = 'FFNAV_HCW_PEO_Params.mat';
save(filename)
fprintf('\n Parameter file %s was created. \n', filename)

%{
%% ========================================================================
% Parameters from Busse Experiment
% =========================================================================
%Definitions for RWOP - Initial COE for Target Spacecraft
theta0_target   = 135*pi/180; %0;                    %True Anomaly (rad)
a_target        = 7017.9957; %7500;                 %Semi-major Axis (km)
e_target        = 0.005; %0;                    %Eccentricity
i_target        = 28.5*pi/180; %pi/6;                 %Inclination (rad)
RAAN_target     = 0;                    %Right-ascension of the Ascending Node (rad)
AoP_target      = -88.76*pi/180; %0;                    %Argument of Perigee (rad)
n               = sqrt(mu/a_target^3);  %Mean orbital motion of the target spacecraft

%Definitions for RWOP - Initial COE for Chaser Spacecraft
theta0_chaser   = 136.13*pi/180; %theta0_target;        %True Anomaly (rad)
a_chaser        = 7017.9886; %a_target;             %Semi-major Axis (km)
e_chaser        = 0.005103; %0.00001;              %Eccentricity - delta e
i_chaser        = i_target;             %Inclination (rad)
RAAN_chaser     = -1.13*pi/180; %RAAN_target;          %Right-ascension of the Ascending Node (rad)
AoP_chaser      = AoP_target;           %Argument of Perigee (rad)

%Additional parameters for initial conditions
p_target        = a_target*(1-e_target^2);                          %Semi-latus rectum  
h_target        = sqrt(mu*p_target);                                %Angular momentum
rt0_target      = p_target/(1+e_target*cos(theta0_target));         %Initial radius
theta_dot0_target = h_target/rt0_target^2;                          %Initial theta_dot 
rt_dot0_target  = sqrt(mu/p_target)*e_target*sin(theta0_target);    %Initial r_dot

%%Initializing parameters for EKF propagation

% Known Initial Conditions for In-Plane Elliptical Formation
x0_chaser       = 0.832187577110917; %-0.075;               %Initial x-position of chaser(km)
y0_chaser       = 16.988800024893138; %0;                    %Initial y-position of chaser(km)
x_dot0_chaser   = -0.017771885771207; %y0_chaser*n/2;        %Initial x-velocity of chaser(km/s)
theta0          = 2.356194490192345;        %Initial true anomaly of target (rad)
rt0             = rt0_target;           %Initial radius of target (km)
y_dot0_chaser   = -0.001206063704430; %-2*n*x0_chaser;       %Initial y-velocity of chaser(km)
z0_chaser       = 45.312543744836360; %0;                    %Initial z-position of chaser(km)
z_dot0_chaser   = -0.051399538490091; %0;                    %Initial z-velocity of chaser(km)
theta_dot0      = theta_dot0_target;    %Initial true anomaly rate (rad/s)
rt_dot0         = rt_dot0_target;       %Initial change in radius (km/s)
                                             
state0 =    [   x0_chaser            % x
                y0_chaser            % y
                z0_chaser            % z
                theta0               % theta
                rt0                  % rt
                x_dot0_chaser        % x_dot
                y_dot0_chaser        % y_dot
                z_dot0_chaser        % z_dot
                theta_dot0           % theta_dot
                rt_dot0    ];        % rt_dot

%Defining the initial state convariance matrix (10 x 10)            
cov0            =10^-11*eye(10);

%Arranging the initial covariance matrix into a vector for the EKF (1 x 100)
cov0            = [ cov0(1,:) cov0(2,:) cov0(3,:) cov0(4,:) cov0(5,:) ...
                    cov0(6,:) cov0(7,:) cov0(8,:) cov0(9,:) cov0(10,:)]; 

%Initial measurements - Assuming measured positions/velocities are correct
x0_m            = x0_chaser;   
y0_m            = y0_chaser;
z0_m            = z0_chaser;
theta0_m        = theta0;
rt0_m           = rt0;
x_dot0_m        = x_dot0_chaser;
y_dot0_m        = y_dot0_chaser;
z_dot0_m        = z_dot0_chaser;

%Arranging the initial measurement vector (7 x 1)
meas0 =     [   x0_m 
                y0_m 
                z0_m 
                theta0_m
                rt0_m
                x_dot0_m 
                y_dot0_m 
                z_dot0_m   ];

%Starting time for the simulation
t0 = 0;

%GPS accuracies - variance values for Random Number Generator block
GPS_r_RMS = 1.2*10^-7; %(km)
GPS_v_RMS = 0.03*10^-8; %(km/s)

%Assembling initial vector to be passed to the EKF simulation (1x118) 
EKF0            = [ t0 state0' cov0 meas0'  ];
filename        = 'FFNAV_Busse_Params.mat';
save(filename)
fprintf('\n Parameter file %s was created. \n', filename)
%}
pause(3)
clear
clc
