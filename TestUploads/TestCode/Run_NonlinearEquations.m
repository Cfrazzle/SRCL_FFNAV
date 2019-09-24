%% Non-linear Equations of Relative Motion ================================
% =========================================================================
% Description: This program uses the user-defined LVLH position and 
% velocity components of the chaser spacecraft and the user-defined COE 
% of the target spacecraft to simulate the relative motion in LVLH using 
% the exact non-linear dynamic equations.
%
% Created by: Cory Fraser - Fall 2015
%% ========================================================================
% Define Initial COE for Target Spacecraft
theta0_target   = 0;          %True Anomaly
a_target    	= 7500;       %Semi-major Axis
e_target        = 0.0;        %Eccentricity
i_target        = pi/6;       %Inclination
RAAN_target     = 0;          %Right-ascension of the Ascending Node
AoP_target      = 0;          %Argument of Perigee
% =========================================================================
%% Define Non-linear Equation Initial Conditions (Units in km)
mu      = 398600;               %Earth's gravitational constant (km^3/s^2)
R_Earth = 6371;                 %Average radius of the Earth
w_Earth = 7.292115*10^-5;       %Mean angular velocity of the Earth
n       = sqrt(mu/a_target^3);  %Mean orbital motion of the target spacecraft

% In-Plane Elliptical Formation
x0_chaser       = 75;               %Initial x-position of chaser(km)
x_dot0_chaser   = .005;             %Initial x-velocity of chaser(km/s)
y0_chaser       = 2*x_dot0_chaser/n;%Initial y-position of chaser(km)
y_dot0_chaser   = -2*n*x0_chaser;   %Initial y-velocity of chaser(km)
z0_chaser       = 0;                %Initial z-position of chaser(km)
z_dot0_chaser   = 0;                %Initial z-velocity of chaser(km)
%}

%Additional Formations
%{
% Along-Track Formation
x0_chaser = 0;
x_dot0_chaser = 0;
y0_chaser = 0.2;
y_dot0_chaser = 0;
z0_chaser = 0; 
z_dot0_chaser = 0;

% In-Track Formation
x0_chaser = 0;
x_dot0_chaser = 0;
y0_chaser = 0.2;
y_dot0_chaser = 0;
z0_chaser = -1*(w_Earth/n)*sin(i_target)*y0_chaser; 
z_dot0_chaser = 0;

% Circular Formation
x0_chaser = 0.2;
x_dot0_chaser = 0.005;
y0_chaser = 2*x_dot0_chaser/n;
y_dot0_chaser = -2*n*x0_chaser;
z0_chaser = sqrt(3)*x0_chaser; 
z_dot0_chaser = sqrt(3)*x_dot0_chaser;

% Projected Circular Formation
x0_chaser = .2;
x_dot0_chaser = 0.005;
y0_chaser = 2*x_dot0_chaser/n;
y_dot0_chaser = -2*n*x0_chaser;
z0_chaser = 2*x0_chaser; 
z_dot0_chaser = 2*x_dot0_chaser;
%}

% Calculate additional initial conditions
p_target = a_target*(1-e_target^2);                 %Semi-latus rectum  
h_target = sqrt(mu*p_target);                       %Angular momentum of target orbit
rt0 = p_target/(1+e_target*cos(theta0_target));     %Initial radius of target
theta_dot0 = h_target/rt0^2;                        %Initial theta_dot of target 
rt_dot0 = sqrt(mu/(a_target*(1-e_target^2)))*e_target*sin(theta0_target); %Initial R_dot of target
% =========================================================================
%% Run Simulink Simulation
sim('Nonlinear_Relative_Motion_Equations');
% =========================================================================
%% Plotting the Results
figure(1)
plot3(x_chaser,y_chaser,z_chaser);
%axis([-2000 2000 -2000 2000 -2000 2000])
title('Relative Position (LVLH) - Chaser to the Target')
xlabel('X Position (km)')
ylabel('Y Position (km)')
zlabel('Z Position (km)')
set(gcf, 'position', [0 540 960 455])
grid on
% =========================================================================