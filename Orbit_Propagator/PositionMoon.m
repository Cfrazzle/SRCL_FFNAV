function [R_Moon] = PositionMoon(JD, R_Earth)
%% Moon Position ==========================================================
% Description: This function calculates the coordinates of the Moon in the 
% ECI reference frame, based on the Modified Julian Day.
%
% Inputs:
%   MJD_TT - Modified Julian Date, for the terrestrial time
%
% Outputs:
%   R_Moon - Position vector of the Moon [m] 
%
% Reference: Curtis, Orbital Mechanics for Engineering Students
%            3rd Edition, Pages 705-711
%
% Created by:  Cory Fraser - JUN 30, 2018
% Latest Edit: Cory Fraser - JUL 15, 2018
%==========================================================================

%Julian Day Number - CAN PASS FROM UPPER FUNCTION
%JD = JulianDate(25,07,2013,08,0,0); %From Example 12.10

if nargin <2
    R_Earth   = 6378.137e3;        %Equatorial Radius - WGS84 [m]
end

%Julian Day number to number of Centruies since J2000
T0 = (JD - 2451545.0)/36525;

%Obliquity of the J2000 Ecliptic [deg]
%n = JD - 2451545.0;            %Days since J2000
%epsilon = 23.439 - n*3.56e-7;  %Equivalent with n = 36525*T0
epsilon = 23.439 - 0.0130042*T0;

%Lunar ecliptic longitude [deg]
lambda_long  = 218.32 + 481267.881*T0           ...
             + 6.29*sind(135.0 + 477198.87*T0)  ...
             - 1.27*sind(259.3 - 413335.36*T0)  ...
             + 0.66*sind(235.7 + 890534.22*T0)  ...
             + 0.21*sind(269.9 + 954397.74*T0)  ...
             - 0.19*sind(357.5 + 35999.05*T0)   ...
             - 0.11*sind(186.5 + 966404.03*T0);
         
lambda_long = AngleRange(lambda_long);

%Lunar ecliptic latitude [deg]
delta_lat   = 5.13*sind( 93.3 + 483202.02*T0)   ...
            + 0.28*sind(228.2 + 960400.89*T0)   ...
            - 0.28*sind(318.3 + 6003.15*T0)     ...
            - 0.17*sind(217.6 - 407332.21*T0);
        
delta_lat = AngleRange(delta_lat);

%Lunar horizontal parallax 
h_par   = 0.9508 ...
        + 0.0518*cosd(135.0 + 477198.87*T0)     ...
        + 0.0095*cosd(259.3 - 413335.36*T0)     ...
        + 0.0078*cosd(235.7 + 890534.22*T0)     ...
        + 0.0028*cosd(269.9 + 954397.74*T0);
    
h_par = AngleRange(h_par);

%Distance of the Sun from the Earth [units match R_Earth] (Note: Angles in degrees)
r_Moon = R_Earth/sind(h_par);

%Unit vector from the Earth to the Moon (Note: Angles are in degrees)
u_x = cosd(delta_lat) * cosd(lambda_long);
u_y = cosd(epsilon)*cosd(delta_lat)*sind(lambda_long) - sind(epsilon)*sind(delta_lat);
u_z = sind(epsilon)*cosd(delta_lat)*sind(lambda_long) + cosd(epsilon)*sind(delta_lat);

u_Moon = [u_x; u_y; u_z];

%Position vector of the Moon in the ECI frame [m]
R_Moon = r_Moon*u_Moon;


end
