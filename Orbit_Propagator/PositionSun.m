function [R_Sun] = PositionSun(JD, AU)
%% Sun Position Vector ====================================================
% Description: This function calculates the coordinates of the Sun in the 
% ECI reference frame, based on the Modified Julian Day.
%
% Inputs:
%   MJD_TT - Modified Julian Date, for the terrestrial time
%
% Outputs:
%   R_Sun - Position vector of the Sun [m] 
%
% Reference: Curtis, Orbital Mechanics for Engineering Students
%            3rd Edition, Pages 698-701
%
% Created by:  Cory Fraser - JUN 30, 2018
% Latest Edit: Cory Fraser - JUL 15, 2018
%==========================================================================

%Julian Day Number - CAN PASS FROM UPPER FUNCTION
%JD = JulianDate(25,07,2013,08,0,0); %From Example 12.9

if nargin <2 
    AU = 149597870699.999988; %Astronomical Unit [m]
end

%Days since J2000
n = JD - 2451545.0;

%Mean solar longitude L and mean anomaly M [deg]
L = 280.459 + 0.98564736*n;  
M = 357.529 + 0.98560023*n;  

L = AngleRange(L); %Restrict 0 < L < 360 degrees
M = AngleRange(M); %Restrict 0 < M < 360 degrees

%Longitude of the Sun [deg] (Note: Angles are in degrees)
lambda = L + 1.915*sind(M) + 0.0200*sind(2*M); 
lambda = AngleRange(lambda); %Restrict 0 < lambda < 360 degrees

%Obliquity of the J2000 Ecliptic [deg]
epsilon = 23.439 - n*3.56e-7;

%Unit vector from the Earth to the Sun (Note: Angles are in degrees)
u_Sun = [cosd(lambda); sind(lambda)*cosd(epsilon); sind(lambda)*sind(epsilon)];

%Distance of the Sun from the Earth (Note: Angles are in degrees)
r_Sun = (1.00014 - 0.01671*cosd(M) - 0.000140*cosd(2*M))*AU; %[Units corresponds w/AU

%Position vector of the Sun in the ECI frame [m]
R_Sun = r_Sun*u_Sun;

end
