%% Orbital Mechanics - Calculate Azimuth and Elevation ====================
% =========================================================================
% Description: This function calculates the azimuth and elevation angles
% between a spacecraft (GPS receiver) and a GPS satellite.
%
% Inputs:
%   r_Rec - Receiver position vector in ECEF frame [m]
%   r_GPS - GPS satellite position vector in ECEF frame [m]
%
% Outputs:
%   azi - Azimuth angle    (A)  [rad]
%   ele - Elevation angle  (E)  [rad]
%
% Reference: Montenbruck & Gill, Page 212
%
% Created by: Cory Fraser - January 29, 2018
% Last Edit : Cory Fraser - January 29, 2018
%% ========================================================================
function [azi, ele] = Calc_AzEl(r_Rec, r_GPS)

% ECEF vector from Reciever to GPS Satellite
R = r_GPS - r_Rec; 

% Calculating lat,long, altitude parameters of the reciever
Rec_gLLA = ECEF_to_gLLA(r_Rec);       

% Rotating position vector to ENZ Frame
R_ENZ = gLLA_to_ENZ(R,Rec_gLLA);

% Components of position unit vector in ENZ Frame
s_East   = R_ENZ(:,1);
s_North  = R_ENZ(:,2);
s_Zenith = R_ENZ(:,3);

% Calculating the elevation and azimuth
azi = atan2(s_East,s_North);
ele  = atan2(s_Zenith, sqrt(s_East.^2+s_North.^2));

end