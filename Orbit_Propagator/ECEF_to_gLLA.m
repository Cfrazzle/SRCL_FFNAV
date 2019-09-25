function geoLLA = ECEF_to_gLLA(r_ECEF, Constants)
%% Orbital Mechanics - ECEF to GEO ========================================
% =========================================================================
% Description: This function obtains the geodetic latitude, longitude and 
% altitude coordinates of a spacecraft, from its Earth-Centered 
% Earth-Fixed (ECEF) coordinates.
%
% Inputs:
%   r_ECEF - Position vector in ECEF frame [m]
%
% Outputs:
%   lat - Geodetic latitude (phi)    [rad]
%   lon - Geodetic longitude (lambda) [rad]
%   alt - Altitude above the referecnce ellipsoid [m] (WGS84)
%
% Reference:
%   Montenbruck - Page 187
%
% Created by: Cory Fraser - JAN 29, 2018
% Last Edit : Cory Fraser - JUN 30, 2018
% Copyright(c) 2018 by Cory Fraser
%% ========================================================================

%r_ECEF = [7000 7000 7000]*1000;

%ECEF Components
x = r_ECEF(1);  y = r_ECEF(2);  z = r_ECEF(3);

%Constants of WGS84 Reference Model
%load(Constants)
R_Earth = Constants.R_Earth;    %Radius of Earth [m]
f       = Constants.f_Earth;    %Flattening Factor (WGS84 Model, NIMA '97) 
ecc     = sqrt(1 -(1-f)^2);     %Ellipsoid Eccentricity

%R_Earth = 6378137;              %Radius of Earth [m]
%f       = 298.257223563^-1;     %Flattening Factor (WGS84 Model, NIMA '97) 
%ecc     = sqrt(1 -(1-f)^2);     %Ellipsoid Eccentricity

delta_z = ecc^2*z;

for i = 1:10
    sin_lat = (z + delta_z)/ sqrt(x^2 + y^2 + (z+ delta_z)^2);
    N       = R_Earth / sqrt(1 - ecc^2 * sin_lat^2);
    delta_z = N*ecc^2 *sin_lat;
end

lon = atan2(y,x);
lat = atan2( z + delta_z, sqrt(x^2 + y^2));
h   = sqrt(x^2 + y^2 + (z+delta_z)^2) - N;

b   = sqrt(R_Earth^2*(1-ecc^2));
k=abs(x)<1 & abs(y)<1;
h(k) = abs(z(k))-b;

geoLLA = [lat,lon,h];

end