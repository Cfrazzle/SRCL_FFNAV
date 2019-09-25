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
%   lat - Geodetic latitude (phi)
%   lon - Geodetic longitude (lambda)
%   alt - Altitude above the referecnce ellipsoid (WGS84)
%
% Reference: Montenbruck & Gill, Page 188
%            Bowring 1985 (Find Reference)
%
% Created by: Cory Fraser - JAN 29, 2018
% Last Edit : Cory Fraser - JAN 29, 2018
% Copyright(c) 2018 by Cory Fraser
%% ========================================================================
function gLLA = ECEF_to_gLLA(r_ECEF)


% ECEF Components
x = r_ECEF(1);  y = r_ECEF(2);  z = r_ECEF(3);

% Constants of WGS84 Reference Model
R_Earth = 6378137;                   %Radius of Earth [m]
f       = 298.257223563^-1;          %Flattening factor (WGS84 Model, NIMA '97) 
ecc     = sqrt(1 -(1-f)^2);          %Ellipsoid eccentricity
b       = sqrt(R_Earth^2*(1-ecc^2)); %Ellipsoid semi-major axis

% =========================================================================
% Check if position is at the origin
origin = (x == 0 & y == 0 & z == 0)
if origin
    lon = 0;
    lat = 0;
    alt = 0;
end

% =========================================================================
% Check if position is near the poles
poles = (abs(x) <= 1 & abs(y) <= 1 & z ~= 0)
if poles 
    lon = 0;
    lat = sign(z)*pi/2;
    alt = abs(z) - R_Earth*(1-f);
end

% =========================================================================
% Check if position is on the equator
equator = (~origin & ~poles & z ==0)
if equator
    lon = atan2(y,x)
    lat = 0
    alt = sqrt(x^2 + y^2) - R_Earth
end

% =========================================================================
% For all other cases
main = (~origin & ~poles & ~equator)
if main
    delta_z = ecc^2*z;
    
    for i = 1:10
        sin_lat = (z + delta_z)/ sqrt(x^2 + y^2 + (z+ delta_z)^2);
        N       = R_Earth / sqrt(1 - ecc^2 * sin_lat^2);
        delta_z = N*ecc^2 *sin_lat;
    end
    
    lon = atan2(y,x);
    lat = atan2( z + delta_z, sqrt(x^2 + y^2));
    alt = sqrt(x^2 + y^2 + (z+delta_z)^2) - N;
end

% =========================================================================
%Output the resulting geodetic components
gLLA = [lat,lon,alt];
end