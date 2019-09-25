%% Orbital Mechanics - GEO to ENZ =========================================
% =========================================================================
% Description: This function transformas a vector to the East-North-Zenith 
% (ENZ) reference frame, given the geodetic latitude, longitude and 
% components.
%
% Inputs:
%   r_ECEF - Position vector in ECEF frame [m]
%   gLLA   - Geodetic latitude, longitude and altitude [rad, rad, m]
%
% Outputs:
%   lat - Geodetic latitude (phi)
%   lon - Geodetic longitude (lambda)
%   alt - Altitude above the referecnce ellipsoid (WGS84)
%
% Reference: Montenbruck & Gill, Page 211
%
% Created by: Cory Fraser - January 29, 2018
% Last Edit : Cory Fraser - January 29, 2018
%% ========================================================================
function r_ENZ = gLLA_to_ENZ(r_ECEF, gLLA)

lat = gLLA(1);
lon = gLLA(2);

% Unit Vectors defining East-North-Zenith directions
e_East =  [ -sin(lon)    
             cos(lon)    
                0       ];
            
e_North = [ -sin(lat)*cos(lon)
            -sin(lat)*sin(lon)
                cos(lat)        ];
            
e_Zenith = [  cos(lat)*cos(lon)
              cos(lat)*sin(lon)
                sin(lat)        ];

%Rotation Matrix from ECEF to ENZ Frame
E = [e_East'
     e_North'
     e_Zenith'];
   
%Rotating the vector to ENZ
r_ENZ = E*r_ECEF';    

%Re-arranging back into column format
r_ENZ = r_ENZ';

end