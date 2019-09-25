function [polar] = PolarAngles(R_vec)
% Polar Angles ================================================
% =========================================================================
% Description: This function calculates polar angles for a given position
% vector
%
% Inputs:
%   R - Position vector [m]
%
% Outputs:
%   DEC  - Declination (delta) [rad]
%   RA   - Right Ascension (alpha) [rad]
%   r    - Magnitude of the position vector [m]
%
% Reference: Curtis, Chapter 4, 2013, pages 193-194
%  
% Created by: Cory Fraser - JUN 30, 2018
% Last Edit : Cory Fraser - JUN 30, 2018
% Copyright(c) 2018 by Cory Fraser
% =========================================================================

%Calculate the magnitude of the vector 
r = norm(R_vec);
x = R_vec(1);
y = R_vec(2);
z = R_vec(3);

% Calculate the direction cosines
l = x/r;
m = y/r;
n = z/r;

%Calculate the declination
DEC = asin(n);

%Calculate the Right Ascension
if x == 0 && y == 0
    RA = 0;
else
    RA = acos(l/cos(DEC));
end
if m <= 0
    RA = 2*pi - RA;
end
RA = AngleRange(RA*180/pi)*pi/180;

polar = [DEC,RA,r];

end