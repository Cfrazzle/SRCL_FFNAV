function nu = AccelSolRadCyl(r_Sat,r_Sun, R_Earth)    
%% Solar Radiation Pressure ===============================================
% Description: This function calculates the fraction of the spacecraft in
% the shadow of the Earth, assuming that the shadow created by the Earth 
% extends as a cylinder.
%
% Inputs:
%   r_Sat       - Spacecraft position vector (ECI) [m] 
%   r_Sun       - Sun position vector (ECI) [m]
%   R_Earth     - Radius of the Earth [m]

% Outputs:
%   nu - Illumination Factor (1 = full illumination, 0 = eclipse)
%
% References:
%   Montenbruck, 'Satellite Orbits', pages 80-81, 2001
%
% Created by:  Cory Fraser - JUN 23, 2018
% Latest Edit: Cory Fraser - JUL 25, 2018
% Copyright(c) 2018 by Cory Fraser
% =========================================================================

%Unit vector pointing from the Earth to the Sun
e_Sun = r_Sun/norm(r_Sun);

%Scalar & vector projections of the spacecraft position on the Sun vector
s_proj = dot(r_Sat, e_Sun);
v_proj = s_proj*e_Sun;

%Illuminated Case 1: Spacecraft is on the Sun-side of Earth
Case1 = s_proj > 0;

%Illuminated Case 2: Spacecraft is outside Earth's cylindrical shadow  
Case2 = norm(r_Sat-v_proj) > R_Earth;

if ( Case1 || Case2)
    nu = 1;
else
    nu = 0;
end

%{
%Example from Curtis
rsat = norm(r_Sat);
rsun = norm(r_Sun);
%...Angle between sun and satellite position vectors:
theta = acosd(dot(r_Sat, r_Sun)/rsat/rsun);
%...Angle between the satellite position vector and the radial to the point
% of tangency with the earth of a line from the satellite:
theta_sat = acosd(R_Earth/rsat);
%...Angle between the sun position vector and the radial to the point
% of tangency with the earth of a line from the sun:
theta_sun = acosd(R_Earth/rsun);
%...Determine whether a line from the sun to the satellite
% intersects the earth:
if theta_sat + theta_sun <= theta
light_switch = 0; %yes
else
light_switch = 1; %no
end
%}

end