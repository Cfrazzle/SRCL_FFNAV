function [rho_atm] = Density_HarrisPriester(R_Sat, R_Sun, Constants)
% Atmospheric Density Model ================================================
% =========================================================================
% Description: This function calculates the density of the atmosphere using
% the Harris-Priester model
%
% Inputs:
%   R_Sat - Position vector of the spacecraft in ECI frame, in TOD
%           system [m]
%   r_Sun - Position vector of the Sun in the ECI frame, referred to the 
%           International Celestial Reference Frame (ICRF) [m]
%
% Outputs:
%   rho_atm - Density of the atmosphere [kg/m^3]
%
% Other Function Calls:
%   ECEF_to_gLLA - Determines satellite geodetic lat, long, and height 
%   CalcPolarAngles - Determines right ascension and declination angles
%
% Reference: Montenbruck & Gill, Chapter 3, 2011, pages 83-104
%  
% Created by: Cory Fraser - JUN 30, 2018
% Last Edit : Cory Fraser - JUN 30, 2018
% Copyright(c) 2018 by Cory Fraser
% =========================================================================

% Convert Satellite Position from ECI Frame to ECEF Frame 
t = 0;
R_Sat_ECEF = ECI_to_ECEF(R_Sat, t);
e_Sat = R_Sat_ECEF/norm(R_Sat_ECEF);

% Calculate satellite altitude [m]
geoLLA = ECEF_to_gLLA(R_Sat_ECEF, Constants);
alt = geoLLA(3);
alt = alt/1000;  %Convert to [km] for HP Model

% Check if the height is within HP Model Limits
if alt >= 1000 || alt <= 100
    %warning('Altitude of %.0f km is outside the limits of the HP Model', alt);
    rho_atm = 0;
    return
end
    
% Calculate the right ascension and declination of the Sun [rads]
[SunAngles] = PolarAngles(R_Sun);
DEC_Sun = SunAngles(1);
RA_Sun  = SunAngles(2);

% Calculate a unit vector to the the diurnal bulge
RA_lag      = 0.523599; % Right ascension lag ~ 30 degrees [rads]

e_bulge(1) = cos(DEC_Sun) * cos(RA_Sun + RA_lag);
e_bulge(2) = cos(DEC_Sun) * sin(RA_Sun + RA_lag);
e_bulge(3) = sin(DEC_Sun);

%% Extract coefficients for HP Model

%Coefficients for the Harriss-Priester Density Model
%Altitude [km]   %Min Density [g/km^3]  %Max Density [g/km^3]
HP_Coefficients = ...
    [  100      497400.0      497400.0
       120       24900.0       24900.0
       130        8377.0        8710.0
       140        3899.0        4059.0
       150        2122.0        2215.0
       160        1263.0        1344.0
       170         800.8         875.8
       180         528.3         601.0
       190         361.7         429.7
       200         255.7         316.2
       210         183.9         239.6
       220         134.1         185.3
       230          99.49        145.5
       240          74.88        115.7
       250          57.09        93.08
       260          44.03        75.55
       270          34.30         61.82
       280          26.97         50.95
       290          21.39         42.26
       300          17.08         35.26
       320          10.99         25.11
       340           7.214        18.19
       360           4.824        13.37
       380           3.274         9.955
       400           2.249         7.492
       420           1.558         5.684
       440           1.091         4.355
       460           0.7701        3.362
       480           0.5474        2.612
       500           0.3916        2.042
       520           0.2819        1.605
       540           0.2042        1.267
       560           0.1488        1.005
       580           0.1092        0.7997
       600           0.08070 	   0.6390
       620           0.06012       0.5123
       640           0.04519       0.4121
       660           0.03430       0.3325
       680           0.02632       0.2691
       700           0.02043       0.2185
       720           0.01607       0.1779
       740           0.01281       0.1452
       760           0.01036       0.1190
       780           0.008496      0.09776
       800           0.007069      0.08059
       840           0.004680      0.05741
       880           0.003200      0.04210   
       920           0.002210      0.03130
       960           0.001560      0.02360
       1000          0.001150      0.01810  ];

h         = HP_Coefficients(:,1);
CoEff_min = HP_Coefficients(:,2);
CoEff_max = HP_Coefficients(:,3); 

%Find the correct height index for interpolation
h_id = 0;
i = 1;
while ~h_id
    if alt >= h(i) && alt< h(i+1)
        h_id = i;
        %break
    end
    i = i+1;
end

% Calculate the minimum and maximum scale heights [km]
H_min = (h(h_id) - h(h_id+1)) / log( CoEff_min(h_id+1) / CoEff_min(h_id) ); 
H_max = (h(h_id) - h(h_id+1)) / log( CoEff_max(h_id+1) / CoEff_max(h_id) );

% Calculate the minimum and maximum density values [g/km^3]
rho_min = CoEff_min(h_id) * exp( (h(h_id) - alt)/H_min);
rho_max = CoEff_max(h_id) * exp( (h(h_id) - alt)/H_max);

% Calculate angles between the satellite and the apex of the diurnal bulge
%n = 2; %(Low inclination orbits)
n = 3; %(PROBA-3 Formation)
%n = 5; %(PRISMA Formation)
%n = 5; %(PEOinLEO Formation)
%n = 6; %(Polar orbits)

cos_npsi = ( 0.5 + 0.5*dot(e_Sat, e_bulge) )^(n/2); 

% Calculate the density, accounting for the dinual variations [g/km^3]
rho_atm = rho_min + ( rho_max - rho_min ) * cos_npsi;

% Convert to [kg/m^3]
rho_atm = rho_atm*1e-12; 

end