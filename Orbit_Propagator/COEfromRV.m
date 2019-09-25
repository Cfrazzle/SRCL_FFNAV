function [COE, Data] = COEfromRV(R, V, mu)
%% Position/Velocity Vectors to Classical Orbital Elements ================
% Description: This function uses the current position and velocity vectors
% of a spacecraft to determine the corresponding orbital elements.
%
% Inputs:
%   R_ECI - Position Vector (ECI Frame)
%   V_ECI - Velocity Vector (ECI Framce)
%   mu    - Gravitational Parameter [m^3/s^2]
%   
% Outputs:
%   COE Vector
%       SMA     - Semi-major Axis [m]
%       ECC     - Eccentricity    [unitless]  
%       INC     - Inclination     [degrees]
%       RAAN    - Right-ascension of the Ascending Node [degrees]
%       AoP     - Argument of Perigee [degrees]
%       TA      - True Amomaly [degrees]
%
%  DATA Vector
%       rp     - Perigee Radius [m]
%       ra     - Apogee Radius [m]
%        T     - Orbital Period [s]
%       EA     - Eccentric Anomaly [deg]
%
% Created by:  Cory Fraser - FEB 11, 2017
% Latest Edit: Cory Fraser - JUN 20, 2018
% ========================================================================

% Definition of ECI Frame
x_ECI   = [1 0 0]; 
y_ECI   = [0 1 0];
z_ECI   = [0 0 1];

%Calculating the psotion/velocity magnitudes
r = norm(R);  v = norm(V);

% Semi-major Axis
SMA = r / (2 - (r*v^2)/mu );

%% Eccentricity Vector and Magnitude
ECC_vec = (1/mu)*( (v^2 - mu/r)*R - dot(R,V)*V );
ECC = norm(ECC_vec); 

ECC_min = 1e-8; %Eccentricity limit for circular orbits

%Angular Momentum
h_vec = cross(R, V);
h = norm(h_vec);

%% Inclination
INC = acosd(h_vec(3)/h);

%Defining the Nodal Vector
n_vec = cross(z_ECI,h_vec/h);
n = norm(n_vec);

%% Calculating the Right-ascension of the Ascending Node
if n ~= 0                       %Inclined Orbit Case
    RAAN = acosd(n_vec(1)/n);
    if n_vec(2) < 0 
        RAAN = 360 - RAAN;
    end
else                            %Equatorial Orbit Case
    RAAN = 0;
end

%% Calculating the Argument of Perigee
if n ~= 0                       %Inclined Orbit Case  
    if ECC > ECC_min            %Elliptical Orbit
        AoP = acosd( dot(n_vec, ECC_vec)/(n*ECC) );
        if ECC_vec(3) < 0
            AoP = 360 - AoP;
        end
    else                        %Circular Orbit Case
        AoP = 0;
    end
else                            %Equatorial Orbit Case
    if ECC > ECC_min            %Elliptical Orbit 
        AoP = acosd(ECC_vec(1)/ECC);
        if ECC_vec(2) < 0
            AoP = 360 - AoP;
        end
    else                        %Circular Orbit 
        AoP = 0;
    end
end
%% Calculating the True Anomaly

if n ~= 0                       %Inclined Orbit Case  
    if ECC > ECC_min            %Elliptical Orbit
        TA = acosd(dot(ECC_vec,R)/(ECC*r));
        EA  = 2*atand(sqrt((1-ECC)/(1+ECC))*tand(TA/2));
        if dot(R,V) < 0
            TA = 360 - TA;
            EA = 360 - EA;
        end    
    else                        %Circular Orbit
        TA = acosd(dot(n_vec,R)/(n*r));
        if R(3) < 0
            TA = 360 - TA;
        end
        EA = TA;
    end
else                            %Equatorial Orbit Case
    if ECC > ECC_min            %Elliptical Orbit
        TA = acosd(dot(ECC_vec,R)/(ECC*r));
        EA  = 2*atand( sqrt((1-ECC)/(1+ECC))*tand(TA/2));
        if dot(R,V) < 0
            TA = 360 - TA;
            EA = 360 - EA;
        end
    else
        TA = acosd(R(1)/r)      %Circular Orbit
        if R(2) < 0
            TA = 360 - TA; 
        end
        EA = TA;
    end
end

%% Calculating other interesting data
rp  = (h^2/mu)*(1/(1+ECC));
ra  = (h^2/mu)*(1/(1-ECC));
T   = (2*pi*SMA^(3/2))/sqrt(mu);

%Assembling the Outputs
COE.SMA  = SMA;
COE.ECC  = ECC;
COE.INC  = INC;
COE.RAAN = RAAN;
COE.AoP  = AoP;
COE.TA   = TA;
Data.EA  = EA;
Data.rp  = rp;
Data.ra  = ra;
Data.T   = T;
end