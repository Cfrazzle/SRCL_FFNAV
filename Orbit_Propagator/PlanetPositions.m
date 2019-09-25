function [r_Mercury,r_Venus,r_Earth,r_Mars,r_Jupiter,r_Saturn,r_Uranus, ...
          r_Neptune,r_Pluto,r_Moon,r_Sun,r_SunSSB] = ...
          PlanetPositions(MJD, JD, Constants)
% Planetary Ephemeris Model ===============================================
% =========================================================================
% Description: This function determines the position of the Sun, the Moon, 
% and nine planets in the ECI frame, along with positions of the Sun in the
% barycentric frame. The position are determined using the epherides
% provided by NASA's Jet Propulsion Lab (JPL).
%
% Inputs:
%   MDJ_TDB      - Modified Julian Day of TDB 
%   Constants.   - Structure containing astrodynamic constants
%   PropOptions. - Structure containing orbit propagation options
%   t_sim        - Current simulation time
%
% Outputs:
%   r_Earth  - Position of the Earth [m]
%   r_SunSSB - Position of the Sun w.r.t the solar system barycenter [m]
%   r_body   - ECI position of a given body [m]   
%   r_Sun    - ECI position of the Sun referred to the International Celestial 
%           Reference Frame (ICRF))
%
% Reference: Montenbruck & Gill, Chapter 3, 2011
%
% Created by: Cory Fraser - JUN 30, 2018
% Last Edit : Cory Fraser - JUN 30, 2018
% Copyright(c) 2018 by Cory Fraser
% =========================================================================

%{
%Load planetary coefficients
PC = Constants.PC;

% Determine relevant planetary coefficients from Julian Date
idx = find(PC(:,1)<=JD & JD<=PC(:,2),1); 
PC_idx = PC(idx,:);
%}

%{
PC_idx = Constants.PC_idx;

%MJD start time of the interval
t1 = PC_idx(1)-2400000.5; 

%Difference in current and start times
delta_t = MJD - t1;

%% Calculate Earth Position

%Extract Coefficients for Earth
temp = (231:13:270);
Cx_Earth = PC_idx(temp(1):temp(2)-1);
Cy_Earth = PC_idx(temp(2):temp(3)-1);
Cz_Earth = PC_idx(temp(3):temp(4)-1);

temp = temp+39;
Cx = PC_idx(temp(1):temp(2)-1);
Cy = PC_idx(temp(2):temp(3)-1);
Cz = PC_idx(temp(3):temp(4)-1);
Cx_Earth = [Cx_Earth,Cx];
Cy_Earth = [Cy_Earth,Cy];
Cz_Earth = [Cz_Earth,Cz];    

if 0 <= delta_t && delta_t <= 16
    j = 0;
    MJD0 = t1;
elseif 16 < delta_t && delta_t <= 32
    j = 1;
    MJD0 = t1+16*j;
end

%Calculate position of the Earth using Chebyshev approximation [km]
r_Earth = Cheb3D(MJD, 13, MJD0, MJD0+16, Cx_Earth(13*j+1:13*j+13),...
                     Cy_Earth(13*j+1:13*j+13), Cz_Earth(13*j+1:13*j+13))';

% Convert position to [m]
r_Earth = r_Earth*1000;
%}

%% Assume Stationary Planets
r_Sun = [ 123994411196.868; 78220316016.5798; 33912542828.4664]; %Distance to the Sun [m]
r_SunSSB = [-86124067.1845607; -756419250.331927; -318619230.43305];
r_Moon = [-358572728.821806;  -32247804.1300803; 20051273.5230284]; %Distance to Moon [m]
r_Mercury = [97870600517.3498; 113049568850.835; 55226267603.4471];
r_Venus = [115845660227.47; 175889800129.773; 78368988662.5045];
r_Earth = [-124080535264.052; -78976735266.9117; -34231162058.8995];
r_Mars = [131723242958.197; 291049535045.856; 131320797354.225];
r_Jupiter = [-143743483334.507; 749045581183.697; 327966282389.645];
r_Saturn = [434488523553.922; 1301047301391.56; 525629876130.815];
r_Uranus = [2590834697311.75; -1458719277428.24;-674131034868.774];
r_Neptune = [2954417982591.55; -3136262953978.92; -1352249734730.92];
r_Pluto = [-965920513727.36; -4234252435980.76; -983381542696.103];
r_Moon = [-358572728.821806; -32247804.1300803; 20051273.5230284];


end

