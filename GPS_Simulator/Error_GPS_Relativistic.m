%% GPS Simulation - GPS Relativistic Error Model ==========================
% =========================================================================
% Description: This function calculates the relativistic error caused in
% the onboard clock of a GPS satellite, 
%
% Inputs:
%   ECC - Eccentricity of the GPS satellite orbit   [-]
%   EA - Eccentric anomaly of the GPS satellite [deg]
%   SMA - Semi-major axis of the orbit [km] 
%
% Outputs:
%   R_rel - Range error due to relativity [m]
%   T_rel - Time error due to relativity  [s]
%
% Reference: Spilker, 2007 - Page 682-683
%            Hecimovic, 2013
%
% Created by: Cory Fraser - FEB 01, 2018
% Last Edit : Cory Fraser - JUN 20, 2018
%% ========================================================================
function [R_rel, t_rel] = Error_GPS_Relativistic(GPS_pos, GPS_vel, mu, c, varargin)

if length(varargin) ~= 2

    % Calculate the relativistic shift using position/velocity
    t_rel = (-2/c^2)*dot(GPS_pos, GPS_vel); 
    
else
    
    % Calculate the relativistic shift using orbital elements
    COE  = varargin{1};
    Data = varargin{2};
    
    ECC = COE.ECC; 
    SMA = COE.SMA;
    EA  = Data.EA;
    
    F = -4.442807633393060e-10;     % F = -2*sqrt(mu)/c^2; [s/m^1/2]
    t_rel = F*ECC*sqrt(SMA)*sind(EA);
    
end

%Calculate the corresponding range shift
R_rel = c*t_rel;

end