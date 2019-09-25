%% GPS Simulation - Klobuchar Ionospheric Error Model =====================
% =========================================================================
% Description: This function calculates the clock time for specific GPS
% satellites given the GPS system time and the clock correction
% coefficients.
%
% Inputs:
%   t_ref   - Satellite clock reference epoch [s]
%   t_GPS   - GPS system time at transmission [s]
%   alpha   - Clock correction coefficients 
%              1) Clock bias        [s]
%              2) Frequency bias    [s/s]
%              3) Frequency drift   [s/s^2]
%
% Outputs:
%   R_ClockOffset - Range error due to ionospheric interference [m]
%   T_ClockOffset - Time error due to ionospheric interference  [s]
%
% Reference: Montenbruck & Gill, Page 217-218
%            Vigneron, Page 43
%            GPS IS-200H, Page 101,per 20.3.3.3.3.1    
%
% Created by: Cory Fraser - JAN 29, 2018
% Last Edit : Cory Fraser - JUN 20, 2018
%% ========================================================================
function [R_ClockOffset, t_ClockOffset] = Error_GPS_Clock_Offset(af, t_GPS, t_ref,c)

%Coefficients present for all measured GPS satellites
if size(af,1) ~=size(t_GPS,2)
    error('Not enough correction coefficients for GPS satellites')
    return
end
%Difference between GPS time and reference, reset every two hours
delta_t = t_GPS-t_ref;
delta_t = mod(delta_t, 7200);

% Calculating the GPS satellite time offset from GPS system time
%t_Sat = t_GPS(:) + alpha(:,1) + alpha(:,2).*(delta_t(:)) + alpha(:,3).*(delta_t(:)).^2;
%t_ClockOffset = t_Sat-t_GPS';

% Calculating clock offset directly
t_ClockOffset = af(:,1) + af(:,2).*(delta_t(:)) + af(:,3).*(delta_t(:)).^2;

% Caclulating the range error due to the clock offset
R_ClockOffset = c*t_ClockOffset;
end