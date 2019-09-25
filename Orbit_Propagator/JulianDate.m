function [JD, T0] = JulianDate(day, month, year, hour, minute, second)
%% Julian Date ============================================================
% Description: This function is used to determine the Julian Day JD based
% on the user-provided date and time. 
%
% Inputs:
%   day - range: 1 - 31
%   month - range: 1 - 12
%   year - range: 1901 - 2099
%   hour - range: 0 - 23 (Universal Time)
%   minute - rage: 0 - 60
%   second - range: 0 - 60
%
% Outputs:
%   JD - Julian day number at specified UT
%   T0 - Number of Julian centuries between date and J2000
%
% Created by:  Cory Fraser - FEB 20, 2017
% Latest Edit: Cory Fraser - JUN 30, 2018
%==========================================================================
if (nargin < 4)
    hour    = 0;
    minute  = 0;
    second  = 0;
end 

UT = hour + minute/60 +second/(60*60);

J0 = 367*year-fix(7*(year+fix((month+9)/12))/4)+fix(275*month/9)+day+1721013.5;

JD = J0 + UT/24;

T0 = (JD - 2451545.0)/36525;
%{
fprintf('Julian Date Conversion \n')
fprintf(' Input data:\n');
fprintf(' Year   = %i \n', year)
fprintf(' Month  = %i \n', month)
fprintf(' Day    = %i \n', day)
fprintf(' Hour   = %i', hour)
fprintf(' Minute = %i', minute)
fprintf(' Second = %i\n', second)
fprintf(' Julian day (midnight) = %11.3f \n', J0);
fprintf(' Julian day number     = %11.3f \n', JD);
fprintf(' Julian centuries      = %11.3f \n', T0);
fprintf('============================================================================== \n')
%}
end
