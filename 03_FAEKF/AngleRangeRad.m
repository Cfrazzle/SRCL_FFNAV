%% AngleRangeRad ==========================================================
% Description: This function limits angles expressed in radians to values
% between 0 and 2pi, to account for the periodicity of orbit calculations
% involving time. 

% Input:  angle_initial
% Output: angle_fixed 
%
% Created by:  Cory Fraser - SEP 11, 2018
% Latest Edit: Cory Fraser - SEP 11, 2018
%==========================================================================
function angle_fixed = AngleRangeRad(angle_initial)

if angle_initial >= 2*pi
    angle_initial = angle_initial - fix(angle_initial/(2*pi))*(2*pi);
elseif angle_initial < 0
    angle_initial = angle_initial - (fix(angle_initial/(2*pi)) - 1)*(2*pi);
end

angle_fixed = angle_initial;

end 
%==========================================================================
