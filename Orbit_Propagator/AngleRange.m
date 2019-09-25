%% AngleRange =============================================================
% Description: This function limits angles expressed in degrees to values
% between 0 and 360, to account for the periodicity of orbit calculations
% involving time. 

% Input:  angle_initial
% Output: angle_fixed 
%
% Created by: Cory Fraser - Feb.20, 2017
% Latest Edit: Cory Fraser - Feb.20, 2017
%==========================================================================
function angle_fixed = AngleRange(angle_initial)

if angle_initial >= 360
    angle_initial = angle_initial - fix(angle_initial/360)*360;
elseif angle_initial < 0
    angle_initial = angle_initial - (fix(angle_initial/360) - 1)*360;
end

angle_fixed = angle_initial;

end 
%==========================================================================
