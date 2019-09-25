function r_ddot = AccelRelativity(R_Sat, V_Sat, mu_Earth, c_light)
%% OrbitPropagator ========================================================
% Description: This function calculates the acceleration caused on a 
% spacecraft due to relativity.
%
% Inputs:
%   R_Sat - Spacecraft position vector 
%   V_Sat - Spacecraft Velocity vector 
%   mu_Earth - Gravitational coefficient of Earth
%
% Outputs:
%   r_ddot - Acceleration of the spacecraft
%
% References:
%   Montenbruck, 'Satellite Orbits', pages 110-111, 2001
%   IERS Technical Note 03: 15. General Relativistic Dynamical Model
%
% Created by:  Cory Fraser - JUN 20, 2018
% Latest Edit: Cory Fraser - JUN 20, 2018
% Copyright(c) 2018 by Cory Fraser
% =========================================================================

r_Sat = norm(R_Sat);
v_Sat = norm(V_Sat);

e_r = R_Sat/r_Sat;
e_v = V_Sat/v_Sat;
    
r_ddot = -(mu_Earth/(r_Sat^2))*...
        (( 4*mu_Earth/(c_light^2*r_Sat)- v_Sat^2/(c_light^2))*e_r...
         + 4*(v_Sat/c_light)^2*dot(e_r,e_v)*e_v);
end
