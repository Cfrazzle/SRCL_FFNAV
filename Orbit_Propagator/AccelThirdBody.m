function r_ddot = AccelThirdBody(R_sat, R_mass, mu)
%% Third-body Acceleration ================================================
% Description: This function calculates the acceleration caused on a 
% spacecraft due a large point mass.
%
% Inputs:
%   R_sat  - Spacecraft position vector 
%   R_mass - Point mass position vector
%   mu     - Gravitational coefficient of point mass
%
% Outputs:
%   r_ddot - Acceleration of the spacecraft
%
%
% References:
%   Montenbruck, 'Satellite Orbits', pages 69-70, 2001
%
% Created by:  Cory Fraser - JUN 20, 2018
% Latest Edit: Cory Fraser - JUN 20, 2018
% Copyright(c) 2018 by Cory Fraser
% =========================================================================

% Relative position of satellite with respect to the mass
R_rel = R_sat - R_mass;

% Acceleration 
r_ddot = -mu*(R_rel/(norm(R_rel)^3)+R_mass/(norm(R_mass)^3));

end
