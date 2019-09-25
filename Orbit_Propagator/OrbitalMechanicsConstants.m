%% Orbital Mechanics Constants ============================================
% Description: This script contains constants values relevant to orbital
% mechanics, spacecraft navigation, and orbit propagation.
%
% Outputs:
%   Constants. - A structure containing all the constants
%
% Other Functions Called:
%   FFNAV_FAEKF_Params.m - Initialization of FAEKF parameters 
%
% References:
%     
% 
% Created by:  Cory Fraser - JUN 20, 2018
% Latest Edit: Cory Fraser - JUN 30, 2018
% Copyright(c) 2018 by Cory Fraser
% =========================================================================

%% Physical Constants
ConstantsModel = 'WGS84 Model';

%Gravitational Parameters - mu, GM, gm
Constants.mu_Earth   = 3.98600435436e14;        % [m^3/s^2]; WGS84
Constants.mu_Sun     = 132712440041.939377e9; % [m^3/s^2]; DE436
Constants.mu_Moon    = 4902.800117e9;         % [m^3/s^2]; DE436
Constants.mu_Mercury = 22031.780000e9;        % [m^3/s^2]; DE436
Constants.mu_Venus   = 324858.592000e9;       % [m^3/s^2]; DE436
Constants.mu_Mars    = 42828.375214e9; 	  	 % [m^3/s^2]; DE436
Constants.mu_Jupiter = 126712764.133446e9;    % [m^3/s^2]; DE436
Constants.mu_Saturn  = 37940585.200000e9;     % [m^3/s^2]; DE436
Constants.mu_Uranus  = 5794556.465752e9;      % [m^3/s^2]; DE436
Constants.mu_Neptune = 6836527.100580e9;      % [m^3/s^2]; DE436
Constants.mu_Pluto   = 975.501176e9;   		 % [m^3/s^2]; DE436

% Earth Parameters
Constants.J2        = 1082.64e-6;        %Earth's J2 Coefficient
Constants.R_Earth   = 6378.137e3;        %Equatorial Radius - WGS84 [m]
Constants.R_Earth_p = 6356.752e3;        %Polar Radius - WGS84 [m]
Constants.w_Earth   = 7.2921151467e-5;   %Mean angular velocity - WGS84 [rad/sec]
Constants.f_Earth   = 1/298.257223563;   %Flattening; WGS-84

% Other Data
Constants.R_Sun     = 696000e3;            % Sun's radius [m]; DE436
Constants.R_Moon    = 1738e3;              % Moon's radius [m]; DE436
Constants.c_light   = 299792457.999999984; % Speed of light  [m/s]; DE436
Constants.AU        = 149597870699.999988; % Astronomical unit [m]; DE436

Constants.MJD_J2000 = 51544.5;             % Modified Julian Date of J2000
Constants.T_B1950   = -0.500002108;        % Epoch B1950

% Solar radiation pressure at 1 AU
Constants.P_Solar = 1367/Constants.c_light; % [N/m^2] (1367 W/m^2); IERS

%Conversion Factors
deg2rad          = pi/180;               % Degrees to radians
rad2deg          = 180/pi;               % Radians to degrees 

%Planetary Ephemeride Coefficients for DE436 
%   load DE436Coeff.mat
%   PC = DE436Coeff;


%% Propagator Time - Universal Coordinated (UTC)
Day     = 29;   %Range: 1 - 31
Month   = 11;   %Range: 1-12
Year    = 2018; %Range 1901-2099

Hour    = 0;    %Range: 0-23
Minute  = 0;    %Range: 0-60
Second  = 0;    %Range: 0-60

Constants.JD = JulianDate(Day,Month,Year,Hour,Minute,Second); %From Example 12.10
Constants.Rvec_Sun = PositionSun(Constants.JD);
Constants.Rvec_Moon = PositionMoon(Constants.JD);
%% Output =================================================================   
filename = 'AstroConstants.mat';
pathname = fileparts('C:\Users\Cory\Desktop\FFNAV Data\');
file_out = fullfile(pathname, filename);
save(file_out, 'Constants', 'ConstantsModel', 'deg2rad','rad2deg')...
    %,'DE436Coeff', 'PC');
% =========================================================================
