function [r_ddot] = Accelerations(R_Sat, V_Sat, Constants, PropOptions, t_sim)
%% Orbit Propagation Force Model ==========================================
% =========================================================================
% Description: This function calculates the accelerations acting on a 
% spacecraft due to various perturbing forces.
%
% Inputs:
%   R_Sat        - Position vector of the spacecraft in ECI frame [m]
%   V_Sat        - Velocity vector of the spacecraft in ECI frame [m/s]
%   Constants.   - Structure containing astrodynamic constants
%   PropOptions. - Structure containing orbit propagation options
%   t_sim        - Current simulation time
%
% Outputs:
%   r_ddot - Acceleration of the spacecraft [m/s^2]
%
% Other Function Calls:
%   AccelThirdBody             - Calculates acceleration due to point mass
%   AccelRelativity            - Calculates acceleration due to relativity
%   AccelSolRad                - Calculates acceleration due to radiation
%   TBD - AccelDrag            - Calculates acceleration due to drag
%   TBD - AccelHarmonic_EEarth - Calculates acceleration w/ Elastic Earth
%   TBD - AccelHarmonic_AEarth - Calculates acceleration w/ Anelastic Earth
%
% Reference: Montenbruck & Gill, Chapter 3, 2011
%  
% Created by: Cory Fraser - JUN 20, 2018
% Last Edit : Cory Fraser - JUL 27, 2018
% Copyright(c) 2018 by Cory Fraser
% =========================================================================
mu_Earth = Constants.mu_Earth;
R_Earth  = Constants.R_Earth;

%% Two-Body Acceleration w/ J2 Option

r_ddot = -(mu_Earth/(norm(R_Sat))^3)*R_Sat;

if (PropOptions.J2)
    J2 = Constants.J2;
    delta_AccelJ2 = (mu_Earth/(norm(R_Sat))^3) * ...
    [ R_Sat(1)*(3/2)*J2*(R_Earth/norm(R_Sat))^2*(5*R_Sat(3)^2/norm(R_Sat)^2 - 1)
      R_Sat(2)*(3/2)*J2*(R_Earth/norm(R_Sat))^2*(5*R_Sat(3)^2/norm(R_Sat)^2 - 1)
      R_Sat(3)*(3/2)*J2*(R_Earth/norm(R_Sat))^2*(5*R_Sat(3)^2/norm(R_Sat)^2 - 3)];
  
  %Alternate fector form from DeRuiter, Page 164
  %delta_AccelJ2 = (3*mu_Earth*J2*R_Earth^2)/(2*norm(R_Sat)^5) * ...
  %  (( 5/norm(R_Sat)^2 * dot(R_Sat, [0; 0; 1])^2 -1 ) *R_Sat - ...
  %  2*dot(R_Sat, [0; 0; 1])*[0; 0; 1]); 
    
  r_ddot = r_ddot + delta_AccelJ2;
end

%% Harmonic Gravity Field Accelerations
%{

global const AuxParam eopdata

MJD_UTC = AuxParam.Mjd_UTC+t/86400;
[x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC] = IERS(eopdata,MJD_UTC,'l');
[UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC] = timediff(UT1_UTC,TAI_UTC);

JD = MJD_UTC+2400000.5;
[year, month, day, hour, minute, sec] = invjday(JD);
[DJMJD0, DATE] = iauCal2jd(year, month, day);
TIME = (60*(60*hour + minute) + sec)/86400;
UTC = DATE + TIME;
TT = UTC + TT_UTC/86400;
TUT = TIME + UT1_UTC/86400;
UT1 = DATE + TUT;

% Polar motion matrix (TIRS->ITRS, IERS 2003)
Pi = iauPom00(x_pole, y_pole, iauSp00(DJMJD0, TT));
% Form bias-precession-nutation matrix
NPB = iauPnm06a(DJMJD0, TT);
% Form Earth rotation matrix
gast = iauGst06(DJMJD0, UT1, DJMJD0, TT, NPB);
Theta  = iauRz(gast, eye(3));
% ICRS to ITRS transformation
E = Pi*Theta*NPB;

% Difference between ephemeris time and universal time
% JD = MJD_UTC+2400000.5;
% [year, month, day, hour, minute, sec] = invjday(JD);
% days = finddays(year, month, day, hour, minute, sec);
% ET_UT = ETminUT(year+days/365.25);
% MJD_ET = MJD_UTC+ET_UT/86400;
% [r_Mercury,r_Venus,r_Earth,r_Mars,r_Jupiter,r_Saturn,r_Uranus, ...
%  r_Neptune,r_Pluto,r_Moon,r_Sun,r_SunSSB] = JPL_Eph_DE436(MJD_ET);

MJD_TDB = Mjday_TDB(TT);
[r_Mercury,r_Venus,r_Earth,r_Mars,r_Jupiter,r_Saturn,r_Uranus, ...
 r_Neptune,r_Pluto,r_Moon,r_Sun,r_SunSSB] = JPL_Eph_DE436(MJD_TDB);

if (PropOptions.GravEl_Earth) %Elastic Earth Model
    r_ddot = AccelHarmonic_ElasticEarth(MJD_UTC,r_Sun,r_Moon,R_Sat,E,UT1_UTC,TT_UTC,x_pole,y_pole);
else                    %Anelastic Earth Model
    r_ddot = AccelHarmonic_AnelasticEarth(MJD_UTC,r_Sun,r_Moon,R_Sat,E,UT1_UTC,TT_UTC,x_pole,y_pole);
end
%}
% =========================================================================
%% Luni-Solar Third-Body Perturbations

% Evaluate position of the planets at current time step
%JD = PropOptions.JD_UTC_start+t_sim/86400;
%MJD = PropOptions.MJD_UTC_start+t_sim/86400;

%[r_Mercury,r_Venus,r_Earth,r_Mars,r_Jupiter,r_Saturn,r_Uranus, ...
% r_Neptune,r_Pluto,r_Moon,r_Sun,r_SunSSB] = PlanetPositions(MJD, JD, Constants);

%Assume fixed Planets, positions evaluated at simulation start time
%{
r_Mercury   = PropOptions.r_Mercury;
r_Venus     = PropOptions.r_Venus;
r_Earth     = PropOptions.r_Earth;
r_Mars      = PropOptions.r_Mars;
r_Jupiter   = PropOptions.r_Jupiter;
r_Saturn    = PropOptions.r_Saturn;
r_Uranus    = PropOptions.r_Uranus;
r_Neptune   = PropOptions.r_Neptune;
r_Pluto     = PropOptions.r_Pluto;
r_Moon      = PropOptions.r_Moon;
r_Sun       = PropOptions.r_Sun;
r_SunSSB    = PropOptions.r_SunSSB;
%}

%r_Sun = [ 123994411196.868; 78220316016.5798; 33912542828.4664]; %Distance to the Sun [m]
r_Sun = Constants.Rvec_Sun;
%r_Moon = [-358572728.821806; -32247804.1300803; 20051273.5230284];
r_Moon = Constants.Rvec_Moon;

r_SunSSB = [-86124067.1845607; -756419250.331927; -318619230.43305];
r_Mercury = [97870600517.3498; 113049568850.835; 55226267603.4471];
r_Venus = [115845660227.47; 175889800129.773; 78368988662.5045];
r_Earth = [-124080535264.052; -78976735266.9117; -34231162058.8995];
r_Mars = [131723242958.197; 291049535045.856; 131320797354.225];
r_Jupiter = [-143743483334.507; 749045581183.697; 327966282389.645];
r_Saturn = [434488523553.922; 1301047301391.56; 525629876130.815];
r_Uranus = [2590834697311.75; -1458719277428.24;-674131034868.774];
r_Neptune = [2954417982591.55; -3136262953978.92; -1352249734730.92];
r_Pluto = [-965920513727.36; -4234252435980.76; -983381542696.103];

if (PropOptions.Sun)
    delta_ThirdBodySun = AccelThirdBody(R_Sat,r_Sun,Constants.mu_Sun);
    r_ddot = r_ddot + delta_ThirdBodySun;
end

if (PropOptions.Moon)
    delta_ThirdBodyMoon = AccelThirdBody(R_Sat,r_Moon,Constants.mu_Moon);
    r_ddot = r_ddot + delta_ThirdBodyMoon;
end              
% =========================================================================
%% Planetary Third-Body Perturbations

if (PropOptions.Planets)
    delta_ThirdBodyMercury = AccelThirdBody(R_Sat,r_Mercury,Constants.mu_Mercury);
    delta_ThirdBodyVenus = AccelThirdBody(R_Sat,r_Venus,Constants.mu_Venus);
    delta_ThirdBodyMars =  AccelThirdBody(R_Sat,r_Mars,Constants.mu_Mars);
    delta_ThirdBodyJupiter = AccelThirdBody(R_Sat,r_Jupiter,Constants.mu_Jupiter);
    delta_ThirdBodySaturn = AccelThirdBody(R_Sat,r_Saturn,Constants.mu_Saturn);
    delta_ThirdBodyUranus = AccelThirdBody(R_Sat,r_Uranus,Constants.mu_Uranus);    
    delta_ThirdBodyNeptune =  AccelThirdBody(R_Sat,r_Neptune,Constants.mu_Neptune);
    delta_ThirdBodyPluto = AccelThirdBody(R_Sat,r_Pluto,Constants.mu_Pluto);
    
    r_ddot  =   r_ddot + delta_ThirdBodyMercury + delta_ThirdBodyVenus + ...
                delta_ThirdBodyMars + delta_ThirdBodyJupiter + ...
                delta_ThirdBodySaturn + delta_ThirdBodyUranus + ...
                delta_ThirdBodyNeptune + delta_ThirdBodyPluto;
end

% =========================================================================
%% Atmospheric drag
%
if (PropOptions.Drag)
    
    % Atmospheric density
    rho_atm = Density_HarrisPriester(R_Sat, r_Sun, Constants); %[kg/m^3]
        
    %Earth's Angular Velocity
    Omega = Constants.w_Earth;
    % Omega = 7292115.8553e-11+4.3e-15*( (MJD_UTC-Constant.MJD_J2000)/36525 ); % [rad/s]
    %Omega = Constant.w_Earth-0.843994809*1e-9*LOD; % IERS [rad/s]
    
    %Bias-Precession-Nutation Transformation
    NPB = eye(3); %Assuming small time frame
    
    delta_Drag = AccelDrag(rho_atm,R_Sat,V_Sat,NPB, ...
        PropOptions.A_Drag,PropOptions.M_Sat,PropOptions.Cd,Omega);
    r_ddot = r_ddot + delta_Drag;
end

% =========================================================================
%% Solar Radiation Pressure

if (PropOptions.SolarRad)
    delta_SolarRad = AccelSolRad(R_Sat,r_Earth,r_Moon,r_Sun,r_SunSSB, ...
        PropOptions.A_Solar,PropOptions.M_Sat,PropOptions.Cr,...
        PropOptions.ShadowModel, Constants.P_Solar, ...
        Constants.AU, Constants.R_Earth);
    
    r_ddot = r_ddot + delta_SolarRad;
end

% =========================================================================
%% Relativistic Effects

if (PropOptions.Relativity) 
    mu_Earth = Constants.mu_Earth;
    c_light  = Constants.c_light;
    
    delta_Relativity = AccelRelativity(R_Sat,V_Sat, mu_Earth, c_light); 
    r_ddot = r_ddot + delta_Relativity;
end

end
