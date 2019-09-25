%% GPS Simulation - Klobuchar Ionospheric Error Model =====================
% =========================================================================
% Description: This 
%
% Inputs:
%   r_Rec - Receiver position vector in ECEF frame      [m]
%   r_GPS - GPS satellite position matrix in ECEF frame [m]
%
% Outputs:
%   R_iono - Range error due to ionospheric interference [m]
%   T_iono - Time error due to ionospheric interference  [s]
%
% Reference: Montenbruck & Gill, Page 212
%            Klobuchar, Page 329
%
% Created by: Cory Fraser - Jan 29, 2018
% Last Edit : Cory Fraser - Feb 01, 2018
%% ========================================================================
function [R_iono, t_iono] = Error_Iono_Klobuchar(r_Rec, r_GPS, alpha, beta, t_GPS, c)

format long g

% Convert receiver position into geodetic components
Rec_gLLA = ECEF_to_gLLA(r_Rec); 

% Convert latitude and longitude to semi-circles [rad]
lat = Rec_gLLA(1)/pi;
lon = Rec_gLLA(2)/pi;

% Compute azimuth and elevation [rad]
[A, E] = Calc_AzEl(r_Rec, r_GPS);

% ========================================================================
% Klobuchar Model
% ========================================================================

% Compute the Earth-centered angle [semi-circles]
psi = (0.0137/(E/pi+0.11) - 0.022)';
%psi = (445/(E*180/pi + 20) - 4)'/180; %Degrees

% Determine the IPP sub-ionospheric latitude [semi-circles]
phi_L = lat + psi.*cos(A);
    if any(phi_L > 0.416)
        phi_L(phi_L > 0.416) = 0.416;
    elseif any(phi_L < -0.416)
        phi_L(phi_L < -0.416) = -0.416;
    end

% Determine the sub-ionospheric longitude [semi-circles]
if abs(cos(phi_L*pi)) < 0.01
    fprintf('\n Lambda calculation almost singular \n')
end
lambda_L = lon + psi.*sin(A)./cos(phi_L*pi);
if isnan(lambda_L)
    fprintf('\n NAN Detected \n')
end

% Calculate the IPP geomagnetic latitude [semi-circles]
phi_M =  phi_L + 0.064*cos((lambda_L - 1.617)*pi);

%Find the local time of day [s]              
t = 43200*lambda_L + t_GPS;
%t = mod(t, 86400);
if any(t > 86400)
    t(t > 86400) = t(t > 86400) - 86400;
elseif any(t < 0)
    t(t < 0) = t(t < 0) + 86400;
end

% Compute the slant factor
F = 1.0 + 16.0*(0.53-E/pi).^3;

% Period of the vTEC cosine approximation [s]
PER = beta(1) + beta(2).*phi_M + beta(3).*phi_M.^2 + beta(4).*phi_M.^3;
if any(PER < 72000)
    PER(PER < 72000) = 72000;
end

% Phase of the vTEC cosine approximation [rad]
x = 2*pi*(t-50400)./PER;

% Amplitude of the vTEC cosine approximation [s]
AMP = alpha(1) + alpha(2).*phi_M + alpha(3).*phi_M.^2 + alpha(4).*phi_M.^3;
if any(AMP < 0)
    AMP(AMP < 0) = 0;
end

% Calculate the signal delay
t_iono = zeros(size(r_GPS,1),1);
phase_out = abs(x) > 1.57;
if any(phase_out)
    t_iono(phase_out) = F(phase_out)*(5*10^-9);
end
if any(~phase_out)
    t_iono(~phase_out) = F(~phase_out) ...
        .* (5*10^-9 + AMP(~phase_out) ...
        .* (1 - x(~phase_out).^2/2 + x(~phase_out).^4/24));
end


% =========================================================================
% Convert time delay to range error

R_iono = c*t_iono;
if any(R_iono < 0)
    R_iono = 0;
end
    

end
% =========================================================================