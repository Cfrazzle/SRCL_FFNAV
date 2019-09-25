function r_ddot = AccelSolRad(r_Sat,r_Earth,r_Moon,r_Sun,r_SunSSB, ...
        A_Solar, M_Sat, Cr, ShadowModel, P_Solar, AU, R_Earth);    
%% Solar Radiation Pressure ===============================================
% Description: This function calculates the worst-case acceleration caused 
% by solar radiation pressure, assumed acting on the full-frontal area of
% the spacecraft .
%
% Inputs:
%   r_sat       - Spacecraft position vector (ECI)
%   r_Earth     - Earth position vector (Barycentric)
%   r_Moon      - Moon position vector (ECI)
%   r_Sun       - Sun position vector (ECI)
%   r_SunSSB    - Sun position vector (Barycentric)
%   A_solar     - Cross-sectional area exposed to sunlight
%   M_sat       - Spacecraft mass
%   Cr          - Coefficient of radiation pressure
%   P_Solar     - Nominal Solar radiation pressure at 1 AU
%   AU          - Astronomical unit of length
%   ShadowModel - Specify type of shadow model to use

%        Moon wrt Earth          pccor      rpc
%        Earth wrt Sun           ccor       rc
%        Moon wrt Sun            pscor      rps   
%        Satellite wrt Earth     sbcor      rsb  
%        Satellite wrt Sun       bcor       rb 
%        Satellite wrt Moon      sbpcor     rsbp
%
% Outputs:
%   r_ddot - Acceleration of the spacecraft
%
% References:
%   Montenbruck, 'Satellite Orbits', pages 77-83, 2001
%
% Created by:  Cory Fraser - JUN 23, 2018
% Latest Edit: Cory Fraser - JUN 23, 2018
% Copyright(c) 2018 by Cory Fraser
% =========================================================================

%pccor = r_Moon;
%ccor = r_Earth-r_SunSSB;
%pscor = r_Moon-r_Sun;
%sbcor = r_Sat;
%bcor = r_Sat-r_Sun;
%sbpcor = r_Sat-r_Moon;

%Unit vector pointing from the Sun to the Satellite
e_SatSun = (r_Sat - r_Sun)/norm(r_Sat - r_Sun);

%Unit vector pointing from the Satellite to the Sun
e_SunSat = -e_SatSun;

%Evaluate shadow effects
if (ShadowModel)    % Geometric Shadow Model
    fprintf('This shadow effect is not ready \n')
    [nu,~] = Shadow(pccor,ccor,pscor,sbcor,bcor,sbpcor);
    return
    nu = 1;
else                %Cylindrical Shadow Model
    nu = AccelSolRadCyl(r_Sat,r_Sun, R_Earth);
    
    %Plot for visual check of Earth/Sun-vector intersection 
    %{
    coder.extrinsic('legend')
    coder.extrinsic('figure')

    if nu == 0
        figure(99)
        clf
        set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
        hold on
        %Plot Earth
        
            equat_rad = R_Earth;
            polar_rad = 6356.752e3;
            [xx yy zz]=ellipsoid (ccor(1),ccor(2),ccor(3),equat_rad, equat_rad, polar_rad);
            load('topo.mat','topo','topomap1');
            topo2 = [topo(:,181:360) topo(:,1:180)];
            pro.FaceColor= 'texturemap';
            pro.EdgeColor = 'none';
            pro.FaceLighting = 'gouraud';     
            pro.Cdata = topo2;
            pro.SpecularStrength = 0.5;      % change the strength of the reflected light
            earth = surface(xx,yy,zz,pro);
            colormap(topomap1) 
            xlabel(' X (km) ')
            ylabel(' Y (km) ')
            zlabel(' Z (km) ')
            zlimits = [min(topo(:)) max(topo(:))];
            demcmap(zlimits);
            view(45,15)
            light('Position',[0.75 0.75 1])     % add a light

            %Plot Spacecraft
            scatter3(bcor(1), bcor(2), bcor(3), 'filled', 'r')
            
            %Plot Sun Vector
            h1 = quiver3(bcor(1), bcor(2), bcor(3), 1000e3*e_SunSat(1), 1000e3*e_SunSat(2), 1000e3*e_SunSat(3));
            set(h1,'AutoScale','on', 'AutoScaleFactor', 10)
        legend('Earth', 'Spacecraft', 'Sun Vector')
    end
    %}
end
   
% Acceleration
r_ddot = -nu*P_Solar*Cr*(A_Solar/M_Sat)*AU*AU*(r_Sun-r_Sat)/norm(r_Sun-r_Sat)^3;

%Version 1 - Overestimates, doesn't account for fluctuation of P0
%r_ddot = nu*Cr*P_Solar*(A_Solar/M_Sat)*e_SatSun;  

%Equation from Montenbruck (assumes distance from Earth-to-Sun)
%r_ddot4 = -nu*Cr*P_Solar*(A_Solar/M_Sat)*AU*AU*r_Sun/norm(r_Sun)^3;


