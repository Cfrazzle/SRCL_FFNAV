%% FFNAV GPS Parameters ===================================================
% ========================================================================
% Description: This script constructs the initial conditions for the 
% GPS constellation
%
% Inputs:
%   
% Created by: Cory Fraser - JAN 22, 2018
% Last Edit : Cory Fraser - JUN 20, 2018
%% ========================================================================
clear
close

%General system constants
mu              = 3.986004418e14;       %Earth's Gravitational Constant - WGS84 [m^3/s^2]
R_Earth         = 6378.137e3;           %Equatorial Radius of the Earth [m]
R_Earth_p       = 6356.752e3;           %Polar Radius of the Earth [m]
w_Earth         = 7.2921151467e-5;      %Mean angular velocity of the Earth - WGS84 [rad/sec]

%General GPS Configuration
ECC_GPS         = 0;                    %Eccentricity of GPS orbit [m]
SMA_GPS         = 26559.710e3;          %Semi-major axis of GPS orbit [m]
INC_GPS         = 55*pi/180;            %Inclination of GPS orbit [rad]
n_GPS           = sqrt(mu/SMA_GPS^3);   %Mean orbital motion [rad/s]
T_GPS           = 2*pi/n_GPS;           %Orbital Period [s]
c               = 2.99792458e8;         %Speed of light in Vacuum [m/s]
deg2rad         = pi/180;               %Coversion for degrees to radians
rad2deg         = 180/pi;               % Coversion for radians to degrees 

%% GPS Parameters
% Data extracted from GPS ICD-200, Epoch of 00:00:00 on Jul 1, 1993
%                M/theta          SMA       ecc       inc           RAAN           AoP 
GPS_COE = [  268.126*deg2rad    SMA_GPS   ECC_GPS   INC_GPS   357.734*deg2rad   0*deg2rad
             161.786*deg2rad    SMA_GPS   ECC_GPS   INC_GPS   357.734*deg2rad   0*deg2rad
              11.676*deg2rad    SMA_GPS   ECC_GPS   INC_GPS   357.734*deg2rad   0*deg2rad
              41.806*deg2rad    SMA_GPS   ECC_GPS   INC_GPS   357.734*deg2rad   0*deg2rad
                %
              94.916*deg2rad    SMA_GPS   ECC_GPS   INC_GPS    57.734*deg2rad   0*deg2rad
              66.356*deg2rad    SMA_GPS   ECC_GPS   INC_GPS    57.734*deg2rad   0*deg2rad
             173.336*deg2rad    SMA_GPS   ECC_GPS   INC_GPS    57.734*deg2rad   0*deg2rad
             309.976*deg2rad    SMA_GPS   ECC_GPS   INC_GPS    57.734*deg2rad   0*deg2rad
             204.376*deg2rad    SMA_GPS   ECC_GPS   INC_GPS    57.734*deg2rad   0*deg2rad
                %
             111.876*deg2rad    SMA_GPS   ECC_GPS   INC_GPS   117.734*deg2rad  0*deg2rad
              11.796*deg2rad    SMA_GPS   ECC_GPS   INC_GPS   117.734*deg2rad  0*deg2rad
             339.666*deg2rad    SMA_GPS   ECC_GPS   INC_GPS   117.734*deg2rad  0*deg2rad
             241.556*deg2rad    SMA_GPS   ECC_GPS   INC_GPS   117.734*deg2rad  0*deg2rad
                %
             135.226*deg2rad    SMA_GPS   ECC_GPS   INC_GPS   177.734*deg2rad  0*deg2rad
             282.676*deg2rad    SMA_GPS   ECC_GPS   INC_GPS   177.734*deg2rad  0*deg2rad
             257.976*deg2rad    SMA_GPS   ECC_GPS   INC_GPS   177.734*deg2rad  0*deg2rad
              35.156*deg2rad    SMA_GPS   ECC_GPS   INC_GPS   177.734*deg2rad  0*deg2rad
             167.356*deg2rad    SMA_GPS   ECC_GPS   INC_GPS   177.734*deg2rad  0*deg2rad
                %
             197.046*deg2rad    SMA_GPS   ECC_GPS   INC_GPS   237.734*deg2rad  0*deg2rad
             302.596*deg2rad    SMA_GPS   ECC_GPS   INC_GPS   237.734*deg2rad  0*deg2rad
              66.066*deg2rad    SMA_GPS   ECC_GPS   INC_GPS   237.734*deg2rad  0*deg2rad
             333.686*deg2rad    SMA_GPS   ECC_GPS   INC_GPS   237.734*deg2rad  0*deg2rad
                %
             238.886*deg2rad    SMA_GPS   ECC_GPS   INC_GPS   297.734*deg2rad  0*deg2rad
               0.456*deg2rad    SMA_GPS   ECC_GPS   INC_GPS   297.734*deg2rad  0*deg2rad
             334.016*deg2rad    SMA_GPS   ECC_GPS   INC_GPS   297.734*deg2rad  0*deg2rad
             105.206*deg2rad    SMA_GPS   ECC_GPS   INC_GPS   297.734*deg2rad  0*deg2rad
             135.346*deg2rad    SMA_GPS   ECC_GPS   INC_GPS   297.734*deg2rad  0*deg2rad ];

%Assembling initial vector to be passed to the EKF simulation (1x118) 
filename        = 'FFNAV_GPS_Params.mat';
save(filename)
fprintf('\n GPS file %s was created. \n', filename)

load('C:\Users\Cory\Desktop\FFNAV Data\FFNAV_PRISMA_Mod_Params');


% GPS Constellation Simulation Setup
time_step   = 1;        % Time step for the simulation
time_end    = T_GPS*1;  % Duration of the simulation (~One Complete GPS Orbit)

% Open the Simulink Model, set simulation parameters
load_system('GPS_Constellation')
if strcmp(get_param('GPS_Constellation', 'shown'), 'off')
    open_system('GPS_Constellation');
end
set_param('GPS_Constellation', 'Solver', 'ode4')
set_param('GPS_Constellation', 'SolverType', 'Fixed-step')
set_param('GPS_Constellation', 'FixedStep', 'time_step')
set_param('GPS_Constellation', 'RelTol','1e-15')
set_param('GPS_Constellation', 'AbsTol','1e-15')
set_param('GPS_Constellation', 'StartTime', '0')
set_param('GPS_Constellation', 'StopTime', 'time_end')

% Run the EKF Simulation
fprintf('\nStarting the GPS simulation...\n')
tic;
sim('GPS_Constellation');
toc;
fprintf('GPS simulation complete.\n')
fprintf('--------------------------------------------------------------\n')


%% Plots for EKF ==========================================================

set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', 20)
set(0,'DefaultTextFontName', 'Times New Roman')
set(0,'DefaultTextFontSize', 24)
set(0,'defaultLineLineWidth', 2)

time = tout/3600;   %Simulation time from seconds to hours, for plots

%Figure: 3D Position Plot for GPS
fig_id = 1;
figure(fig_id); fig_id = fig_id + 1;    
set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
grid on
hold on
plot3(GPS_A1(:,1),GPS_A1(:,2),GPS_A1(:,3));
plot3(r_target_ECI(:,1),r_target_ECI(:,2),r_target_ECI(:,3), '-.r');
view(150,10)
axis equal
lim_ECI = 1.25*SMA_GPS;
title('GPS Constellation Orbits - ECI Frame')
set(gca, 'XLim', 1.1*[-lim_ECI lim_ECI], 'YLim', ...
    1.1*[-lim_ECI lim_ECI], 'ZLim', ...
    1.1*[-lim_ECI lim_ECI]);
legend('GPS A1', 'Target') 

%Other Planes
%{
plot3(GPS_B1F(:,1),GPS_B1F(:,2),GPS_B1F(:,3));
plot3(GPS_C1(:,1),GPS_C1(:,2),GPS_C1(:,3));
plot3(GPS_D1(:,1),GPS_D1(:,2),GPS_D1(:,3));
plot3(GPS_E1(:,1),GPS_E1(:,2),GPS_E1(:,3));
plot3(GPS_F1(:,1),GPS_F1(:,2),GPS_F1(:,3));
legend('A', 'B', 'C', 'D', 'E', 'F') 
%}

%Plot Earth
equat_rad = R_Earth;
polar_rad = R_Earth_p;
[xx yy zz]=ellipsoid (0,0,0,equat_rad, equat_rad, polar_rad);
load('topo.mat','topo','topomap1');
topo2 = [topo(:,181:360) topo(:,1:180)];
pro.FaceColor= 'texture';
pro.EdgeColor = 'none';
pro.FaceLighting = 'phong';     
pro.Cdata = topo2;
earth= surface(xx,yy,zz,pro,'HandleVisibility','off');
colormap(topomap1) 
xlabel(' X (km) ')
ylabel(' Y (km) ')
zlabel(' Z (km) ')

line([0 lim_ECI],[0 0],[0 0],'Color', 'black', 'Marker','.','LineStyle','-','LineWidth', 2.5, 'HandleVisibility','off')
line([0 0],[0 lim_ECI],[0 0],'Color', 'black', 'Marker','.','LineStyle','-','LineWidth', 2.5, 'HandleVisibility','off')
line([0 0],[0 0],[0 lim_ECI],'Color', 'black', 'Marker','.','LineStyle','-','LineWidth', 2.5, 'HandleVisibility','off')
set(gca,'xticklabel',num2str(get(gca,'xtick')'/1000))
set(gca,'yticklabel',num2str(get(gca,'ytick')'/1000))
set(gca,'zticklabel',num2str(get(gca,'ztick')'/1000))   

%All Other Satellites
%{
plot3(GPS_A2(:,1),GPS_A2(:,2),GPS_A2(:,3));
plot3(GPS_A3(:,1),GPS_A3(:,2),GPS_A3(:,3));
plot3(GPS_A4(:,1),GPS_A4(:,2),GPS_A4(:,3));

plot3(GPS_B1F(:,1),GPS_B1F(:,2),GPS_B1F(:,3));
plot3(GPS_B1A(:,1),GPS_B1A(:,2),GPS_B1A(:,3));
plot3(GPS_B2(:,1),GPS_B2(:,2),GPS_B2(:,3));
plot3(GPS_B3(:,1),GPS_B3(:,2),GPS_B3(:,3));
plot3(GPS_B4(:,1),GPS_B4(:,2),GPS_B4(:,3));

plot3(GPS_C1(:,1),GPS_C1(:,2),GPS_C1(:,3));
plot3(GPS_C2(:,1),GPS_C2(:,2),GPS_C2(:,3));
plot3(GPS_C3(:,1),GPS_C3(:,2),GPS_C3(:,3));
plot3(GPS_C4(:,1),GPS_C4(:,2),GPS_C4(:,3));

plot3(GPS_D1(:,1),GPS_D1(:,2),GPS_D1(:,3));
plot3(GPS_D2F(:,1),GPS_D2F(:,2),GPS_D2F(:,3));
plot3(GPS_D2A(:,1),GPS_D2A(:,2),GPS_D2A(:,3));
plot3(GPS_D3(:,1),GPS_D3(:,2),GPS_D3(:,3));
plot3(GPS_D4(:,1),GPS_D4(:,2),GPS_D4(:,3));

plot3(GPS_E1(:,1),GPS_E1(:,2),GPS_E1(:,3));
plot3(GPS_E2(:,1),GPS_E2(:,2),GPS_E2(:,3));
plot3(GPS_E3(:,1),GPS_E3(:,2),GPS_E3(:,3));
plot3(GPS_E4(:,1),GPS_E4(:,2),GPS_E4(:,3));

plot3(GPS_F1(:,1),GPS_F1(:,2),GPS_F1(:,3));
plot3(GPS_F2F(:,1),GPS_F2F(:,2),GPS_F2F(:,3));
plot3(GPS_F2A(:,1),GPS_F2A(:,2),GPS_F2A(:,3));
plot3(GPS_F3(:,1),GPS_F3(:,2),GPS_F3(:,3));
plot3(GPS_F4(:,1),GPS_F4(:,2),GPS_F4(:,3));
%}

%% =============================================================
%Plot Pseudorange Errors
figure(fig_id); fig_id = fig_id + 1;
    
    % Figure: Ionospheric Errors
        subplot(3,1,1)
        hold on
        plot(time, R_iono)   
        ylabel('Ionospheric Errors (m)')
        set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
        grid on
        xlim([0 inf])
     
    % Figure: Relativity Errors
        subplot(3,1,2)
        hold on
        plot(time, R_rel) 
        ylabel('Relativity Errors (m)')
        set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
        grid on
        xlim([0 inf])
        
    % Figure: Clock Offset Errors
        subplot(3,1,3)
        hold on
        plot(time, R_ClockOff) 
        ylabel('Clock Offset Errors (m)')
        xlabel('Time (hr)')
        set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
        grid on
        xlim([0 inf])
        
     %Total Errors   
        figure(fig_id); fig_id = fig_id + 1;
        hold on
        plot(time, R_total)   
        title('Total Pseudorange Error')
        xlabel('Time (hrs)')
        ylabel('Range Error (m)')
        set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
        grid on
        xlim([0 inf])
        
%% =============================================================
% Orbit Animations
%{
figure(30)
title('GPS Positions', ... 
    'Fontweight', 'bold', 'FontSize', 12)
axis equal
curve1 = animatedline('color','r', 'LineWidth', 0.5);
curve2 = animatedline('color','b', 'LineWidth', 0.5);
curve3 = animatedline('color','g', 'LineWidth', 0.5);
curve4 = animatedline('color','k', 'LineWidth', 0.5);
curve5 = animatedline('color','r', 'LineWidth', 0.5);
curve6 = animatedline('color','b', 'LineWidth', 0.5);
curve7 = animatedline('color','g', 'LineWidth', 0.5);
curve8 = animatedline('color','k', 'LineWidth', 0.5);
curve9 = animatedline('color','c', 'LineWidth', 0.5);

set(gca, 'XLim', [-30000 30000], 'YLim', ...
    [-30000 30000], 'ZLim', [-30000 30000]); 
xlabel(' x (m) ')
ylabel(' y (m) ')
zlabel(' z (m) ')
view(-40, 30);
set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
grid on
hold on;
GPS_A1path  = scatter3(GPS_A1(1,1), GPS_A1(1,2), GPS_A1(1,3), 'MarkerEdgeColor', 'r');
GPS_A2path  = scatter3(GPS_A2(1,1), GPS_A2(1,2), GPS_A2(1,3), 'MarkerEdgeColor', 'b');
GPS_A3path  = scatter3(GPS_A3(1,1), GPS_A3(1,2), GPS_A3(1,3), 'MarkerEdgeColor', 'g');
GPS_A4path  = scatter3(GPS_A4(1,1), GPS_A4(1,2), GPS_A4(1,3), 'MarkerEdgeColor', 'k');

GPS_B1Fpath  = scatter3(GPS_B1F(1,1), GPS_B1F(1,2), GPS_B1F(1,3), 'MarkerEdgeColor', 'r');
GPS_B1Apath  = scatter3(GPS_B1A(1,1), GPS_B1A(1,2), GPS_B1A(1,3), 'MarkerEdgeColor', 'b');
GPS_B2path   = scatter3(GPS_B2(1,1), GPS_B2(1,2), GPS_B2(1,3), 'MarkerEdgeColor', 'g');
GPS_B3path   = scatter3(GPS_B3(1,1), GPS_B3(1,2), GPS_B3(1,3), 'MarkerEdgeColor', 'k');
GPS_B4path   = scatter3(GPS_B4(1,1), GPS_B4(1,2), GPS_B4(1,3), 'MarkerEdgeColor', 'c');
    equat_rad=R_Earth;
    polar_rad=6356.7523142;
    [xx yy zz]=ellipsoid (0,0,0,equat_rad, equat_rad, polar_rad);
    load('topo.mat','topo','topomap1');
    topo2 = [topo(:,181:360) topo(:,1:180)];
    pro.FaceColor= 'texture';
    pro.EdgeColor = 'none';
    pro.FaceLighting = 'phong';     
    pro.Cdata = topo2;
    earth= surface(xx,yy,zz,pro);
    colormap(topomap1) 
    %rotate (earth, [0 0 1], GMST_deg(i));

for i = 1:1000:length(time)
    delete(GPS_A1path); delete(GPS_A2path); delete(GPS_A3path); delete(GPS_A4path);
    delete(GPS_B1Fpath); delete(GPS_B1Apath); delete(GPS_B2path); 
    delete(GPS_B3path); delete(GPS_B4path);
    %{
    addpoints(curve1, GPS_A1(i,1), GPS_A1(i,2), GPS_A1(i,3));
    addpoints(curve2, GPS_A2(i,1), GPS_A2(i,2), GPS_A2(i,3));
    addpoints(curve3, GPS_A3(i,1), GPS_A3(i,2), GPS_A3(i,3));
    addpoints(curve4, GPS_A4(i,1), GPS_A4(i,2), GPS_A4(i,3));
    addpoints(curve5, GPS_B1F(i,1), GPS_B1F(i,2), GPS_B1F(i,3));
    %}
    GPS_A1path  = scatter3(GPS_A1(i,1), GPS_A1(i,2), GPS_A1(i,3), 'filled', 'r');
    GPS_A2path  = scatter3(GPS_A2(i,1), GPS_A2(i,2), GPS_A2(i,3), 'filled', 'b');
    GPS_A3path  = scatter3(GPS_A3(i,1), GPS_A3(i,2), GPS_A3(i,3), 'filled', 'g');
    GPS_A4path  = scatter3(GPS_A4(i,1), GPS_A4(i,2), GPS_A4(i,3), 'filled', 'k');
    GPS_B1Fpath  = scatter3(GPS_B1F(i,1), GPS_B1F(i,2), GPS_B1F(i,3), 'filled', 'r');
    GPS_B1Apath  = scatter3(GPS_B1A(i,1), GPS_B1A(i,2), GPS_B1A(i,3), 'filled', 'b');
    GPS_B2path  = scatter3(GPS_B2(i,1), GPS_B2(i,2), GPS_B2(i,3), 'filled', 'g');
    GPS_B3path  = scatter3(GPS_B3(i,1), GPS_B3(i,2), GPS_B3(i,3), 'filled', 'k');
    GPS_B4path  = scatter3(GPS_B4(i,1), GPS_B4(i,2), GPS_B4(i,3), 'filled', 'c');
    
    drawnow
    delete(GPS_A1path); delete(GPS_A2path); delete(GPS_A3path); delete(GPS_A4path);
    delete(GPS_B1Fpath); delete(GPS_B1Apath); delete(GPS_B2path); 
    delete(GPS_B3path); delete(GPS_B4path);
end
%}

    
  