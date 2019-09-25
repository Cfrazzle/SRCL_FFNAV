%% High Fidelity Orbit Propagator =========================================
% =========================================================================
% Description: This script simulates an orbit with various perturbations.
%
% Inputs:
%   
% Created by: Cory Fraser - JUN 20, 2018
% Last Edit : Cory Fraser - JUN 30, 2018
% Copyright(c) 2018 by Cory Fraser
% =========================================================================

%% Initialization
clear
clc
close all

fprintf('\n--------------------------------------------------------------\n')
fprintf(' \t\t\t\t FFNAV ORBIT PROPAGATOR \t\t             \n')
fprintf('--------------------------------------------------------------\n')

% =========================================================================
%% General Constants

%Add path to orbit propagator, load orbital mechanics constants
addpath('C:\Users\Cory\OneDrive\FFNAV Project\MATLAB Code\Orbit Propagator')
OrbitalMechanicsConstants

%Add path to formation parameters, load formation configuration
pathname = fileparts('C:\Users\Cory\Desktop\FFNAV Data\');
Parameter_filename = 'PRISMA_Mod';

    % PRISMA_Mod = Modified PRISMA (LEO, low eccentricty)
    % BUSSE      = Busse Article   (LEO, 
    % PROBA-3    = PROBA-3 Mission (LEO, high eccentricity)
    % PROBA-3Mod = PROBA-3 Mission (LEO, high eccentricity)
    % Kuiack33   = Figure 3.3 in Thesis (In-plane Elliptical Orbit, low e)
    % Kuiack34   = Figure 3.4 in Thesis (In-plane Elliptical Orbit, e=0.2)
    % Kuiack35   = Figure 3.5 in Thesis (In-plane Elliptical Orbit, e=0.4)
        %Current Results - QR-FAEKF>QR-MLE-AEKF>EKF
    % Kuiack511   = Figure 5.11 in Thesis (Highly Elliptical Orbit, e=0.8)
    % Kuiack_J2invariant = Figure 5.15 in Thesis ( Elliptical Orbit, e=0.8)
    % PRISMA_Mod_Scaled
    % PRISMA_ModTemp = %Includes 1 SDPR
    % PEO = Early version of projected-elliptical orbit
    
file_in = fullfile(pathname, ['FFNAV_', Parameter_filename, '_Params.mat']);
load(file_in)
fprintf('\nParameter file loaded : %s  \n', Parameter_filename)

% =========================================================================
%% Orbit Simulation Setup

time_step   = 1;                   % Time step for the simulation
orbit_num   = 2;                   % Number of orbits to simulate
time_end    = orbit_num*T_target;  % Duration of the simulation
Kepler_comp = 1;                   % Compare Purturbed w/ Keplerian Motion

PropOptions.J2              = 1;
PropOptions.GravEl_Earth    = 0; %NOT COMPLETE
PropOptions.Sun             = 1; %Assuming fixed position
PropOptions.Moon            = 1; %Assuming fixed position
PropOptions.Planets         = 1; %Assuming fixed positions
PropOptions.Relativity      = 1;

PropOptions.Drag            = 1;     %Harris-Priester Model
    PropOptions.Cd          = 2;     %Coefficient of Drag
    PropOptions.A_Drag      = 110.5; %Area exposed to drag forces
    
PropOptions.SolarRad        = 1;    
    PropOptions.A_Solar     = 100;     %Area exposed to solar radiation
    PropOptions.M_Sat       = 8000;    %Spacecraft mass
    PropOptions.Cr          = 2;       %Coefficient of Reflection
    PropOptions.ShadowModel = 0;    %0 = Cylindrical Shadow, 1 = Geometric 


fprintf('Astrodynamic constants: %s \n', ConstantsModel)
fprintf('Orbital Perturbations : ')
if (PropOptions.J2 || PropOptions.GravEl_Earth || PropOptions.Sun || ...
        PropOptions.Moon ||  PropOptions.Planets || PropOptions.Relativity...
        || PropOptions.Drag || PropOptions.SolarRad )
    fprintf('Activated \n')
else
    fprintf('None (Two-Body Motion) \n')
end
if (PropOptions.J2)
    fprintf('\t    -Earth''s Oblateness (J2) \n')
end
if (PropOptions.Sun)
    fprintf('\t    -Solar Third-Body \n')
end
if (PropOptions.Moon)
    fprintf('\t    -Lunar Third-Body \n')
end
if (PropOptions.Planets)
    fprintf('\t    -Planetary Third-Body \n')
end    
if (PropOptions.SolarRad) 
    fprintf('\t    -Solar Radiation Pressure \n')    
end
if (PropOptions.Drag)
    fprintf('\t    -Atmospheric Drag \n')
end
if (PropOptions.Relativity)
    fprintf('\t    -Relativistic Effects \n')
end
if (PropOptions.GravEl_Earth)
    fprintf('\t    -Harmonic Gravity (Elastic Earth) \n')
end

% Time Definitions
year = 2002;
month = 04;
day = 24;
hour = 21;
minute = 55;
second = 28;

[JD_start,~] = JulianDate(day, month, year, hour, minute, second);
MJD_start = JD_start - 2400000.5; 
PropOptions.JD_UTC_start = JD_start;
PropOptions.MJD_UTC_start = MJD_start;

%Load planetary coefficients relevant from Julian Date
%idx = find(PC(:,1)<=JD_start & JD_start<=PC(:,2),1); 
%Constants.PC_idx = PC(idx,:);

%Evaluate positions of the planets, Sun and Moon
%{
[r_Mercury,r_Venus,r_Earth,r_Mars,r_Jupiter,r_Saturn,r_Uranus, ...
 r_Neptune,r_Pluto,r_Moon,r_Sun,r_SunSSB] = ...
    PlanetPositions(MJD_start, JD_start, Constants);

%Assume fixed Planets, positions evaluated at simulation start time
PropOptions.r_Mercury = r_Mercury;
PropOptions.r_Venus = r_Venus;
PropOptions.r_Earth = r_Earth;
PropOptions.r_Mars = r_Mars;
PropOptions.r_Jupiter = r_Jupiter;
PropOptions.r_Saturn = r_Saturn;
PropOptions.r_Uranus = r_Uranus;
PropOptions.r_Neptune = r_Neptune;
PropOptions.r_Pluto = r_Pluto; 
PropOptions.r_Moon = r_Moon;
PropOptions.r_Sun = r_Sun;
PropOptions.r_SunSSB = r_SunSSB;
%}
% =========================================================================
%% Open the Simulink Model, set simulation parameters
load_system('Orbit_Propagator')
if strcmp(get_param('Orbit_Propagator', 'shown'), 'off')
    open_system('Orbit_Propagator');
end
set_param('Orbit_Propagator', 'Solver', 'ode4')
set_param('Orbit_Propagator', 'SolverType', 'Fixed-step')
set_param('Orbit_Propagator', 'FixedStep', 'time_step')
set_param('Orbit_Propagator', 'RelTol','1e-9')
set_param('Orbit_Propagator', 'AbsTol','1e-9')
set_param('Orbit_Propagator', 'StartTime', '0')
set_param('Orbit_Propagator', 'StopTime', 'time_end')

% Run the EKF Simulation
fprintf('Simulation duration   : Orbital period x %i \n', orbit_num)
fprintf('                      : %.0f seconds \n', time_end)
fprintf('\nStarting the simulation...\n')
tic;
sim('Orbit_Propagator');
toc;
fprintf('Orbit simulation complete.\n')
fprintf('--------------------------------------------------------------\n')

% ========================================================================
%% Plots for EKF ==========================================================

set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', 20)
set(0,'DefaultTextFontName', 'Times New Roman')
set(0,'DefaultTextFontSize', 24)
set(0,'defaultLineLineWidth', 2)

time = tout/3600;   %Simulation time from seconds to hours, for plots

%% Figure: 3D Position Plot for GPS
fig_id = 1;
figure(fig_id); fig_id = fig_id + 1;    
set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
grid on
hold on
plot3(r_target_ECI(:,1),r_target_ECI(:,2),r_target_ECI(:,3), '--r');
view(150,10)
%plot3(r_chaser_ECI(:,1),r_chaser_ECI(:,2),r_chaser_ECI(:,3), ':b');
lim_ECI = 1.1*max(max(abs(r_target_ECI)));
title('ECI Frame')
set(gca, 'XLim', 1.1*[-lim_ECI lim_ECI], 'YLim', ...
    1.1*[-lim_ECI lim_ECI], 'ZLim', ...
    1.1*[-lim_ECI lim_ECI]);
legend('Target Orbit') 
set(gca,'xticklabel',num2str(get(gca,'xtick')'/1000))
set(gca,'yticklabel',num2str(get(gca,'ytick')'/1000))
set(gca,'zticklabel',num2str(get(gca,'ztick')'/1000)) 

%Plot Earth
equat_rad = Constants.R_Earth;
polar_rad = Constants.R_Earth_p;
[xx yy zz]=ellipsoid (0,0,0,equat_rad, equat_rad, polar_rad);
load('topo.mat','topo','topomap1');
topo2 = [topo(:,181:360) topo(:,1:180)];
pro.FaceColor= 'texturemap';
pro.EdgeColor = 'none';
pro.FaceLighting = 'gouraud';     
pro.Cdata = topo2;
pro.SpecularStrength = 0.5;      % change the strength of the reflected light
earth = surface(xx,yy,zz,pro,'HandleVisibility','off');
colormap(topomap1) 
xlabel(' X (km) ')
ylabel(' Y (km) ')
zlabel(' Z (km) ')
zlimits = [min(topo(:)) max(topo(:))];
demcmap(zlimits);
light('Position',[0.75 0.75 1])     % add a light
axis square

line([0 lim_ECI],[0 0],[0 0],'Color', 'black', 'Marker','.','LineStyle','-','LineWidth', 2.5, 'HandleVisibility','off')
line([0 0],[0 lim_ECI],[0 0],'Color', 'black', 'Marker','.','LineStyle','-','LineWidth', 2.5, 'HandleVisibility','off')
line([0 0],[0 0],[0 lim_ECI],'Color', 'black', 'Marker','.','LineStyle','-','LineWidth', 2.5, 'HandleVisibility','off')
 

%% Supblots of ECI Positions
 figure(fig_id); fig_id = fig_id + 1;
    h(1) = subplot(3,1,1);
    plot(time, r_target_ECI(:,1)) 
    ylabel('X (m)')
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    grid on
    hold on
    xlim([0 inf])

    h(2) = subplot(3,1,2);
    plot(time, r_target_ECI(:,2)) 
    ylabel('Y (m)')
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    grid on
    hold on
    xlim([0 inf])

    h(3) = subplot(3,1,3);
    plot(time, r_target_ECI(:,3)) 
    ylabel('Z (m)')
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    grid on
    hold on
    xlim([0 inf])
    xlabel('Time (hrs)')
        
%% Figure: 3D LVLH Position Plot
figure(fig_id); fig_id = fig_id + 1;    
set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
grid on
hold on
view(-90,90)
p1 = plot3(r_clean(:,1),r_clean(:,2),r_clean(:,3), '-.','color',[0,0,0]);
p1.Color(4) = 0.75;
axis equal
lim_ECI = 1.1*max(max(abs(r_clean)));
title('LVLH Frame')
set(gca, 'XLim', 1.1*[-lim_ECI lim_ECI], 'YLim', ...
    1.1*[-lim_ECI lim_ECI], 'ZLim', ...
    1.1*[-lim_ECI lim_ECI]);
xlabel(' X (m) ')
ylabel(' Y (m) ')
zlabel(' Z (m) ')

%% Comparing Two-body Motion w/ Perturbed Motion
if (Kepler_comp)

PropOptions.J2              = 0;
PropOptions.GravEl_Earth    = 0; 
PropOptions.Sun             = 0; 
PropOptions.Moon            = 0; 
PropOptions.Planets         = 0; 
PropOptions.Drag            = 0; 
PropOptions.SolarRad        = 0;
PropOptions.Relativity      = 0;

    
    fprintf('\nStarting the Two-Body simulation...\n')
    tic;
    sim('Orbit_Propagator');
    toc;
    fprintf('Orbit simulation complete.\n')
    fprintf('--------------------------------------------------------------\n')
    
    figure(fig_id-1);
    p1 = plot3(r_clean(:,1),r_clean(:,2),r_clean(:,3), '-.r');
    p1.Color(4) = 1;
    legend('Perturbed Motion', 'Two-Body Motion')
end