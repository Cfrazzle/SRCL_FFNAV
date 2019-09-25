%% FFNAV_FAEKF_RunSim =====================================================
% FFNAV Fuzzy Adaptive Extended Kalman Filter Implementation ==============
% Description: This script executes a fuzzy adaptive extended kalman filter
% for relative navigation of a spacecraft formation. The exact nonlinear
% equations of relative motion are utilized in the dynamics model, and GPS
% measurements are simulated using an orbit propogator with various
% perturbations and noise sources.
%
% Inputs:
%   Parameter_filename  - Spacecraft and formation parameters filename
%   EKF_flag            - Select filter (EKF, MLE-EKF, or FAEKF)
%   Q_adapt             - Choose Q-adaptations (ON/OFF)
%   R_adapt             - Choose R-adaptations (ON/OFF)
%   Smooth_flag         - Choose residual fixed-window smoothing (ON/OFF)
%   N_window            - Fixed-window size for smoothing
%   PropOptions.        - Structure for orbit propagator settings
%   time_start          - Start time for the simulation
%   time_step           - Time step of the simulation
%   orbit_num           - Number of orbits to simulate
%   time_end            - End time of the simulation
%
%   save_flag           - Choose to save output to a .mat file
%   post_flag           - Choose to post-process the data
%   n_start             - Starting point for post-processing
%   n_end               - Ending point for post-processing
%   plot_flag           - Choose to create output plots
%   print_flag          - Choose to print plots to .eps files
%   profiler_flag       - Choose to show performance profiler
%
% Outputs:
%   Results.mat file
%   Figures & animations (optional)
%   RMS Error Analysis (optional)
%
% Other Functions Called:
%   FFNAV_EKF_MakeParams.m      - Create FAEKF parameters (if none found)
%   OrbitalMechanicsConstants   - Initialize astrodynamics constants
%   EKF_Sim_FAEKF.slx           - Simulink simulation of the FAEKF
%
% Created by:  Cory Fraser - OCT 31, 2017
% Latest Edit: Cory Fraser - OCT 21, 2018
% Copyright(c) 2018 by Cory Fraser
% =========================================================================

%% Initialization
clear
close all

fprintf('\n--------------------------------------------------------------\n')
fprintf(' \t\t FFNAV EXTENDED KALMAN FILTER SIMULATION \t\t             \n')
fprintf('--------------------------------------------------------------\n')

% Define path to parameters/output data
pathname = fileparts('C:\Users\Cory\Desktop\FFNAV Data\');

% Check for parameter files
if isempty(dir([pathname, '\*Params.mat']))
    fprintf('Initialization parameters were not found in this directory...\n')
    pause(2)
    fprintf('Generating initial parameters... \n')
    
    mydir  = pwd;                       %Get current directory
    idcs   = strfind(mydir,filesep);    %Index file separations
    newdir = mydir(1:idcs(end)-1);      %Path to higher folder (1-level up)
    addpath(newdir)                     %Add path to location of MakeParams
    FFNAV_EKF_MakeParams;               %Create a parameter file
end

% =========================================================================
%% Select Formation Configuration Parameter File

Parameter_filename = 'PRISMA';

%Thesis Test Cases
    % PRISMA     = PRISMA Mission (LEO, low eccentricty)
    % PROBA-3    = PROBA-3 Mission (HEO, high eccentricity)
    % PEOinLEO   = Modified PRISMA, with inclination change (LEO,low e,delta i)
%Other Test Cases
    % PROBA-3Mod = PROBA-3 Mission (HEO, high eccentricity)
    % Kuiack33   = Figure 3.3  in Thesis (In-plane Elliptical Orbit, low e)
    % Kuiack34   = Figure 3.4  in Thesis (In-plane Elliptical Orbit, e=0.2)
    % Kuiack35   = Figure 3.5  in Thesis (In-plane Elliptical Orbit, e=0.4)
    % Kuiack511  = Figure 5.11 in Thesis (Highly Elliptical Orbit, e=0.8)
    % Kuiack_J2invariant = Figure 5.15 in Thesis (Elliptical Orbit, e=0.8)
    % LF_Hold    = Leader-follower Along-track holding position
    % BUSSE      = Busse Article   (LEO)
    
file_in = fullfile(pathname, ['FFNAV_', Parameter_filename, '_Params.mat']);
load(file_in)
fprintf('\nParameter file loaded : %s  \n', [Parameter_filename, '_Params.mat'])

% =========================================================================
%% Select Orbit Propagator Options (1 = ON, 0 = OFF)
PropOptions.J2              = 1;
PropOptions.Sun             = 1;        %Assuming fixed position of Sun
PropOptions.Moon            = 1;        %Assuming fixed position of Moon
PropOptions.Planets         = 0;        %Assuming fixed positions of planets
PropOptions.Drag            = 1;        %Using Harris-Priester Model  
PropOptions.SolarRad        = 1;        
    PropOptions.ShadowModel = 0;        %0=Cylindrical Shadow, 1=Geometric 
PropOptions.Relativity      = 0;        
PropOptions.GravEl_Earth    = 0;        %NOT COMPLETE
    
%Add path to orbit propagator, load orbital mechanics constants
addpath('C:\Users\Cory\OneDrive\FFNAV Project\MATLAB Code\Orbit Propagator')
OrbitalMechanicsConstants
    
fprintf('Astrodynamic constants: %s \n', ConstantsModel)
fprintf('Orbital Perturbations : ')
if (PropOptions.J2 || PropOptions.GravEl_Earth || PropOptions.Sun || ...
        PropOptions.Moon ||  PropOptions.Planets || PropOptions.Relativity...
        || PropOptions.Drag || PropOptions.SolarRad )
    fprintf('On \n')
else
    fprintf('None (Two-Body Motion) \n')
end
if (PropOptions.J2)
    fprintf('    - Earth''s Oblateness (J2) \n')
end
if (PropOptions.Sun)
    fprintf('    - Solar Third-Body Gravity\n')
end
if (PropOptions.Moon)
    fprintf('    - Lunar Third-Body Gravity \n')
end
if (PropOptions.Planets)
    fprintf('    - Planetary Third-Body Gravity \n')
end    
if (PropOptions.SolarRad) 
    fprintf('    - Solar Radiation Pressure \n')    
end
if (PropOptions.Drag)
    fprintf('    - Atmospheric Drag \n')
end
if (PropOptions.Relativity)
    fprintf('    - Relativistic Effects \n')
end
if (PropOptions.GravEl_Earth)
    fprintf('    - Harmonic Gravity (Elastic Earth) \n')
end

% Create PropOptions for target spacecraft (same for both spacecraft)
PropOptions_target.J2              = PropOptions.J2;
PropOptions_target.Sun             = PropOptions.Sun;
PropOptions_target.Moon            = PropOptions.Moon;
PropOptions_target.Planets         = PropOptions.Planets;
PropOptions_target.Drag            = PropOptions.Drag;        
PropOptions_target.SolarRad        = PropOptions.SolarRad;        
    PropOptions_target.ShadowModel = PropOptions.ShadowModel; 
PropOptions_target.Relativity      = PropOptions.Relativity;        
PropOptions_target.GravEl_Earth    = PropOptions.GravEl_Earth;

% Create PropOptions for chaser spacecraft (same for both spacecraft)
PropOptions_chaser.J2              = PropOptions.J2;
PropOptions_chaser.Sun             = PropOptions.Sun;
PropOptions_chaser.Moon            = PropOptions.Moon;
PropOptions_chaser.Planets         = PropOptions.Planets;
PropOptions_chaser.Drag            = PropOptions.Drag;        
PropOptions_chaser.SolarRad        = PropOptions.SolarRad;        
    PropOptions_chaser.ShadowModel = PropOptions.ShadowModel; 
PropOptions_chaser.Relativity      = PropOptions.Relativity;        
PropOptions_chaser.GravEl_Earth    = PropOptions.GravEl_Earth;

% =========================================================================
%% Select Kalman Filter Options

% 0 = EKF           Extended Kalman Filter        
% 1 = MLE-EKF       MLE Adaptive Extended Kalman Filter
% 2 = FAEKF         Fuzzy Adaptive Extended Kalman Filter 
EKF_flag = 0;

% Select Adaptation Options - MLE-AEKF & FAEKF (1 = ON, 0 = OFF)
Q_adapt = 0; 
R_adapt = 0; 

% Select Residuals Smoothing Options - MLE-AEKF & FAEKF (1 = ON, 0 = OFF)
Smooth_flag = '1';
N_window    = 30;  %Running-Average Window Size 
%N=5 seems to work better for MLE-AEKF

% =========================================================================
%% Select Simulation Options
time_start   = 0;                           % Start time of the simulation
time_step    = 1;                           % Time step for the simulation
orbit_num    = 5;                           % Number of orbits to simulate
time_end     = ceil(orbit_num*T_target);    % End time of the simulation

% Post Processing Options
save_flag     = 'on';     % Save data to a .mat file
post_flag     = 'on';     % Post-process the simulation data
    n_start   = floor(time_end*0.5); 
    n_end     = time_end+1;  
    cov_flag  = 'on';    % Perform analysis of covariances
plot_flag     = 'on';     % Create output plots
print_flag    = 'off';    % Print plots to .eps files
profiler_flag = 'off';    % Show profile simulink performance

%PEOinLEO Thesis Case
if strcmp(FF_Config, 'PEOinLEO')
    Q0 = Q0;
    R0 = 5*R0;
    GPS_r_STD = 10*GPS_r_STD; %(m)  
    GPS_v_STD = 10*GPS_v_STD; %(m/s)

end

%PROBA-3 Thesis Case
if strcmp(FF_Config, 'PROBA-3')
     Q0 = 10*Q0;
     R0 = R0;
end

% ========================================================================
% Open the Simulink Model, set simulation parameters

Sim_filename = 'EKF_Sim_FAEKF';

load_system(Sim_filename)
if strcmp(get_param(Sim_filename, 'shown'), 'off')
    open_system(Sim_filename);
end
set_param(Sim_filename, 'Solver', 'ode4')
set_param(Sim_filename, 'SolverType', 'Fixed-step')
set_param(Sim_filename, 'FixedStep', 'time_step')
set_param(Sim_filename, 'RelTol','1e-6')
set_param(Sim_filename, 'AbsTol','1e-6')
set_param(Sim_filename, 'StartTime', 'time_start')
set_param(Sim_filename, 'StopTime', 'time_end')
set_param(Sim_filename, 'Profile', profiler_flag)
fprintf('\nSimulink model loaded : %s.slx \n', Sim_filename)

% =========================================================================
%% Run the EKF Simulation

switch EKF_flag
    
    %Extended Kalman Filter
    case 0
        EKF_name = 'EKF';
        fprintf('Filtering algorithm   : EKF \n')
        File_specifier = 'NOadapt_5orbits.mat';

    %MLE Adaptive Extended Kalman Filter
    case 1 
        EKF_name = 'MLE_EKF';
        fprintf('Filtering algorithm   : MLE Adaptive EKF \n')   
        fprintf('RTS Smoother          : On, N = %i \n', N_window)

        %Set Q and R Adaptation Switches
        if (Q_adapt) 
            fprintf('Q-adaptations         : On \n')
        else
            fprintf('Q-adaptations         : Off \n')
        end
        if (R_adapt)
            fprintf('R-adaptations         : On \n')
        else
            fprintf('R-adaptations         : Off \n')
        end
             
    %Fuzzy Adaptive Extended Kalman Filter
    case 2 
        EKF_name = 'FAEKF';
        fprintf('Filtering algorithm   : Fuzzy Adaptive EKF \n')
        
        % Initialize Fuzzy Inference System
        fprintf('Fuzzy system loaded   :')
        MSFplot_flag = 0; % 1 = Plot Membership Functions
        FUZZY_props = FFNAV_FAEKF_FuzzyInitialize(MSFplot_flag);
        
        %Set Smoothing Settings
        if Smooth_flag == '1'
            fprintf('Residual smoothing    : On, N = %i \n', N_window)
        else
            fprintf('Residual smoothing    : Off \n')
        end
        set_param('EKF_Sim_FAEKF/Kalman Filter/FAEKF/Adaptations/Switch: Smoothing On//Off', 'sw', Smooth_flag)
        
        %Set Q and R Adaptation Switches
        if (Q_adapt)
            fprintf('Q-adaptations         : On \n')
        else
            fprintf('Q-adaptations         : Off \n')
        end
        if (R_adapt)
            fprintf('R-adaptations         : On \n')
        else
            fprintf('R-adaptations         : Off \n')
        end
end

%Set output filename based on adaptation scheme
if (EKF_flag~=0) && (Q_adapt) && (R_adapt)
    File_specifier = 'QRadapt.mat';
elseif (EKF_flag~=0) && (Q_adapt)
    File_specifier = 'Qadapt.mat';
elseif (EKF_flag~=0) && (R_adapt)
    File_specifier = 'Radapt.mat';
else
    File_specifier = 'NOadapt.mat';
end
% Start the simulation
fprintf('Simulation duration   : Orbital period x %i \n', orbit_num)
fprintf('                      : %.0f seconds \n', time_end)
fprintf('.............................................\n')
fprintf('....... running the simulation ..............\n')
fprintf('.............................................\n')

if strcmp(profiler_flag,'on')
    profile on
end

tic;
sim(Sim_filename);
t_run = toc;

if strcmp(profiler_flag,'on')
    profile off
end

fprintf('EKF simulation complete \n')
fprintf('Runtime = %f \n', t_run)
fprintf('..............\\(-_-)/........................\n')

% =========================================================================
%% Clean-up and save data
clearvars -except   EKF_output   GPS_output     NERM_output             ...
                    EKF_flag     EKF_name       time        t_run       ...
                    r_clean      v_clean        R_Earth     mu          ...
                    r_target_ECI v_target_ECI                           ...
                    r_chaser_ECI v_chaser_ECI                           ...
                    theta_target rt_clean       rt_dot_clean            ...
                    theta_clean  theta_dot_clean                        ...
                    Pv_obs       Pv_smoothed    Pr_theo                 ...
                    P_error      Pdot_error     P0                      ...
                    Q_output     R_output       Q0          R0          ...
                    lambda       NIS            deltaDOD    DOD_true    ...
                    save_flag    post_flag      plot_flag   print_flag  ...
                    pathname     Parameter_filename         Sim_filename...
                    File_specifier Constants...
                    profiler_flag  cov_flag Q_adapt R_adapt n_start n_end;
                    
if strcmp(save_flag, 'on') %(1:end-10) for line below
    file_out = strcat(Parameter_filename, ['_', EKF_name,...
                            '_',File_specifier]);
    clear save_output
    file_out = fullfile(pathname, file_out);
    save(file_out);        
end

% =========================================================================
%% Post Processing  
if strcmp(post_flag, 'on')
    tic;
    FFNAV_FAEKF_PostProcess(file_out)
    t_process = toc;
    fprintf('Post-processing Time = %f \n', t_process)
end

load gong
sound(y,Fs) 

% Plotting
if strcmp(plot_flag, 'on')           
    FFNAV_FAEKF_Plotter(file_out)
end

% Computation Profile Viewer
if strcmp(profiler_flag,'on')
    profile viewer
end

fprintf('--------------------------------------------------------------\n')
% =========================================================================