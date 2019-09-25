%% FFNAV_EKF_RunSim ================================================================
% FFNAV Extended Kalman Filter Implementation =============================
% Description: This script performs extended Kalman filtering on a noisy 
% set of GPS data using the exact nonlinear equations of relative motion. 
%
% Inputs:
%   .mat Parameter File - Contain initial conditions for EKF
%   time_step - Defines size of the time step
%   time_end - Defines duration of the simulation
%   make_plots - Turns plot output on/off
%   print_plots - Turns .eps plot printing on/off
%   show_RMS - Turns RMS claculations on/off
%
% Outputs:
%   .mat Results file
%   Result plots (optional)
%   RMS Error Analysis (optional)
%
% Other Functions Called:
%   FFNAV_EKF_Params.m
%           Initialization of the parameters for the EKF
%   EKF_Sim_MLE.slx
%           Simulink simulation of the EKF

% Created by: Cory Fraser - October 2015
% Latest Edit: Cory Fraser - September 20, 2017
%% ========================================================================
%Initialization
clear
%clc
fprintf('\n--------------------------------------------------------------\n')
fprintf('   FFNAV EKF (NONLINEAR RELATIVE DYNAMICS) SIMULATION         \n')
fprintf('--------------------------------------------------------------\n')

%Load initial conditions file
if isempty(dir('C:\Users\Cory\Desktop\FFNAV Data\*Params.mat'))
    fprintf('Initialization parameters were not found in this directory...\n')
    pause(2)
    fprintf('Generating initial parameters... \n')
    FFNAV_EKF_MakeParams;
    pause(2)
end
pathname = fileparts('C:\Users\Cory\Desktop\FFNAV Data\');
filename = 'FFNAV_PEO_Params.mat';
file_in = fullfile(pathname, filename);

load(file_in)
fprintf('\nParameter file %s has been loaded...\n', filename)

 %% ======================================================================
% MLE EKF Simulation Setup
time_step   = 1;        % Time step for the simulation
time_end    = 13000;  % Duration of the simulation (~Two Complete Orbits)
save_output = 'on';     % 'on' to save output data to a .mat file

% Open the Simulink Model, set simulation parameters
load_system('EKF_Sim_MLE')
if strcmp(get_param('EKF_Sim_MLE', 'shown'), 'off')
    open_system('EKF_Sim_MLE');
end
set_param('EKF_Sim_MLE', 'Solver', 'ode4')
set_param('EKF_Sim_MLE', 'SolverType', 'Fixed-step')
set_param('EKF_Sim_MLE', 'FixedStep', 'time_step')
set_param('EKF_Sim_MLE', 'RelTol','1e-15')
set_param('EKF_Sim_MLE', 'AbsTol','1e-15')
set_param('EKF_Sim_MLE', 'StartTime', '0')
set_param('EKF_Sim_MLE', 'StopTime', 'time_end')

% Select Adaptation Options (1= ON, 0 = OFF)
Q_adapt     = '1'; 
R_adapt     = '1'; 
%Smoothing   = '1';
N_window    = 30; 

set_param('EKF_Sim_MLE/Q Selection/Switch: Q Adaptation On//Off', 'sw', Q_adapt)
set_param('EKF_Sim_MLE/R Selection/Switch: R Adaptation On//Off', 'sw', R_adapt)

% Run the MLE EKF Simulation
fprintf('\nStarting the MLE-EKF simulation...\n')
tic;
sim('EKF_Sim_MLE');
runtime = toc;
fprintf('MLE-EKF simulation complete.\n')
fprintf('--------------------------------------------------------------\n')

%Store_output = squeeze(Storage_output)';
% ========================================================================
%Cleaning up and saving the data
clearvars -except   EKF_output time runtime Store_output...
                    save_output file_in filename pathname...
                    r_clean v_clean...
                    GPS_output ...
                    theta_target ...
                    rt_clean rt_dot_clean...
                    theta_clean theta_dot_clean...
                    NERM_output...
                    Pv_obs Pv_smoothed Pr_theo P_error Pdot_error...
                    Q_output R_output Q0 R0...
                    L;
                    
if strcmp(save_output, 'on')
    file_out = strcat(filename(1:end-10), 'NoAdaptResults_Temp.mat');
    clear save_output
    file_out = fullfile(pathname, file_out);
    save(file_out);
    FFNAV_EKF_PostProcess(file_out)
end

%pause
%{
%% ========================================================================
% EKF Simulation Setup
pathname = fileparts('C:\Users\Cory\Desktop\FFNAV Data\');
filename = 'FFNAV_PEO_Params.mat';
file_in = fullfile(pathname, filename);

load(file_in)
fprintf('\nParameter file %s has been loaded...\n', filename)


time_step   = 1;        % Time step for the simulation
time_end    = 13000;  % Duration of the simulation (~Two Complete Orbits)
save_output = 'on';     % 'on' to save output data to a .mat file

% Run the EKF Simulation 
path(path,'C:\Users\Cory\OneDrive\FFNAV Project\MATLAB Code\NERM EKF Simulation\EKF')
load_system('EKF_Sim')
if strcmp(get_param('EKF_Sim', 'shown'), 'off')
    open_system('EKF_Sim');
end
fprintf('\nStarting the EKF simulation...\n')
tic;
sim('EKF_Sim');
toc;
fprintf('EKF simulation complete.\n')
fprintf('--------------------------------------------------------------\n')

%% ========================================================================
%Cleaning up and saving the data
clearvars -except   EKF_output time save_output file_in filename pathname...
                    r_clean v_clean...
                    GPS_output ...
                    theta_target ...
                    rt_clean rt_dot_clean...
                    theta_clean theta_dot_clean...
                    NERM_output...
                    Pv_obs Pv_smoothed Pr_theo P_error Pdot_error...
                    Q_output R_output Q0 R0...
                    L;
                    
if strcmp(save_output, 'on')
    file_out = strcat(filename(1:end-10), 'EKF_Results.mat');
    clear save_output
    file_out = fullfile(pathname, file_out);
    save(file_out);
    FFNAV_EKF_PostProcess(file_out)
end
%}