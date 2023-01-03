%% FFNAV_EKF_RunSim ================================================================
% FFNAV Extended Kalman Filter Implementation =============================
% Description: This script performs extended Kalman filtering on a noisy 
% set of GPS data using the exact non-linear equations of relative motion. 
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
%   EKF_Sim.slx
%           Simulink simulation of the EKF
% Created by: Cory Fraser - October 2015
% Latest Edit: Cory Fraser - August 8, 2017
%% ========================================================================
%Initialization
clear
%clc
fprintf('\n--------------------------------------------------------------\n')
fprintf('   FFNAV EKF (NONLINEAR RELATIVE DYNAMICS) SIMULATION         \n')
fprintf('--------------------------------------------------------------\n')
%Load initial conditions file

%top_level_directory = 'C:\Users\Cory\Desktop\FFNAV Data\';
top_level_directory = 'C:\Users\CF4715\Desktop\FFNAV Data\';

if isempty(dir([top_level_directory '*Params.mat']))
    fprintf('Initialization parameters were not found in this directory...\n')
    pause(2)
    fprintf('Generating initial parameters... \n')
    FFNAV_EKF_MakeParams(top_level_directory);
    pause(2)
end
pathname = fileparts(top_level_directory);
filename = 'FFNAV_PEO_Params.mat';
file_in = fullfile(pathname, filename);

load(file_in)
fprintf('\nParameter file %s has been loaded...\n', filename)

 %% ======================================================================
% EKF Simulation Setup
time_step   = 1;        % Time step for the simulation
time_end    = 13000;  % Duration of the simulation (~Two Complete Orbits)
save_output = 'on';     % 'on' to save output data to a .mat file

% Open the Simulink Model, set simulation parameters
load_system('EKF_Sim')
if strcmp(get_param('EKF_Sim', 'shown'), 'off')
    open_system('EKF_Sim');
end
set_param('EKF_Sim', 'Solver', 'ode4')
set_param('EKF_Sim', 'SolverType', 'Fixed-step')
set_param('EKF_Sim', 'FixedStep', 'time_step')
set_param('EKF_Sim', 'RelTol','1e-15')
set_param('EKF_Sim', 'AbsTol','1e-15')
set_param('EKF_Sim', 'StartTime', '0')
set_param('EKF_Sim', 'StopTime', 'time_end')

% Run the EKF Simulation
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
                    Pv_obs Pv_smoothed Pr_theo P_error Pdot_error;
                    
if strcmp(save_output, 'on')
    file_out = strcat(filename(1:end-10), 'Results.mat');
    clear save_output
    file_out = fullfile(pathname, file_out);
    save(file_out);
    FFNAV_EKF_PostProcess(file_out)
end