%% FFNAV HCW EKF RunSim ===================================================
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
%   CW_EKF_Sim.slx
%           Simulink simulation of the EKF
%   FFNAV_CW_EKF_Params.m
%           Initialization of the parameters for the EKF

% Created by: Cory Fraser - October 2015
% Latest Edit: Cory Fraser - July 14, 2017
%% ========================================================================
%Initialization
fprintf('\n--------------------------------------------------------------\n')
fprintf('         FFNAV EKF (CW RELATIVE DYNAMICS) SIMULATION          \n')
fprintf('--------------------------------------------------------------\n')
%Load initial conditions file
if isempty(dir('*Params.mat')) 
    fprintf('Initialization parameters were not found in this directory...\n')
    pause(2)
    fprintf('Generating initial parameters... \n')
    FFNAV_HCW_EKF_MakeParams;
    pause(2)
end
file_in = 'FFNAV_HCW_PEO_Params.mat';
load(file_in)
fprintf('\nParameter file %s has been loaded...\n', file_in)

 %% ======================================================================
% EKF Simulation Setup
time_step       = 1;        % Time step for the simulation
time_end        = 30000;    % Duration of the simulation
save_output     = 'on';     % 'on' to save output data to a .mat file

% Open the Simulink Model, set simulation parameters
open_system('HCW_EKF_Sim');
set_param('HCW_EKF_Sim', 'Solver', 'ode4')
set_param('HCW_EKF_Sim', 'SolverType', 'Fixed-step')
set_param('HCW_EKF_Sim', 'FixedStep', 'time_step')
set_param('HCW_EKF_Sim', 'RelTol','1e-12')
set_param('HCW_EKF_Sim', 'AbsTol','1e-12')
set_param('HCW_EKF_Sim', 'StartTime', '0')
set_param('HCW_EKF_Sim', 'StopTime', 'time_end')

%Select Adaptation Mode (1= ON, 0 = OFF)
Q_adapt = '0'; 
R_adapt = '0'; 

set_param('EKF_Sim/Q Adaptation/Switch: Q Adaptation On//Off', 'sw', Q_adapt)
set_param('EKF_Sim/R Adaptation/Switch: R Adaptation On//Off', 'sw', R_adapt)

% Run the Simulation
fprintf('\nStarting the EKF simulation...\n')
tic;
sim('HCW_EKF_Sim');
toc;
fprintf('EKF simulation complete.\n')
fprintf('--------------------------------------------------------------\n')

%% ========================================================================
% Cleaning up and saving the data
clearvars -except   time save_output file_in...
                    r_clean v_clean theta_target ...
                    EKF_output HCW_output GPS_output;

% Create .mat file of output data
if strcmp(save_output, 'on')                 
    file_out = strcat(file_in(1:end-10), 'Results.mat');
    clear save_output
    save(file_out)
    FFNAV_HCW_EKF_PostProcess(file_out)
end