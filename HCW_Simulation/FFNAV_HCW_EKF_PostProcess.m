function [] = FFNAV_HCW_EKF_PostProcess(filename)
%% FFNAV_EKF_PostProcess===============================================
% FFNAV EKF Data Processing and Analysis ==================================
% Description: This script takes a .mat results file created by the EKF
% simulation, and computes and plots results.
%
% Inputs:
%   .mat Parameter File - Contain initial conditions for EKF
%   make_plots - Turns plot output on/off
%   print_plots - Turns .eps plot printing on/off
%   show_RMS - Turns RMS claculations on/off
%
% Outputs:
%   Result plots (optional)
%   RMS Error Analysis (optional)
%
% Other Functions Called:
%   HCW_EKF_Plotter
%
% Created by: Cory Fraser - July 2016
% Latest Edit: Cory Fraser - July 11, 2017
%% ========================================================================
%Initialization
fprintf('\n--------------------------------------------------------------\n')
fprintf('   FFNAV EKF SIMULATION POST-PROCESSING                       \n')
fprintf('--------------------------------------------------------------\n')
load(filename)
%% ========================================================================
% Selection of operations to be performed
show_RMS        = 'on';     % 'on' prints results of RMS error analysis 
make_plots      = 'on';     % 'on' to create various output plots
print_plots     = 'off';    % 'on' saves plots to .eps files after simulation

%% ========================================================================
%RMS Analysis of Position and Velocity Errors
errors_x         = EKF_output(:,1)  - r_clean(:,1); 
errors_y         = EKF_output(:,2)  - r_clean(:,2); 
errors_z         = EKF_output(:,3)  - r_clean(:,3); 
errors_xdot      = EKF_output(:,4)  - v_clean(:,1); 
errors_ydot      = EKF_output(:,5)  - v_clean(:,2); 
errors_zdot      = EKF_output(:,6)  - v_clean(:,3); 


%Examination of the RMS error of the Nosiy GPS and EKF signals
if strcmp(show_RMS, 'on') 
    format long
    fprintf(' Root-Mean-Squared Errors for Noisy GPS and EKF Output [km, km/s] \n')
    fprintf(' Row: 1 = x, 2 = y, 3 = z, 4 = xdot, 5 = ydot, 6 = zdot  \n\n')
    fprintf(' HCW Dynamics                Noisy GPS               EKF Output \n')    
    
    GPS_error = zeros(size(EKF_output,1), 1);
    EKF_error = zeros(size(EKF_output,1), 1);
    HCW_error = zeros(size(EKF_output,1), 1);
    
    for k = 1:1:6
        if k < 4
            GPS_error   = GPS_output(:,k) - r_clean(:,k);
            EKF_error   = EKF_output(:,k) - r_clean(:,k);
            HCW_error  = HCW_output(:,k) - r_clean(:,k);
        else 
            GPS_error   = GPS_output(:,k) - v_clean(:,k-3);
            EKF_error   = EKF_output(:,k) - v_clean(:,k-3);
            HCW_error  = HCW_output(:,k) - v_clean(:,k-3);
        end
        EKF_RMSE(k) = sqrt(mean(EKF_error.^2));
        GPS_RMSE(k) = sqrt(mean(GPS_error.^2));
        HCW_RMSE(k) = sqrt(mean(HCW_error.^2));

        fprintf(' %18.12f     %18.12f     %18.12f \n', HCW_RMSE(k), GPS_RMSE(k), EKF_RMSE(k))
    end   
        
    GPS_RMSE_Avg = mean(GPS_RMSE);
    EKF_RMSE_Avg = mean(EKF_RMSE);
    HCW_RMSE_Avg = mean(HCW_RMSE);
    
    GPS_RMSE_Avg_pos = mean(GPS_RMSE(1:3));
    EKF_RMSE_Avg_pos = mean(EKF_RMSE(1:3));
    HCW_RMSE_Avg_pos = mean(HCW_RMSE(1:3));
    
    GPS_RMSE_Avg_vel = mean(GPS_RMSE(4:6));
    EKF_RMSE_Avg_vel = mean(EKF_RMSE(4:6));
    HCW_RMSE_Avg_vel = mean(HCW_RMSE(4:6));
    
    GPS_x_std = cov(GPS_output(:,1)-r_clean(:,1))^0.5;
    GPS_y_std = cov(GPS_output(:,2)-r_clean(:,2))^0.5;
    GPS_z_std = cov(GPS_output(:,3)-r_clean(:,3))^0.5;  
    
    GPS_xdot_std = cov(GPS_output(:,4)-v_clean(:,1))^0.5;
    GPS_ydot_std = cov(GPS_output(:,5)-v_clean(:,2))^0.5;
    GPS_zdot_std = cov(GPS_output(:,6)-v_clean(:,3))^0.5;  
    
    %fprintf('--------------------------------------------------------------\n')
    %fprintf(' %18.12f     %18.12f     %18.12f   <- Pos&Vel Avg \n', NERM_RMSE_Avg, GPS_RMSE_Avg, EKF_RMSE_Avg)
    %fprintf(' %12.6f x 10^0    %12.6f x 10^-3   <- RMSE EKF Improvement \n',...
    %                    (NERM_RMSE_Avg - EKF_RMSE_Avg)*10^0, (GPS_RMSE_Avg - EKF_RMSE_Avg)*10^3)
    fprintf('--------------------------------------------------------------\n')
    fprintf(' %18.12f     %18.12f     %18.12f   <- Position RMSE Avg \n', HCW_RMSE_Avg_pos, GPS_RMSE_Avg_pos, EKF_RMSE_Avg_pos)
    fprintf(' %18.12f     %18.12f     %18.12f   <- Velocity RMSE Avg \n', HCW_RMSE_Avg_vel, GPS_RMSE_Avg_vel, EKF_RMSE_Avg_vel)
    fprintf('--------------------------------------------------------------\n')
    fprintf(' %12.6f x 10^0    %12.6f x 10^-3   <- Position RMSE EKF Improvement \n',...
                        (HCW_RMSE_Avg_pos - EKF_RMSE_Avg_pos)*10^0, (GPS_RMSE_Avg_pos - EKF_RMSE_Avg_pos)*10^3)
    fprintf(' %12.6f x 10^0    %12.6f x 10^-3   <- Velocity RMSE EKF Improvement \n\n',...
                        (HCW_RMSE_Avg_vel - EKF_RMSE_Avg_vel)*10^0, (GPS_RMSE_Avg_vel - EKF_RMSE_Avg_vel)*10^3)
    %Having the differnce in RMSE_Avg be positive is good in this case,
    %since it shows that the EKF error is lower than the other method
end

%% ========================================================================
%Covariance Analysis
P_start     = 7; %Index for first column of covariance output data
P_x         = EKF_output(:,P_start); 
P_y         = EKF_output(:,P_start+7); 
P_z         = EKF_output(:,P_start+14); 
P_xdot      = EKF_output(:,P_start+21); 
P_ydot      = EKF_output(:,P_start+28); 
P_zdot      = EKF_output(:,P_start+35); 

sig_x         = sqrt(P_x);
sig_y         = sqrt(P_y);
sig_z         = sqrt(P_z);
sig_xdot      = sqrt(P_xdot);
sig_ydot      = sqrt(P_ydot);
sig_zdot      = sqrt(P_zdot);

%% ========================================================================
clearvars -except   time filename make_plots print_plots ...
                    EKF_output GPS_output HCW_output...
                    EKF_RMSE EKF_RMSE_Avg EKF_RMSE_Avg_pos EKF_RMSE_Avg_vel... 
                    GPS_RMSE GPS_RMSE_Avg GPS_RMSE_Avg_pos GPS_RMSE_Avg_vel...
                    HCW_RMSE HCW_RMSE_Avg HCW_RMSE_Avg_pos HCW_RMSE_Avg_vel...
                    errors_x errors_y errors_z ...
                    errors_xdot errors_ydot errors_zdot ...
                    P_x P_y P_z P_xdot P_ydot P_zdot ...
                    sig_x sig_y sig_z sig_xdot sig_ydot sig_zdot...
                    GPS_x_std GPS_y_std GPS_z_std ...
                    GPS_xdot_std GPS_ydot_std GPS_zdot_std...
                    r_clean v_clean...
                    state_priori;
                    
save(filename)
%% ========================================================================
%Plotting the Results
% =========================================================================
if strcmp(make_plots, 'on')           
    FFNAV_HCW_EKF_Plotter(filename)
end
% =========================================================================
