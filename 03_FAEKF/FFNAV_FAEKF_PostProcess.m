function [] = FFNAV_FAEKF_PostProcess(file_data)
% FFNAV FAEKF Data Processing and Analysis ================================
% Description: This script takes a .mat results file created by the FAEKF
% simulation, computes state errors and RMSE.
%
% Inputs:
%   file_data - .mat containing data from the EKF simulation
%   show_RMS  - Turns RMS calculations on/off
%
% Outputs:
%   RMS Error Analysis 
%   Result plots (optional)
%
% Other Functions Called:
%
% Created by:  Cory Fraser - OCT 31, 2017
% Latest Edit: Cory Fraser - OCT 17, 2018
% Copyright(c) 2018 by Cory Fraser
% =========================================================================

%% Initialization
fprintf('\n--------------------------------------------------------------\n')
fprintf(' \t\t FFNAV SIMULATION POST-PROCESSING                     \n')
fprintf('--------------------------------------------------------------\n')
load(file_data)

labels = {'x', 'y', 'z', 'xdot', 'ydot', 'zdot'};
% =========================================================================
%% RMS Analysis of Position and Velocity Errors
fprintf(' Processing Window: %i-%i of %i \n', n_start, n_end, size(EKF_output,1)) 

errors_x         = r_clean(:,1)       - EKF_output(:,1); 
errors_y         = r_clean(:,2)       - EKF_output(:,2); 
errors_z         = r_clean(:,3)       - EKF_output(:,3); 
errors_theta     = theta_clean(:)     - EKF_output(:,4); 
errors_rt        = rt_clean(:)        - EKF_output(:,5); 
errors_xdot      = v_clean(:,1)       - EKF_output(:,6); 
errors_ydot      = v_clean(:,2)       - EKF_output(:,7); 
errors_zdot      = v_clean(:,3)       - EKF_output(:,8); 
errors_theta_dot = theta_dot_clean(:) - EKF_output(:,9); 
errors_rt_dot    = rt_dot_clean(:)    - EKF_output(:,10); 

%Min and max distance between spacecraft
dis_min = min(sqrt(r_clean(:,1).^2+r_clean(:,2).^2+r_clean(:,3).^2));
dis_max = max(sqrt(r_clean(:,1).^2+r_clean(:,2).^2+r_clean(:,3).^2));


%Min and max speed between spacecraft
spd_min = min(sqrt(v_clean(:,1).^2+v_clean(:,2).^2+v_clean(:,3).^2));
spd_max = max(sqrt(v_clean(:,1).^2+v_clean(:,2).^2+v_clean(:,3).^2));

%Examination of the RMS error of the Nosiy GPS and EKF signals
format bank
fprintf(' Root Mean Squared Errors [m, m/s] \n\n')
fprintf('                   Noisy GPS               EKF Output \n')

%Position and velocity errors, for points n_start-to-n_end
%GPS_error = zeros(n_end-n_sta
for k = 1:1:6
    if k < 4
        GPS_error(n_start:n_end,k)   = r_clean(n_start:n_end,k) - GPS_output(n_start:n_end,k);
        EKF_error(n_start:n_end,k)   = r_clean(n_start:n_end,k) - EKF_output(n_start:n_end,k);
%        NERM_error(n_start:n_end,k)  = r_clean(n_start:n_end,k) - NERM_output(n_start:n_end,k);
    else
        GPS_error(n_start:n_end,k)   = v_clean(n_start:n_end,k-3) - GPS_output(n_start:n_end,k+1);
        EKF_error(n_start:n_end,k)   = v_clean(n_start:n_end,k-3) - EKF_output(n_start:n_end,k+2);
%        NERM_error(n_start:n_end,k)  = v_clean(n_start:n_end,k-3) - NERM_output(n_start:n_end,k+2);
    end
    EKF_RMSE(k) = sqrt(mean(EKF_error(n_start:n_end,k).^2));
    GPS_RMSE(k) = sqrt(mean(GPS_error(n_start:n_end,k).^2));
%    NERM_RMSE(k) = sqrt(mean(NERM_error(n_start:n_end,k).^2));
    
    fprintf('            %4s  %8.4f       %18.4f \n', char(labels(k)), GPS_RMSE(k), EKF_RMSE(k))
end

% Averaged Root Mean Square (RMS) Errors
GPS_RMSE_Avg = mean(GPS_RMSE);
EKF_RMSE_Avg = mean(EKF_RMSE);
%NERM_RMSE_Avg = mean(NERM_RMSE);

% Position Average Root Mean Square (RMS) Errors
GPS_RMSE_Avg_pos = mean(GPS_RMSE(1:3));
EKF_RMSE_Avg_pos = mean(EKF_RMSE(1:3));
%NERM_RMSE_Avg_pos = mean(NERM_RMSE(1:3));

% Velocity Average Root Mean Square (RMS) Errors
GPS_RMSE_Avg_vel = mean(GPS_RMSE(4:6));
EKF_RMSE_Avg_vel = mean(EKF_RMSE(4:6));
%NERM_RMSE_Avg_vel = mean(NERM_RMSE(4:6));

%============================
%3D RMS Error Calculations

% Position 3D Root Mean Square (3D-RMS) Errors
GPS_3D_RMS_pos = sqrt(sum(GPS_error(n_start:n_end,1).^2+GPS_error(n_start:n_end,2).^2+GPS_error(n_start:n_end,3).^2)/length(GPS_error(n_start:n_end,k).^2));
EKF_3D_RMS_pos = sqrt(sum(EKF_error(n_start:n_end,1).^2+EKF_error(n_start:n_end,2).^2+EKF_error(n_start:n_end,3).^2)/length(EKF_error(n_start:n_end,k).^2));
%NERM_3D_RMS_pos = sqrt(sum(NERM_error(n_start:n_end,1).^2+NERM_error(n_start:n_end,2).^2+NERM_error(n_start:n_end,3).^2)/length(NERM_error(n_start:n_end,k).^2));

% Position 3D Root Mean Square (3D-RMS) Percent Errors
EKF_3D_RMS_pos_ErrorPercent = EKF_3D_RMS_pos/dis_min*100;
GPS_3D_RMS_pos_ErrorPercent = GPS_3D_RMS_pos/dis_min*100;

% Velocity 3D Root Mean Square (3D-RMS) Errors
GPS_3D_RMS_vel = sqrt(sum(GPS_error(n_start:n_end,4).^2+GPS_error(n_start:n_end,5).^2+GPS_error(n_start:n_end,6).^2)/length(GPS_error(n_start:n_end,k).^2));
EKF_3D_RMS_vel = sqrt(sum(EKF_error(n_start:n_end,4).^2+EKF_error(n_start:n_end,5).^2+EKF_error(n_start:n_end,6).^2)/length(EKF_error(n_start:n_end,k).^2));
%NERM_3D_RMS_vel = sqrt(sum(NERM_error(n_start:n_end,4).^2+NERM_error(n_start:n_end,5).^2+NERM_error(n_start:n_end,6).^2)/length(NERM_error(n_start:n_end,k).^2));

% Velocity 3D Root Mean Square (3D-RMS) Percent Errors
EKF_3D_RMS_vel_ErrorPercent = EKF_3D_RMS_vel/spd_min*100;
GPS_3D_RMS_vel_ErrorPercent = GPS_3D_RMS_vel/spd_min*100;



fprintf('--------------------------------------------------------------\n')
fprintf(' Position 3D RMS %9.4f (%5.2f %%)         %6.4f (%5.2f %%)   \n', GPS_3D_RMS_pos, GPS_3D_RMS_pos_ErrorPercent, EKF_3D_RMS_pos, EKF_3D_RMS_pos_ErrorPercent)
fprintf(' Velocity 3D RMS %9.4f (%5.2f %%)         %6.4f (%5.2f %%)    \n\n', GPS_3D_RMS_vel, GPS_3D_RMS_vel_ErrorPercent,  EKF_3D_RMS_vel, EKF_3D_RMS_vel_ErrorPercent)
fprintf('Position Avg RMS %9.4f       %18.4f    \n', GPS_RMSE_Avg_pos, EKF_RMSE_Avg_pos)
fprintf('Velocity Avg RMS %9.4f       %18.4f    \n', GPS_RMSE_Avg_vel, EKF_RMSE_Avg_vel)
fprintf('--------------------------------------------------------------\n')
fprintf('       EKF Position Improvement (Avg RMS) -> %6.4f  (+ Good)   \n',...
    (GPS_RMSE_Avg_pos - EKF_RMSE_Avg_pos))
fprintf('       EKF Velocity Improvement (Avg RMS) -> %6.4f  (+ Good)  \n',...
    (GPS_RMSE_Avg_vel - EKF_RMSE_Avg_vel))
fprintf('--------------------------------------------------------------\n')

% =========================================================================
%% Covariance Analysis

if strcmp(cov_flag, 'on')
    
    P_start     = 11; %Index for first column of covariance output data
    P_x         = EKF_output(:,P_start); 
    P_y         = EKF_output(:,P_start+11); 
    P_z         = EKF_output(:,P_start+22); 
    P_theta     = EKF_output(:,P_start+33); 
    P_rt        = EKF_output(:,P_start+44); 
    P_xdot      = EKF_output(:,P_start+55); 
    P_ydot      = EKF_output(:,P_start+66); 
    P_zdot      = EKF_output(:,P_start+77); 
    P_theta_dot = EKF_output(:,P_start+88);  
    P_rt_dot    = EKF_output(:,P_start+99); 

    sig_x         = sqrt(P_x);
    sig_y         = sqrt(P_y);
    sig_z         = sqrt(P_z);
    sig_theta     = sqrt(P_theta);
    sig_rt        = sqrt(P_rt);
    sig_xdot      = sqrt(P_xdot);
    sig_ydot      = sqrt(P_ydot);
    sig_zdot      = sqrt(P_zdot);
    sig_theta_dot = sqrt(P_theta_dot);
    sig_rt_dot    = sqrt(P_rt_dot);
    
    % Calculating standard deviation of GPS Errors
    %{
    GPS_cov = diag(cov(GPS_error));
    GPS_std = sqrt(GPS_cov);

    %Creating the running standard deviation
    GPS_cov_sample = zeros(length(GPS_error), 6);
    GPS_std_sample = GPS_cov_sample;
    for k = 2:length(GPS_error)
        GPS_cov_sample(k,:) = (diag(cov(GPS_error(1:k,:))))';
        GPS_std_sample(k,:) = sqrt(GPS_cov_sample(k,:));
    end
    %}

end
%save(file_data)
fprintf('--------------------------------------------------------------\n')
end
% =========================================================================
