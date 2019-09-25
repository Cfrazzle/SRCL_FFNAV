function [] = FFNAV_EKF_PostProcess(file_data)
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
%   EKF_Plotter
%
% Created by: Cory Fraser - July 2016
% Latest Edit: Cory Fraser - June 26, 2017
%% ========================================================================
%Initialization
fprintf('\n--------------------------------------------------------------\n')
fprintf('   FFNAV EKF SIMULATION POST-PROCESSING                       \n')
fprintf('--------------------------------------------------------------\n')
load(file_data)
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
errors_theta     = EKF_output(:,4)  - theta_clean; 
errors_rt        = EKF_output(:,5)  - rt_clean; 
errors_xdot      = EKF_output(:,6)  - v_clean(:,1); 
errors_ydot      = EKF_output(:,7)  - v_clean(:,2); 
errors_zdot      = EKF_output(:,8)  - v_clean(:,3); 
errors_theta_dot = EKF_output(:,9)  - theta_dot_clean; 
errors_rt_dot    = EKF_output(:,10) - rt_dot_clean; 

%Examination of the RMS error of the Nosiy GPS and EKF signals
if strcmp(show_RMS, 'on') 
    format long
    fprintf(' Root-Mean-Squared Errors for Noisy GPS and EKF Output [km, km/s] \n')
    fprintf(' Row: 1 = x, 2 = y, 3 = z, 4 = xdot, 5 = ydot, 6 = zdot  \n\n')
    fprintf(' Nonlinear Dynamics           Noisy GPS               EKF Output \n')    
    
    for k = 1:1:6
        if k < 4
            GPS_error   = GPS_output(:,k) - r_clean(:,k);
            EKF_error   = EKF_output(:,k) - r_clean(:,k);
            NERM_error  = NERM_output(:,k) - r_clean(:,k);
        else 
            GPS_error   = GPS_output(:,k+1) - v_clean(:,k-3);
            EKF_error   = EKF_output(:,k+2) - v_clean(:,k-3);
            NERM_error  = NERM_output(:,k+2) - v_clean(:,k-3);
        end
        EKF_RMSE(k) = sqrt(mean(EKF_error.^2));
        GPS_RMSE(k) = sqrt(mean(GPS_error.^2));
        NERM_RMSE(k) = sqrt(mean(NERM_error.^2));

        fprintf(' %18.12f     %18.12f     %18.12f \n', NERM_RMSE(k), GPS_RMSE(k), EKF_RMSE(k))
    end   
        
    GPS_RMSE_Avg = mean(GPS_RMSE);
    EKF_RMSE_Avg = mean(EKF_RMSE);
    NERM_RMSE_Avg = mean(NERM_RMSE);
    
    GPS_RMSE_Avg_pos = mean(GPS_RMSE(1:3));
    EKF_RMSE_Avg_pos = mean(EKF_RMSE(1:3));
    NERM_RMSE_Avg_pos = mean(NERM_RMSE(1:3));
    
    GPS_RMSE_Avg_vel = mean(GPS_RMSE(4:6));
    EKF_RMSE_Avg_vel = mean(EKF_RMSE(4:6));
    NERM_RMSE_Avg_vel = mean(NERM_RMSE(4:6));
    
    
    GPS_x_std = cov(GPS_output(:,1)-r_clean(:,1))^0.5;
    GPS_y_std = cov(GPS_output(:,2)-r_clean(:,2))^0.5;
    GPS_z_std = cov(GPS_output(:,3)-r_clean(:,3))^0.5;  
    
    GPS_xdot_std = cov(GPS_output(:,5)-v_clean(:,1))^0.5;
    GPS_ydot_std = cov(GPS_output(:,6)-v_clean(:,2))^0.5;
    GPS_zdot_std = cov(GPS_output(:,7)-v_clean(:,3))^0.5;  
    
    %fprintf('--------------------------------------------------------------\n')
    %fprintf(' %18.12f     %18.12f     %18.12f   <- Pos&Vel Avg \n', NERM_RMSE_Avg, GPS_RMSE_Avg, EKF_RMSE_Avg)
    %fprintf(' %12.6f x 10^0    %12.6f x 10^-3   <- RMSE EKF Improvement \n',...
    %                    (NERM_RMSE_Avg - EKF_RMSE_Avg)*10^0, (GPS_RMSE_Avg - EKF_RMSE_Avg)*10^3)
    fprintf('--------------------------------------------------------------\n')
    fprintf(' %18.12f     %18.12f     %18.12f   <- Position RMSE Avg \n', NERM_RMSE_Avg_pos, GPS_RMSE_Avg_pos, EKF_RMSE_Avg_pos)
    fprintf(' %18.12f     %18.12f     %18.12f   <- Velocity RMSE Avg \n', NERM_RMSE_Avg_vel, GPS_RMSE_Avg_vel, EKF_RMSE_Avg_vel)
    fprintf('--------------------------------------------------------------\n')
    fprintf(' %12.6f x 10^0    %12.6f x 10^-3   <- Position RMSE EKF Improvement \n',...
                        (NERM_RMSE_Avg_pos - EKF_RMSE_Avg_pos)*10^0, (GPS_RMSE_Avg_pos - EKF_RMSE_Avg_pos)*10^3)
    fprintf(' %12.6f x 10^0    %12.6f x 10^-3   <- Velocity RMSE EKF Improvement \n\n',...
                        (NERM_RMSE_Avg_vel - EKF_RMSE_Avg_vel)*10^0, (GPS_RMSE_Avg_vel - EKF_RMSE_Avg_vel)*10^3)
    %Having the differnce in RMSE_Avg be positive is good in this case,
    %since it shows that the EKF error is lower than the other method
end

%% ========================================================================
%Covariance Analysis
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

L_sum(1) = L(1);
for i = 1:length(time)
    trace_P(i) = sum([P_x(i) P_y(i) P_z(i) P_theta(i) P_rt(i) P_xdot(i) ...
        P_ydot(i) P_zdot(i) P_theta_dot(i) P_rt_dot(i)]);
    if i > 1
        L_sum(i) = L(i) + L(i-1);
    end
end

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


%% ========================================================================
clearvars -except   EKF_output NERM_output GPS_output time ...
                    file_data pathname...
                    r_clean v_clean...
                    make_plots print_plots ...
                    theta_target rt_clean ...
                    state_priori...
                    EKF_RMSE EKF_RMSE_Avg EKF_RMSE_Avg_pos EKF_RMSE_Avg_vel... 
                    GPS_RMSE GPS_RMSE_Avg GPS_RMSE_Avg_pos GPS_RMSE_Avg_vel...
                    NERM_RMSE NERM_RMSE_Avg NERM_RMSE_Avg_pos NERM_RMSE_Avg_vel...
                    errors_x errors_y errors_z errors_theta errors_rt...
                    errors_xdot errors_ydot errors_zdot errors_theta_dot errors_rt_dot...
                    P_x P_y P_z P_theta P_rt P_xdot P_ydot P_zdot P_theta_dot P_rt_dot...
                    sig_x sig_y sig_z sig_theta sig_rt sig_xdot sig_ydot sig_zdot...
                    sig_theta_dot sig_rt_dot...
                    GPS_x_std GPS_y_std GPS_z_std ...
                    GPS_xdot_std GPS_ydot_std GPS_zdot_std ...
                    theta_clean theta_dot_clean rt_clean rt_dot_clean...
                    Pv_obs Pv_smoothed Pr_theo P_error Pdot_error...
                    K_pe K_Ie K_e...
                    Q0 R0 Q_output R_output trace_P...
                    L L_sum; 
save(file_data)
%% ========================================================================
%Plotting the Results
% =========================================================================
if strcmp(make_plots, 'on')           
    FFNAV_EKF_Plotter(file_data)
end
% =========================================================================
