function [] = FFNAV_FAEKF_Plotter(output_file)
%% FFNAV_EKF_plots=========================================================
% Description: This function contstructs plots for the FFNAV FAEKF
% simulation, based on a passed .mat file that can be created by the
% FAEKF_Run.m program.
%
% Inputs:
%   output_file - .mat file containing simulation results of the EKF
%   time_step   - Time step of the simulation
%   time_end    - Duration of the simulation
%   make_plots  - Turns plot output on/off
%   print_flag  - Turns .eps plot printing on/off
%   show_RMS    - Turns RMS claculations on/off
%   EKF_flag    - EKF algorithm used (0 = EKF, 1 = MLE_AEKF, 2 = FAEKF) 
%
% Outputs:
%   Result plots (optional)
%   
% Other Functions Called:
%   None
%
% Created by:  Cory Fraser - OCT 31, 2017
% Latest Edit: Cory Fraser - OCT 24, 2018
% Copyright(c) 2018 by Cory Fraser
% =========================================================================

%% Initialization, set default characteristics
load(output_file)
mydir  = pwd;                       %Get current directory
idcs   = strfind(mydir,filesep);    %Index file separations
newdir = mydir(1:idcs(end)-1);      %Path to higher folder (1-level up)
addpath(newdir)                     %Add path to location of MakeParams
    
TF_size       = 24;  %Title font size
AF_size       = 20;  %Axes label font size
AT_size       = 20;  %Tick label font size
Line_width    = 2;   %Line width

set(0,'DefaultTextFontName', 'Times New Roman')
set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', AF_size)
set(0,'DefaultTextFontSize', TF_size)
set(0,'defaultLineLineWidth', Line_width)


fig_id = 1; %Counter for figures

plot_start  = 1;                %Starting time for plots
plot_end    = length(time);     %Ending time for plots
time_plot   = time/3600;        %Scale plot time from seconds to hours

%% Select Plots ===========================================================

% 3D Position Figures
EKF_3D_plot     = 0;        % 3D Positions from EKF estimation
GPS_3D_plot     = 0;        % 3D Positions from GPS Measurements
NERM_3D_plot    = 0;        % 3D Positions from Dynamics Propagation
Truth_3D_plot   = 1;        % 3D Positions from Truth Orbit Propagation

% EKF State Estimates Time History
EKF_subplot      = 0;        % Subplots of the EKF state  
EKF_x_plot       = 0;        % EKF x-position state 
EKF_y_plot       = 0;        % EKF y-position state 
EKF_z_plot       = 0;        % EKF z-position state 
EKF_theta_plot   = 0;        % EKF theta state 
EKF_rt_plot      = 0;        % EKF rt state 
EKF_xd_plot      = 0;        % EKF x-velocity state 
EKF_yd_plot      = 0;        % EKF y-velocity state 
EKF_zd_plot      = 0;        % EKF z-velocity state 
EKF_thetad_plot  = 0;        % EKF theta-dot state 
EKF_rtd_plot     = 0;        % EKF rt-dot state 

% EKF State Error Time History
EKF_error_subplot    = 0;        % Subplots of the EKF state errors 
EKF_x_err_plot       = 0;        % EKF x-position state errors
EKF_y_err_plot       = 0;        % EKF y-position state errors
EKF_z_err_plot       = 0;        % EKF z-position state errors
EKF_theta_err_plot   = 0;        % EKF theta state errors
EKF_rt_err_plot      = 0;        % EKF rt state errors
EKF_xd_err_plot      = 0;        % EKF x-velocity state errors
EKF_yd_err_plot      = 0;        % EKF y-velocity state errors
EKF_zd_err_plot      = 0;        % EKF z-velocity state errors
EKF_thetad_err_plot  = 0;        % EKF theta-dot state errors
EKF_rtd_err_plot     = 0;        % EKF rt-dot state errors

%Covariance Plot Time History
Q_cov_subplot       = 0;         % Subplots of the Q Covariances
R_cov_subplot       = 0;         % Subplots of the R Covariances


%Other Time History Plots
RMSbar_plot         = 0;         % Bar graphs of RMS errors
DOM_subplot         = 0;         % FAEKF Innovations Mismatch
deltaDOD_plot       = 0;         % FAEKF Degree of Divergence Difference
NIS_plot            = 0;         % FAEKF Normalized Innovations Squared

%Animations
RelativeMotion_Movie       = 1;  % 3D Positions w/ Earth
RelativePosition_Animation = 0;  % 3D LVLH Position
RelativeVelocity_Animation = 0;  % 3D LVLH Velocity

%Add covariances to plots
Cov_plots = strcmp(cov_flag, 'on');

% EKF State Error Histograms
EKF_x_err_histo       = 0;        % EKF x-position state errors
EKF_y_err_histo       = 0;        % EKF y-position state errors
EKF_z_err_histo       = 0;        % EKF z-position state errors
EKF_theta_err_histo   = 0;        % EKF theta state errors
EKF_rt_err_histo      = 0;        % EKF rt state errors
EKF_xd_err_histo      = 0;        % EKF x-velocity state errors
EKF_yd_err_histo      = 0;        % EKF y-velocity state errors
EKF_zd_err_histo      = 0;        % EKF z-velocity state errors
EKF_thetad_err_histo  = 0;        % EKF theta-dot state errors
EKF_rtd_err_histo     = 0;        % EKF rt-dot state errors

% EKF States and Errors, Side-by-side
EKF_StatesError_pos = 0;
EKF_StatesError_vel = 0;

%% Plotting Routines ======================================================

%Figure: 3D Position Plot for EKF
if (EKF_3D_plot)
    figure(fig_id); fig_id = fig_id + 1;
        plot3(EKF_output(:,1), EKF_output(:,2), EKF_output(:,3));
        title('EKF Relative Position Estimates')
        xlab = xlabel('X Position (m)');
        ylabel('Y Position (m)')
        zlabel('Z Position (m)')
        set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
        axis([-1.2*max(abs(GPS_output(:,2))) 1.2*max(abs(GPS_output(:,2))) -1.2*max(abs(GPS_output(:,2)))...
            1.2*max(abs(GPS_output(:,2))) -1.2*max(abs(GPS_output(:,3))) 1.2*max(abs(GPS_output(:,3)))]);
        view(-35,50)
        grid on
        hold on
        Target   = scatter3(0, 0, 0, 'filled', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'MarkerEdgeColor', 'k', 'SizeData', 5^2);
        plot3(r_clean(:,1), r_clean(:,2), r_clean(:,3));
        legend('Chaser Trajectory', 'Target Spacecraft', 'Truth Trajectory','Location','NorthEast')
        if strcmp(print_flag, 'on')
            print -deps EKF_Output
        end
end

%Figure: 3D Position Plot for GPS
if (GPS_3D_plot)
    figure(fig_id); fig_id = fig_id + 1;
        plot3(GPS_output(:,1),GPS_output(:,2),GPS_output(:,3));
        title('GPS Relative Position Measurements')
        xlab = xlabel('x Position (m)');
        ylabel('y Position (m)')
        zlabel('z Position (m)')
        set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
        axis([-1.2*max(abs(GPS_output(:,2))) 1.2*max(abs(GPS_output(:,2))) -1.2*max(abs(GPS_output(:,2)))...
            1.2*max(abs(GPS_output(:,2))) -1.2*max(abs(GPS_output(:,3))) 1.2*max(abs(GPS_output(:,3)))]);
        view(-35,50)
        grid on
        hold on
        Target   = scatter3(0, 0, 0, 'filled', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'MarkerEdgeColor', 'k', 'SizeData', 5^2);
        legend('Chaser Trajectory', 'Target Spacecraft', 'Location','NorthEast')
        %file_print = fullfile(pathname, 'Fig_NoisyGPS');
        %print(file_print,  '-deps');   
end

%Figure: 3D Position Plot NERMs
if (NERM_3D_plot)
    figure(fig_id); fig_id = fig_id + 1;
        plot3(NERM_output(:,1), NERM_output(:,2), NERM_output(:,3));
        title('NERM Dynamics Relative Position Estimates')
        xlab = xlabel('X Position (m)');
        ylabel('Y Position (m)')
        zlabel('Z Position (m)')
        set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
        lim_NERM = 1.2*max(max(abs(NERM_output)));
       % axis([-lim_NERM lim_NERM -lim_NERM lim_NERM -lim_NERM lim_NERM]);
        axis([-1.2*max(abs(GPS_output(:,2))) 1.2*max(abs(GPS_output(:,2))) -1.2*max(abs(GPS_output(:,2)))...
            1.2*max(abs(GPS_output(:,2))) -1.2*max(abs(GPS_output(:,3))) 1.2*max(abs(GPS_output(:,3)))]);
        view(-35,50)
        grid on
        hold on
        Target   = scatter3(0, 0, 0, 'filled', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'MarkerEdgeColor', 'k', 'SizeData', 5^2);
        legend('Chaser Trajectory', 'Target Spacecraft', 'Location','NorthEast')
        if strcmp(print_flag, 'on')
            print -deps NERM_Prop
        end
end

%Figure: Truth 3D Position Plot
if (Truth_3D_plot)
    figure(fig_id); fig_id = fig_id + 1;
        plot3(r_clean(:,1), r_clean(:,2), r_clean(:,3));
        title('Orbit Propagator - Relative Position Truth')
        xlab = xlabel('$x$ Position (m)');
        ylab = ylabel('$y$ Position (m)');
        zlab = zlabel('$z$ Position (m)');
        set(gcf, 'units','normalized','outerposition',[0 0 0.75 1]);
        axis([-1.2*max(abs(GPS_output(:,2))) 1.2*max(abs(GPS_output(:,2))) -1.2*max(abs(GPS_output(:,2)))...
            1.2*max(abs(GPS_output(:,2))) -1.2*max(abs(GPS_output(:,3))) 1.2*max(abs(GPS_output(:,3)))]);
        view(-35,50)
        grid on
        hold on
        Target   = scatter3(0, 0, 0, 'filled', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'MarkerEdgeColor', 'k', 'SizeData', 5^2);
        set(xlab,'Interpreter','latex');
        set(ylab,'Interpreter','latex');
        set(zlab,'Interpreter','latex');
        legend('Chaser Trajectory', 'Target Spacecraft', 'Location','NorthEast')
        %NiceLegend({'Chaser Trajectory', 'Target Spacecraft'});
        if strcmp(print_flag, 'on')
            print -deps Truth_Relative Trajectory
        end
        
    %3-View Plots
    figure(fig_id); fig_id = fig_id + 1;
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);    
    subplot(2,2,1)
        plot(r_clean(:,1), r_clean(:,2));
        %axis([-1.2*max(abs(GPS_output(:,2))) 1.2*max(abs(GPS_output(:,2))) -1.2*max(abs(GPS_output(:,2)))...
        %    1.2*max(abs(GPS_output(:,2))) -1.2*max(abs(GPS_output(:,3))) 1.2*max(abs(GPS_output(:,3)))]);
        grid on
        hold on
        axis equal
        Target   = scatter(0, 0, 'filled', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'MarkerEdgeColor', 'k', 'SizeData', 5^2);
        xlab = xlabel('$x$ Position (m)');
        ylab = ylabel('$y$ Position (m)');
        set(xlab,'Interpreter','latex');
        set(ylab,'Interpreter','latex');
        
    subplot(2,2,2)
        plot3(r_clean(:,1), r_clean(:,2), r_clean(:,3));
        xlab = xlabel('$x$ (m)');
        ylab = ylabel('$y$ (m)');
        zlab = zlabel('$z$ (m)');
        set(xlab,'Interpreter','latex');
        set(ylab,'Interpreter','latex');
        set(zlab,'Interpreter','latex');
        %set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
        %axis([-1.2*max(abs(GPS_output(:,2))) 1.2*max(abs(GPS_output(:,2))) -1.2*max(abs(GPS_output(:,2)))...
        %    1.2*max(abs(GPS_output(:,2))) -1.2*max(abs(GPS_output(:,3))) 1.2*max(abs(GPS_output(:,3)))]);
        view(-45,45)
        grid on
        hold on
        Target   = scatter3(0, 0, 0, 'filled', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'MarkerEdgeColor', 'k', 'SizeData', 5^2);
        %legend('Chaser Trajectory', 'Target Spacecraft', 'Location','NorthEast')
        NiceLegend({'Chaser Trajectory', 'Target Spacecraft'});
        
    subplot(2,2,3)
        plot(r_clean(:,1), r_clean(:,3));
        xlab = xlabel('$x$ Position (m)');
        ylab = ylabel('$z$ Position (m)');
        set(xlab,'Interpreter','latex');
        set(ylab,'Interpreter','latex');
        axis equal
        % ylim([-50 50])
        %axis([-1.2*max(abs(GPS_output(:,2))) 1.2*max(abs(GPS_output(:,2))) -1.2*max(abs(GPS_output(:,2)))...
        %    1.2*max(abs(GPS_output(:,2))) -1.2*max(abs(GPS_output(:,3))) 1.2*max(abs(GPS_output(:,3)))]);
        grid on
        hold on
        Target   = scatter3(0, 0, 0, 'filled', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'MarkerEdgeColor', 'k', 'SizeData', 5^2);        
        
    subplot(2,2,4)
        plot(r_clean(:,2), r_clean(:,3));
        xlab = xlabel('$y$ Position (m)');
        ylab = ylabel('$z$ Position (m)');
        set(xlab,'Interpreter','latex');
        set(ylab,'Interpreter','latex');
        axis equal
        %ylim([-50 50])
        %axis([-1.2*max(abs(GPS_output(:,2))) 1.2*max(abs(GPS_output(:,2))) -1.2*max(abs(GPS_output(:,2)))...
        %    1.2*max(abs(GPS_output(:,2))) -1.2*max(abs(GPS_output(:,3))) 1.2*max(abs(GPS_output(:,3)))]);
        grid on
        hold on
        Target   = scatter(0, 0, 'filled', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'MarkerEdgeColor', 'k', 'SizeData', 5^2);
        if strcmp(print_flag, 'on')
            print -deps Truth_Relative_Orbit3View
        end
end

%% ========================================================================  
%Figure w/ Subfigures: EKF State Outputs
if (EKF_subplot)

    figure(fig_id); fig_id = fig_id + 1;
    %title('EKF Output States')    
        h(1) = subplot(5,2,1);
        plot(time_plot, EKF_output(:,1)) 
        set(gca,'FontSize',AT_size);
        title('EKF Output States')
        ylabel('X(m)')
        set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
        set(gca,'XTick',[]);
        %set(gca,'YTick',[]);
        grid on
        hold on
        xlim([0 inf])

        h(2) = subplot(5,2,3);
        plot(time_plot, EKF_output(:,2)) 
        set(gca,'FontSize',AT_size);
        ylabel('Y(m)')
        set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
        set(gca,'XTick',[]);
        %set(gca,'YTick',[]);
        grid on
        hold on
        %axis([0 time(length(time)) -200 200])
        xlim([0 inf])

        h(3) = subplot(5,2,5);
        plot(time_plot, EKF_output(:,3)) 
        set(gca,'FontSize',AT_size);
        ylabel('Z (m)')
        set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
        set(gca,'XTick',[]);
        %set(gca,'YTick',[]);
        grid on
        hold on
        xlim([0 inf])

        h(4) = subplot(5,2,7);
        plot(time_plot, EKF_output(:,4)) 
        set(gca,'FontSize',AT_size);
        ylabel('\theta (rad)')
        set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
        set(gca,'XTick',[]);
        %set(gca,'YTick',[]);
        grid on
        hold on
        xlim([0 inf])

        h(5) = subplot(5,2,9);
        plot(time_plot, EKF_output(:,5)) 
        set(gca,'FontSize',AT_size);
        xlab = xlabel('Time (hrs)');
        ylabel('r_t (m)')
        set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
        grid on
        hold on
        xlim([0 inf])

        h(6) = subplot(5,2,2);
        plot(time_plot, EKF_output(:,6)) 
        set(gca,'FontSize',AT_size);
        ylabel('X_d_o_t (m/s)')
        set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
        set(gca,'XTick',[]);
        %set(gca,'YTick',[]);
        grid on
        hold on
        xlim([0 inf])

        h(7) = subplot(5,2,4);
        plot(time_plot, EKF_output(:,7)) 
        set(gca,'FontSize',AT_size);
        ylabel('Y_d_o_t (m/s)')
        set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
        set(gca,'XTick',[]);
        %set(gca,'YTick',[]);
        grid on
        hold on
        xlim([0 inf])

        h(8) = subplot(5,2,6); 
        plot(time_plot, EKF_output(:,8)) 
        set(gca,'FontSize',AT_size);
        ylabel('Z_d_o_t (m/s)')
        set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
        set(gca,'XTick',[]);
        %set(gca,'YTick',[]);
        grid on
        hold on
        xlim([0 inf])

        h(9) = subplot(5,2,8);
        plot(time_plot, EKF_output(:,9)) 
        set(gca,'FontSize',AT_size);
        ylabel('\theta_d_o_t (rad/s)')
        set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
        set(gca,'XTick',[]);
        %set(gca,'YTick',[]);
        grid on
        hold on
        xlim([0 inf])

        h(10) = subplot(5,2,10);
        plot(time_plot, EKF_output(:,10)) 
        set(gca,'FontSize',AT_size);
        ylabel('r_t_,_d_o_t')
        set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
        %set(gca,'XTick',[]);
        %set(gca,'YTick',[]);
        grid on
        hold on
        xlim([0 inf])
        xlab = xlabel('Time (hrs)');

        for j = 1:10
            p(j,:) = get(h(j), 'pos');
            p(j,4) = p(j,4) + 0.015;
            set(h(j), 'pos', p(j,:));
            set(get(h(j),'YLabel'),'Rotation',0)
        end   
end

%% ========================================================================
% Figures: EKF State Outputs
if EKF_x_plot == 1
    figure(fig_id); fig_id = fig_id + 1;
    plot(time_plot, EKF_output(:,1), '-', time_plot, NERM_output(:,1), '-.', time_plot, GPS_output(:,1), '--')
    legend('EKF Output', 'Dynamics', 'GPS Measurements')
    title('EKF X Position vs Time')
    ylabel('X Position (m)')
    xlab = xlabel('Time (hr)');
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    grid on  
end
    
if EKF_y_plot == 1
    figure(fig_id); fig_id = fig_id + 1;
    plot(time_plot, EKF_output(:,2), '-', time_plot, NERM_output(:,2), '-.', time_plot, GPS_output(:,2), '--')
    legend('EKF Output', 'Dynamics', 'GPS Measurements')
    title('EKF Output Y vs Time')
    xlab = xlabel('Time (hr)');
    ylabel('Y Position (m)')
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    grid on  
end

if EKF_z_plot == 1  
    figure(fig_id); fig_id = fig_id + 1;
    plot(time_plot, EKF_output(:,3), '-', time_plot, NERM_output(:,3), '-.', time_plot, GPS_output(:,3), '--')
    legend('EKF Output', 'Dynamics', 'GPS Measurements')
    title('EKF Output Z vs Time')
    xlab = xlabel('Time (hrs)');
    ylabel('Z Position (m)')
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    grid on  
end

if EKF_theta_plot == 1
    figure(fig_id); fig_id = fig_id + 1;
    plot(time_plot, EKF_output(:,4), '-', time_plot, NERM_output(:,4),'-.', time_plot, theta_target, '--')
    title('EKF Output \theta vs Time')
    xlab = xlabel('Time (hrs)');
    ylabel('\theta (rad)')
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    %axis([0 time(length(time)) -1 20]);
    legend('EKF Output', 'Dynamics', 'Virtual Measurement') 
    grid on  
end
   
if EKF_rt_plot == 1
    figure(fig_id); fig_id = fig_id + 1;
    plot(time_plot, EKF_output(:,5), '-', time_plot, NERM_output(:,5), '-.', time_plot, rt_clean, '--')
    title('EKF Output R_t vs Time')
    xlab = xlabel('Time (hrs)');
    ylabel('R_t (m)')
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    legend('EKF Output', 'Dynamics', 'Virtual Measurement')
    grid on  
end

if EKF_xd_plot == 1
    figure(fig_id); fig_id = fig_id + 1;
    plot(time_plot, EKF_output(:,6)*10^3, '-', time_plot, NERM_output(:,6)*10^3, '-.', time_plot, GPS_output(:,5)*10^3, '--')
    title('EKF Output X-dot vs Time')
    xlab = xlabel('Time (hrs)');
    ylabel('X_d_o_t (mm/s)')
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    legend('EKF Output', 'Dynamics', 'GPS Measurement')
    grid on  
end

if EKF_yd_plot == 1
    figure(fig_id); fig_id = fig_id + 1;
    plot(time_plot, EKF_output(:,7)*10^3, '-', time_plot, NERM_output(:,7), '-.', time_plot, GPS_output(:,6)*10^3, '--')
    title('EKF Output Y-dot vs Time')
    xlab = xlabel('Time (hrs)');
    ylabel('Y_d_o_t (mm/s)')
    legend('EKF Output', 'Dynamics', 'GPS Measurement')
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    grid on  
end

if EKF_zd_plot == 1
    figure(fig_id); fig_id = fig_id + 1;
    plot(time_plot, EKF_output(:,8)*10^3, '-', time_plot, NERM_output(:,8), '-.', time_plot, GPS_output(:,7)*10^3, '--')
    hold on
    title('EKF Output Z-dot vs Time')
    xlab = xlabel('Time (hrs)');
    ylabel('Z_d_o_t (mm/s)')
    legend('EKF Output', 'Dynamics', 'GPS Measurement')
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    grid on  
end

if EKF_thetad_plot == 1
    figure(fig_id); fig_id = fig_id + 1;
    plot(time_plot, EKF_output(:,9), time_plot, NERM_output(:,9), '-.')
    title('EKF Output \theta-dot vs Time')
    xlab = xlabel('Time (hrs)');
    ylabel('\theta_dot (rad/s)')
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    legend('EKF Output', 'Dynamics')
    grid on  
end

if EKF_rtd_plot == 1
    figure(fig_id); fig_id = fig_id + 1;
    plot(time_plot, EKF_output(:,10), time_plot, NERM_output(:,10), '-.')
    title('EKF Output R_t-dot vs Time')
    xlab = xlabel('Time (hrs)');
    ylabel('R_t-dot (m/s)')
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    legend('EKF Output', 'Dynamics')
    grid on  
end  

%% ========================================================================  
% Figure w/ Subfigures: EKF State Errors
if EKF_error_subplot == 1
    figure(fig_id); fig_id = fig_id + 1;
    
    % Figure: EKF X-Position Error
        subplot(5,2,1)
        hold on
        plot(time_plot, errors_x)
        if (Cov_plots)
            plot(time_plot, 3*sig_x, '-.k') 
            plot(time_plot, -3*sig_x,'-.k') 
        end
        set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
        grid on
        xlim([0 inf])
     
    % Figure: EKF Y-Position Error
        subplot(5,2,3)
        hold on
        plot(time_plot, errors_y) 
        if (Cov_plots)
                plot(time_plot, 3*sig_y, '-.k') 
                plot(time_plot, -3*sig_y, '-.k') 
        end
        ylabel('Y Position (m)')
        set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
        grid on
        xlim([0 inf])
        
    % Figure: EKF Z-Position Error
        subplot(5,2,5)
        hold on
        plot(time_plot, errors_z) 
        if (Cov_plots)
            plot(time_plot, 3*sig_z, '-.k')
            plot(time_plot, -3*sig_z, '-.k') 
        end
        ylabel('Z Position (m)')
        set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
        grid on
        xlim([0 inf])
    
    % Figure: EKF Theta Error
    subplot(5,2,7)
        plot(time_plot, errors_theta) 
        hold on
        if (Cov_plots)
            plot(time_plot, 3*sig_theta, '-.k') 
            plot(time_plot, -3*sig_theta, '-.k')
        end
        ylabel('\theta (rad)')
        set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
        grid on
        xlim([0 inf])

    % Figure: EKF rt Error
    subplot(5,2,9)
        plot(time_plot, errors_rt) 
        hold on
        if (Cov_plots)
            plot(time_plot, 3*sig_rt, '-.k')
            plot(time_plot, -3*sig_rt, '-.k')
        end
        xlab = xlabel('Time (hours)');
        ylabel('r_t (m)')
        set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
        grid on
        xlim([0 inf])

    % Figure: EKF X-Velocity Error
    subplot(5,2,2)
        plot(time_plot, errors_xdot) 
        hold on
        if (Cov_plots)
            plot(time_plot, 3*sig_xdot, '-.k') 
            plot(time_plot, -3*sig_xdot, '-.k') 
            legend('EKF Error', '3-\sigma Bounds')
        else
            legend('EKF Error')
        end
        ylabel('X Velocity (m/s)')
        set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
        grid on
        xlim([0 inf])

    % Figure: EKF Y-Velocity Error
    subplot(5,2,4)
        plot(time_plot, errors_ydot) 
        hold on
        if (Cov_plots)
                plot(time_plot, 3*sig_ydot, '-.k') 
                plot(time_plot, -3*sig_ydot, '-.k')
        end
        ylabel('Y Velocity (m/s)')
        set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
        grid on
        xlim([0 inf])

    % Figure: EKF Z-Velocity Error
    subplot(5,2,6)
        plot(time_plot, errors_zdot) 
        hold on
        if (Cov_plots)
            plot(time_plot, 3*sig_zdot, '-.k') 
            plot(time_plot, -3*sig_zdot, '-.k')
        end
        ylabel('Z Velocity (m/s)')
        set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
        grid on
        xlim([0 inf])

    % Figure: EKF Theta_dot Error
    subplot(5,2,8)
        plot(time_plot, errors_theta_dot) 
        hold on
        if (Cov_plots)
            plot(time_plot, 3*sig_theta_dot, '-.k') 
            plot(time_plot, -3*sig_theta_dot, '-.k')
        end
        ylabel('\theta dot(rad/s)')
        set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
        grid on
        xlim([0 inf])

    % Figure: EKF rt_dot Error
    subplot(5,2,10)
        plot(time_plot, errors_rt_dot) 
        hold on
        if (Cov_plots)
            plot(time_plot, 3*sig_rt_dot, '-.k')
            plot(time_plot, -3*sig_rt_dot, '-.k')
        end
        xlab = xlabel('Time (hours)');
        ylabel('r_t dot (m/s)')
        set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
        grid on
        xlim([0 inf])
end       

%% =========================================================================  
% Figures: EKF State Errors

% Figure: EKF X-Position Error
if EKF_x_err_plot == 1
    figure(fig_id); fig_id = fig_id + 1;
        if (Cov_plots)
            fill([time_plot(n_start:n_end)', fliplr(time_plot(n_start:n_end)')], [3*sig_x', fliplr(-3*sig_x')], [0.9 0.9 0.9], 'FaceAlpha', 0.25, 'LineStyle', '-.')
        end
        hold on
        plot(time_plot, errors_x)   
        GPS_err_dots = plot(time_plot(n_start:20:n_end), GPS_error(n_start:20:n_end,1),'.k');
        if (Cov_plots)
            legend('EKF 3-\sigma Bounds', 'EKF Error', 'GPS Errors')
        else
            legend('EKF Error', 'GPS Errors')
        end
        title('EKF-X Error')
        xlab = xlabel('Time (hrs)');
        ylabel('X Position (m)')
        set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
        grid on
        xlim([0 inf])
        if strcmp(print_flag, 'on')
            print -deps EKF_Errors_X
        end
end
%==========================================================================   
% Figure: EKF Y-Position Error
if EKF_y_err_plot == 1

    figure(fig_id); fig_id = fig_id + 1;
        if (Cov_plots)
             fill([time_plot(n_start:n_end)', fliplr(time_plot(n_start:n_end)')], [3*sig_y', fliplr(-3*sig_y')], [0.9 0.9 0.9], 'FaceAlpha', 0.25, 'LineStyle', '-.')
        end
        hold on
        plot(time_plot, errors_y) 
        plot(time_plot(n_start:20:n_end), GPS_error(n_start:20:n_end,2),'.k')
        if (Cov_plots)
            legend('EKF 3-\sigma Bounds', 'EKF Error', 'GPS Errors')
        else
            legend('EKF Error', 'GPS Errors')
        end
                title('EKF-Y Error')
        xlab = xlabel('Time (hrs)');
        ylabel('Y Position (m)')
        set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
        grid on
        xlim([0 inf])
        if strcmp(print_flag, 'on')
            print -deps EKF_Errors_Y
        end
end
%==========================================================================  
% Figure: EKF Z-Position Error
if EKF_z_err_plot == 1

    figure(fig_id); fig_id = fig_id + 1;
        if (Cov_plots)
            fill([time_plot(n_start:n_end)', fliplr(time_plot(n_start:n_end)')], [3*sig_z', fliplr(-3*sig_z')], [0.9 0.9 0.9], 'FaceAlpha', 0.25, 'LineStyle', '-.')
        end
        hold on
        plot(time_plot, errors_z) 
        plot(time_plot(n_start:20:n_end), GPS_error(n_start:20:n_end,3),'.k')
        if (Cov_plots)
            legend('EKF 3-\sigma Bounds', 'EKF Error', 'GPS Errors')
        else
            legend('EKF Error', 'GPS Errors')
        end
        title('EKF-Z Error')
        xlab = xlabel('Time (hrs)');
        ylabel('Z Position (m)')
        set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
        grid on
        xlim([0 inf])
        if strcmp(print_flag, 'on')
            print -deps EKF_Errors_Z
        end
end

%==========================================================================  
% Figure: EKF Theta Error
if EKF_theta_err_plot == 1

    figure(fig_id); fig_id = fig_id + 1;
        plot(time_plot, errors_theta) 
        hold on
        plot(time_plot, 3*sig_theta) 
        plot(time_plot, -3*sig_theta) 
        legend('EKF Error', '3-\sigma Bounds')
        title('EKF \theta Error')
        ylabel('\theta (rad)')
        set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
        grid on
        xlim([0 inf])
end
%==========================================================================  
% Figure: EKF rt Error
if EKF_rt_err_plot == 1

    figure(fig_id); fig_id = fig_id + 1;
        plot(time_plot, errors_rt) 
        hold on
        plot(time_plot, 3*sig_rt)
        plot(time_plot, -3*sig_rt)
        legend('EKF Error', '3-\sigma Bounds')
        title('EKF r_t Error')
        ylabel('r_t (m)')
        set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
        grid on
        xlim([0 inf])
end
%==========================================================================  
% Figure: EKF X-Velocity Error
if EKF_xd_err_plot == 1

    figure(fig_id); fig_id = fig_id + 1;
        plot(time_plot, errors_xdot) 
        hold on
        plot(time_plot, 3*sig_xdot) 
        plot(time_plot, -3*sig_xdot) 
        legend('EKF Error', '3-\sigma Bounds')
        title('EKF X-Velocity Error')
        ylabel('X Velocity (m/s)')
        set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
        grid on
        xlim([0 inf])
end
%==========================================================================   
% Figure: EKF Y-Velocity Error
if EKF_yd_err_plot == 1

    figure(fig_id); fig_id = fig_id + 1;
        plot(time_plot, errors_ydot) 
        hold on
        plot(time_plot, 3*sig_ydot) 
        plot(time_plot, -3*sig_ydot) 
        legend('EKF Error', '3-\sigma Bounds')
        title('EKF Y-Velocity Error')
        ylabel('Y Velocity (m/s)')
        set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
        grid on
        xlim([0 inf])
end
%==========================================================================  
% Figure: EKF Z-Velocity Error
if EKF_zd_err_plot == 1

    figure(fig_id); fig_id = fig_id + 1;
        plot(time_plot, errors_zdot) 
        hold on
        plot(time_plot, 3*sig_zdot) 
        plot(time_plot, -3*sig_zdot) 
        legend('EKF Error', '3-\sigma Bounds')
        title('EKF Z-Velocity Error')
        ylabel('Z Velocity (m/s)')
        set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
        grid on
        xlim([0 inf])
end
%==========================================================================  
% Figure: EKF Theta_dot Error
if EKF_thetad_err_plot == 1

    figure(fig_id); fig_id = fig_id + 1;
        plot(time_plot, errors_theta_dot) 
        hold on
        plot(time_plot, 3*sig_theta_dot) 
        plot(time_plot, -3*sig_theta_dot) 
        legend('EKF Error', '3-\sigma Bounds')
        title('EKF \theta dot Error')
        ylabel('\theta dot(rad/s)')
        set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
        grid on
        xlim([0 inf])
end
%==========================================================================  
% Figure: EKF rt_dot Error
if EKF_rtd_err_plot == 1

    figure(fig_id); fig_id = fig_id + 1;
        plot(time_plot, errors_rt_dot) 
        hold on
        plot(time_plot, 3*sig_rt_dot)%*10^-4) 
        plot(time_plot, -3*sig_rt_dot)%*10^-4) 
        legend('EKF Error', '3-\sigma Bounds')
        title('EKF r_t dot Error')
        ylabel('r_t dot (m/s)')
        set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
        grid on
        xlim([0 inf])
end
    
%==========================================================================  
% Figure: EKF State Covariances
%{
    figure(fig_id); fig_id = fig_id + 1;
    clf
    num = 11;
    plot(time, EKF_output(:, num), time, EKF_output(:, num+11), time, EKF_output(:, num+22),...
        time, EKF_output(:, num+33), time, EKF_output(:, num+44), time, EKF_output(:, num+55),...
        time, EKF_output(:, num+66), time, EKF_output(:, num+77),time, EKF_output(:, num+88), time, EKF_output(:, num+99))
    title('State Covariances vs Time')
    xlab = xlabel('Time (hr)');
    ylabel('State Covariance')
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    legend('x', 'y', 'z', '\theta', 'rt', 'x_d_o_t', 'y_d_o_t', 'z_d_o_t', '\theta_d_o_t', 'rt_d_o_t')
    grid on  
    axis([0 time(length(time)) -20 20])
    
% 
    figure(fig_id); fig_id = fig_id + 1;
    plot(time, EKF_output(:, 55))
    title('r_t Covariance vs Time')
    xlab = xlabel('Time (hr)');
    ylabel('r_t Covariance')
    %set(gca,'FontName', 'Times')
    grid on  
 %} 
%==========================================================================  
% Figure: Kalman Gains
%{
    figure(fig_id); fig_id = fig_id + 1;
    clf
    num = 119; %Index for first column of output data cxontaining gain data
    plot(time, EKF_output(:, num), time, EKF_output(:, num+1), time, EKF_output(:, num+2),...
        time, EKF_output(:, num+3), time, EKF_output(:, num+4), time, EKF_output(:, num+5),...
        time, EKF_output(:, num+6), time, EKF_output(:, num+7))
    set(gca,'FontSize',AT_size);
    title('X-Kalman Gains vs Time')
    xlab = xlabel('Time (hr)');
    ylabel('Kalman Gain')
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    legend('x', 'y', 'z', '\theta', 'rt', 'x_d_o_t', 'y_d_o_t', 'z_d_o_t')
    grid on  

    figure(17)
    clf
    num = 127;
    plot(time, EKF_output(:, num), time, EKF_output(:, num+1), time, EKF_output(:, num+2),...
        time, EKF_output(:, num+3), time, EKF_output(:, num+4), time, EKF_output(:, num+5),...
        time, EKF_output(:, num+6), time, EKF_output(:, num+7))
    set(gca,'FontSize',AT_size);
    title('Y-Kalman Gains vs Time')
    xlab = xlabel('Time (hr)')
    ylabel('Kalman Gain')
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    legend('x', 'y', 'z', '\theta', 'rt', 'x_d_o_t', 'y_d_o_t', 'z_d_o_t')
    grid on  
    
    figure(fig_id); fig_id = fig_id + 1;
    clf
    num = 135;
    plot(time, EKF_output(:, num), time, EKF_output(:, num+1), time, EKF_output(:, num+2),...
        time, EKF_output(:, num+3), time, EKF_output(:, num+4), time, EKF_output(:, num+5),...
        time, EKF_output(:, num+6), time, EKF_output(:, num+7))
    set(gca,'FontSize',AT_size);
    title('Z-Kalman Gains vs Time')
    xlab = xlabel('Time (hr)')
    ylabel('Kalman Gain')
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    legend('x', 'y', 'z', '\theta', 'rt', 'x_d_o_t', 'y_d_o_t', 'z_d_o_t')
    grid on  
    
    figure(fig_id); fig_id = fig_id + 1;
    clf
    num = 143;
    plot(time, EKF_output(:, num), time, EKF_output(:, num+1), time, EKF_output(:, num+2),...
        time, EKF_output(:, num+3), time, EKF_output(:, num+4), time, EKF_output(:, num+5),...
        time, EKF_output(:, num+6), time, EKF_output(:, num+7))
    set(gca,'FontSize',AT_size);
    title('\theta-Kalman Gains vs Time')
    xlab = xlabel('Time (hr)')
    ylabel('Kalman Gain')
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    legend('x', 'y', 'z', '\theta', 'rt', 'x_d_o_t', 'y_d_o_t', 'z_d_o_t')
    grid on  
    
    figure(fig_id); fig_id = fig_id + 1;
    clf
    num = 151;
    plot(time, EKF_output(:, num), time, EKF_output(:, num+1), time, EKF_output(:, num+2),...
        time, EKF_output(:, num+3), time, EKF_output(:, num+4), time, EKF_output(:, num+5),...
        time, EKF_output(:, num+6), time, EKF_output(:, num+7))
    set(gca,'FontSize',AT_size);
    title('r_t-Kalman Gains vs Time')
    xlab = xlabel('Time (hrs)');
    ylabel('Kalman Gain')
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    legend('x', 'y', 'z', '\theta', 'rt', 'x_d_o_t', 'y_d_o_t', 'z_d_o_t')
    grid on  
    
    figure(fig_id); fig_id = fig_id + 1;
    clf
    num = 159;
    plot(time, EKF_output(:, num), time, EKF_output(:, num+1), time, EKF_output(:, num+2),...
        time, EKF_output(:, num+3), time, EKF_output(:, num+4), time, EKF_output(:, num+5),...
        time, EKF_output(:, num+6), time, EKF_output(:, num+7))
    set(gca,'FontSize',AT_size);
    title('x_d_o_t-Kalman Gains vs Time')
    xlab = xlabel('Time (hrs)');
    ylabel('Kalman Gain')
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    legend('x', 'y', 'z', '\theta', 'rt', 'x_d_o_t', 'y_d_o_t', 'z_d_o_t')
    grid on  
    
    figure(fig_id); fig_id = fig_id + 1;
    clf
    num = 167;
    plot(time, EKF_output(:, num), time, EKF_output(:, num+1), time, EKF_output(:, num+2),...
        time, EKF_output(:, num+3), time, EKF_output(:, num+4), time, EKF_output(:, num+5),...
        time, EKF_output(:, num+6), time, EKF_output(:, num+7))
    set(gca,'FontSize',AT_size);
    title('y_d_o_t-Kalman Gains vs Time')
    xlab = xlabel('Time (hrs)');
    ylabel('Kalman Gain')
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    legend('x', 'y', 'z', '\theta', 'rt', 'x_d_o_t', 'y_d_o_t', 'z_d_o_t')
    grid on  
    
    figure(fig_id); fig_id = fig_id + 1;
    clf
    num = 175;
    plot(time, EKF_output(:, num), time, EKF_output(:, num+1), time, EKF_output(:, num+2),...
        time, EKF_output(:, num+3), time, EKF_output(:, num+4), time, EKF_output(:, num+5),...
        time, EKF_output(:, num+6), time, EKF_output(:, num+7))
    set(gca,'FontSize',AT_size);
    title('z_d_o_t-Kalman Gains vs Time')
    xlab = xlabel('Time (hrs)');
    ylabel('Kalman Gain')
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    legend('x', 'y', 'z', '\theta', 'rt', 'x_d_o_t', 'y_d_o_t', 'z_d_o_t')
    grid on  
    
    figure(fig_id); fig_id = fig_id + 1;
    clf
    num = 183;
    plot(time, EKF_output(:, num), time, EKF_output(:, num+1), time, EKF_output(:, num+2),...
        time, EKF_output(:, num+3), time, EKF_output(:, num+4), time, EKF_output(:, num+5),...
        time, EKF_output(:, num+6), time, EKF_output(:, num+7))
    set(gca,'FontSize',AT_size);
    title('\theta_d_o_t-Kalman Gains vs Time')
    xlab = xlabel('Time (hrs)');
    ylabel('Kalman Gain')
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    legend('x', 'y', 'z', '\theta', 'rt', 'x_d_o_t', 'y_d_o_t', 'z_d_o_t')
    grid on  
    
    figure(fig_id); fig_id = fig_id + 1;
    clf
    num = 191;
    plot(time, EKF_output(:, num), time, EKF_output(:, num+1), time, EKF_output(:, num+2),...
        time, EKF_output(:, num+3), time, EKF_output(:, num+4), time, EKF_output(:, num+5),...
        time, EKF_output(:, num+6), time, EKF_output(:, num+7))
    set(gca,'FontSize',AT_size);
    title('r_t_d_o_t-Kalman Gains vs Time')
    xlab = xlabel('Time (hrs)');
    ylabel('Kalman Gain')
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    legend('x', 'y', 'z', '\theta', 'rt', 'x_d_o_t', 'y_d_o_t', 'z_d_o_t')
    grid on  
    
    figure(fig_id); fig_id = fig_id + 1;
    clf
    num = 199;
    plot(time, EKF_output(:, num), time, EKF_output(:, num+1), time, EKF_output(:, num+2),...
        time, EKF_output(:, num+3), time, EKF_output(:, num+4), time, EKF_output(:, num+5),...
        time, EKF_output(:, num+6), time, EKF_output(:, num+7))
    set(gca,'FontSize',AT_size);
    title('EKF Residuals vs Time')
    xlab = xlabel('Time (hrs)');
    ylabel('Residual')
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    legend('x', 'y', 'z', '\theta', 'rt', 'x_d_o_t', 'y_d_o_t', 'z_d_o_t')
    grid on  
    
    figure(fig_id); fig_id = fig_id + 1;
    clf
    num = 11;
    plot(time, EKF_output(:, num), time, EKF_output(:, num+11), time, EKF_output(:, num+22),...
        time, EKF_output(:, num+33), '--', time, EKF_output(:, num+44), time, EKF_output(:, num+55),...
        time, EKF_output(:, num+66), time, EKF_output(:, num+77), time, EKF_output(:, num+88), ...
        time, EKF_output(:, num+99))
    set(gca,'FontSize',AT_size);
    title('EKF a Posteriori State Covariance Diagonal Elements vs Time')
    xlab = xlabel('Time (hrs)');
    ylabel('Diagonal Covariance Element Value')
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    axis([0 time(length(time)) 0 0.1])
    legend('x', 'y', 'z', '\theta',' r_t', 'x_d_o_t', 'y_d_o_t', 'z_d_o_t', '\theta_d_o_t', 'rt_d_o_t')
    grid on  
    %}
    
%==========================================================================  
%% Animation: Relative Positions
if RelativePosition_Animation == 1
    figure(fig_id); fig_id = fig_id + 1;
    title('Relative Positions', ... 
        'Fontweight', 'bold', 'FontSize', 12)
    curve1 = animatedline('color','r', 'LineWidth', 0.5);
    curve2 = animatedline('color','b', 'LineWidth', 0.5);
    curve3 = animatedline('color','g', 'LineWidth', 0.5);
    %set(gca, 'XLim', [-200 200], 'YLim', ...
    %    [-200 200], 'ZLim', [-2*10^-7 2*10^-7]); 
    xlab = xlabel(' x (m) ');
    ylabel(' y (m) ')
    zlabel(' z (m) ')
    view(-40, 30);
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    grid on
    hold on;
    EKF  = scatter3(EKF_output(1,1), EKF_output(1,2), EKF_output(1,3), 'MarkerEdgeColor', 'r');
    NERM  = scatter3(NERM_output(1,1), NERM_output(1,2), NERM_output(1,3), 'MarkerEdgeColor', 'b');
    GPS  = scatter3(GPS_output(1,1), GPS_output(1,2), GPS_output(1,3), 'MarkerEdgeColor', 'b');
    legend('EKF', 'Dynamics', 'GPS')
    for i = 1:100:plot_end
        addpoints(curve1, EKF_output(i,1), EKF_output(i,2), EKF_output(i,3));
        addpoints(curve2, NERM_output(i,1), NERM_output(i,2), NERM_output(i,3));
        addpoints(curve3, GPS_output(i,1), GPS_output(i,2), GPS_output(i,3));
        EKF  = scatter3(EKF_output(i,1), EKF_output(i,2), EKF_output(i,3), 'MarkerEdgeColor', 'r');
        NERM  = scatter3(NERM_output(i,1), NERM_output(i,2), NERM_output(i,3), 'MarkerEdgeColor', 'b');
        GPS  = scatter3(GPS_output(i,1), GPS_output(i,2), GPS_output(i,3), 'MarkerEdgeColor', 'b');
        drawnow
        delete(EKF);
        delete(GPS);
        delete(NERM);
    end
end
%==========================================================================  
%% Animation: Relative Velocities

if RelativeVelocity_Animation == 1

    figure(fig_id); fig_id = fig_id + 1;
    title('Relative Velocities Velocities', ... 
        'Fontweight', 'bold', 'FontSize', 12)
    curve1 = animatedline('color','r', 'LineWidth', 0.5);
    curve2 = animatedline('color','b', 'LineWidth', 0.5);
    curve3 = animatedline('color','g', 'LineWidth', 0.5);
    %set(gca, 'XLim', [-.2 .2], 'YLim', ...
    %    [-.2 .2], 'ZLim', [-2*10^-7 2*10^-7]); 
    xlab = xlabel(' x Velocity (m/s) ');
    ylabel(' y Velocity (m/s) ')
    zlabel(' z Velocity (m/s) ')
    view(-40, 30);
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    grid on
    hold on;
    EKF  = scatter3(EKF_output(1,4), EKF_output(1,5), EKF_output(1,6), 'MarkerEdgeColor', 'r');
    NERM  = scatter3(NERM_output(1,4), NERM_output(1,5), NERM_output(1,6), 'MarkerEdgeColor', 'b');
    GPS  = scatter3(GPS_output(1,4), GPS_output(1,5), GPS_output(1,6), 'MarkerEdgeColor', 'b');
    legend('EKF', 'Dynamics', 'Measurements')
    for i = 1:100:30000;
        addpoints(curve1, EKF_output(i,4), EKF_output(i,5), EKF_output(i,6));
        addpoints(curve2, NERM_output(i,4), NERM_output(i,5), NERM_output(i,6));
        addpoints(curve3, GPS_output(i,4), GPS_output(i,5), GPS_output(i,6));
        EKF  = scatter3(EKF_output(i,4), EKF_output(i,5), EKF_output(i,6), 'MarkerEdgeColor', 'r');
        NERM  = scatter3(NERM_output(i,4), NERM_output(i,5), NERM_output(i,6), 'MarkerEdgeColor', 'b');
        GPS  = scatter3(GPS_output(i,4), GPS_output(i,5), GPS_output(i,6), 'MarkerEdgeColor', 'b');
        drawnow
        delete(EKF);
        delete(GPS);
        delete(NERM);
    end
end

%==========================================================================  
%% Animation: 3D ECI & LVLH Positions w/ Earth

if RelativeMotion_Movie == 1
    
    video = VideoWriter([pwd '\animations\RelativeMotionMonday.mpeg']); 
    video.FrameRate = 10;
    open(video)

    figure(fig_id); fig_id = fig_id + 1;
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    
    %Main Figure - ECI Frame
    lim_ECI = 1.25*max(max(abs(r_target_ECI)));
    sub1 = subplot(1,2,1);
    title('ECI Positions')
    hold on;
    set(gca, 'XLim', 1.1*[-lim_ECI lim_ECI], 'YLim', ...
        1.1*[-lim_ECI lim_ECI], 'ZLim', ...
        1.1*[-lim_ECI lim_ECI]); 

    %Plot Earth
    axes(sub1);
    equat_rad = Constants.R_Earth;
    polar_rad = 6356752.3142;
    [xx yy zz]=ellipsoid (0,0,0,equat_rad, equat_rad, polar_rad);
    load('topo.mat','topo','topomap1');
    topo2 = [topo(:,181:360) topo(:,1:180)];
    pro.FaceColor= 'texture';
    pro.EdgeColor = 'none';
    pro.FaceLighting = 'phong';     
    pro.Cdata = topo2;
    earth= surface(xx,yy,zz,pro,'HandleVisibility','off');
    colormap(topomap1) 

    line([0 lim_ECI],[0 0],[0 0],'Color', 'black', 'Marker','.','LineStyle','-','LineWidth', 2.5, 'HandleVisibility','off')
    line([0 0],[0 lim_ECI],[0 0],'Color', 'black', 'Marker','.','LineStyle','-','LineWidth', 2.5, 'HandleVisibility','off')
    line([0 0],[0 0],[0 lim_ECI],'Color', 'black', 'Marker','.','LineStyle','-','LineWidth', 2.5, 'HandleVisibility','off')

    curve1 = animatedline('color','k', 'LineWidth', 1.5, 'LineStyle', '-');
    curve2 = animatedline('color','r', 'LineWidth', 1.5, 'LineStyle', '-.');
    xlab = xlabel(' X (km) ');
    ylabel(' Y (km) ')
    zlabel(' Z (km) ')
    %legend('Target', 'Chaser', 'Location', 'southeast')
    llab = legend({'Target', 'Chaser'},'orientation','vertical', 'box','off','color','none',...
        'Position',[0.75 0.178 0.1 0.08], 'Interpreter', 'Latex');
    view(130, 30);
    Target  = scatter3(r_target_ECI(1,1), r_target_ECI(1,2), r_target_ECI(1,3), 5,'filled', 'MarkerFaceColor', 'k', 'LineWidth', 2.5, 'MarkerEdgeColor', 'k', 'HandleVisibility','off');
    Chaser  = scatter3(r_chaser_ECI(1,1), r_chaser_ECI(1,2), r_chaser_ECI(1,3), 100, 'LineWidth', 2.5,'MarkerEdgeColor', 'r', 'HandleVisibility','off');
    sub1.Color = 'none';

    set(gca,'xticklabel',num2str(get(gca,'xtick')'/1000))
    set(gca,'yticklabel',num2str(get(gca,'ytick')'/1000))
    set(gca,'zticklabel',num2str(get(gca,'ztick')'/1000))   
    
    %Secondary Figure - LVLH Frame
    lim_LVLH = 1.1*max(max(abs(r_clean)));
    sub2 = subplot(1,2,2);
    title('LVLH Positions', ... 
        'Fontweight', 'bold', 'FontSize', 12)
    set(gca, 'XLim', [-lim_LVLH lim_LVLH], 'YLim', ...
        [-lim_LVLH lim_LVLH], 'ZLim', ...
        [-lim_LVLH lim_LVLH]); 
    curve3 = animatedline('color','r', 'LineWidth', 0.5);
    sub2.Position = [ 0.6 0.5 0.35 0.45];
    sub2.Color = 'none';
    xlab = xlabel(' x (m) ');
    ylabel(' y (m) ')
    zlabel(' z (m) ')
    hold on;
    grid on
    box on
    %view(45, 45);
    view(145, 15);
    Target_LVLH  = scatter3(0, 0, 0, 'filled', 'MarkerFaceColor', 'k','HandleVisibility','off');
    Chaser_LVLH  = scatter3(r_clean(1,1), r_clean(1,2), r_clean(1,3), 'MarkerEdgeColor', 'r');

    sub1.Position = [ 0.1 0.1 0.6 0.85];    %Readjust placement of main figure

    for i = 1:100:plot_end
        axes(sub1);
        delete(Target);
        delete(Chaser);
        addpoints(curve1, r_target_ECI(i,1), r_target_ECI(i,2), r_target_ECI(i,3));
        addpoints(curve2, r_chaser_ECI(i,1), r_chaser_ECI(i,2), r_chaser_ECI(i,3));      

        Target  = scatter3(r_target_ECI(i,1), r_target_ECI(i,2), r_target_ECI(i,3), 5,'filled', 'MarkerFaceColor', 'k', 'LineWidth', 2.5, 'MarkerEdgeColor', 'k', 'HandleVisibility','off');
        Chaser  = scatter3(r_chaser_ECI(i,1), r_chaser_ECI(i,2), r_chaser_ECI(i,3), 100, 'LineWidth', 2.5,'MarkerEdgeColor', 'r', 'HandleVisibility','off');

        axes(sub2);
        delete(Chaser_LVLH);
        Chaser_LVLH  = scatter3(r_clean(i,1), r_clean(i,2), r_clean(i,3), 'MarkerEdgeColor', 'r');
        addpoints(curve3, r_clean(i,1), r_clean(i,2), r_clean(i,3));   

        drawnow
        f = getframe(gcf);
        writeVideo(video, f);    

    end
    close(video)
end

%==========================================================================  
%% Bar Graphs: Relative Position and Velocity RMS
if RMSbar_plot == 1
    filter_names = {'GPS', 'NERM EKF'};
    RMS_pos_values = [GPS_RMSE_Avg_pos EKF_RMSE_Avg_pos]; 
    RMS_vel_values = [GPS_RMSE_Avg_vel EKF_RMSE_Avg_vel]; 

    figure(fig_id); fig_id = fig_id + 1;
    bar(RMS_pos_values)
    set(gca,'xticklabel',filter_names);
    ylabel('Position Error RMS (m)')
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);

    figure(fig_id); fig_id = fig_id + 1;
    bar(RMS_vel_values)
    set(gca,'xticklabel',filter_names);
    ylabel('Velocity Error RMS (m/s)')
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
end

%==========================================================================
% Plots for Adaptive EKFs =================================================
if EKF_flag ~= 0 
%% Figure w/ Subfigures: Q Covariance Plots
    if (Q_cov_subplot)
        figure(fig_id); fig_id = fig_id + 1;
        set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
        subplot(5,2,1)
        %subplot(3,1,1)
        plot(time_plot(plot_start:plot_end), squeeze(Q_output(1,1,plot_start:plot_end)))
        set(gca,'xticklabel',[])
        ylabel('q_{11}')
        grid on
        xlim([0 inf])
        hold on


        subplot(5,2,2)
        plot(time_plot(plot_start:plot_end), squeeze(Q_output(6,6,plot_start:plot_end)))
        ylabel('q_{66}')
        set(gca,'xticklabel',[])
        grid on
        xlim([0 inf])
        hold on


        subplot(5,2,3)
        %subplot(3,1,2)
        plot(time_plot(plot_start:plot_end), squeeze(Q_output(2,2,plot_start:plot_end)))
        ylabel('q_{22}')
        set(gca,'xticklabel',[])
        grid on
        xlim([0 inf])
        hold on


        subplot(5,2,4)
        plot(time_plot(plot_start:plot_end), squeeze(Q_output(7,7,plot_start:plot_end)))
        ylabel('q_{77}')
        set(gca,'xticklabel',[])
        grid on
        xlim([0 inf])
        hold on

        subplot(5,2,5)
        %subplot(3,1,3)
        plot(time_plot(plot_start:plot_end), squeeze(Q_output(3,3,plot_start:plot_end)))
        ylabel('q_{33}')
        set(gca,'xticklabel',[])
        grid on
        xlim([0 inf])
        hold on

        subplot(5,2,6)
        plot(time_plot(plot_start:plot_end), squeeze(Q_output(8,8,plot_start:plot_end)))
        ylabel('q_{88}')
        set(gca,'xticklabel',[])
        grid on
        xlim([0 inf])
        hold on

        subplot(5,2,7)
        plot(time_plot(plot_start:plot_end), squeeze(Q_output(4,4,plot_start:plot_end)))
        ylabel('q_{44}')
        set(gca,'xticklabel',[])
        grid on
        xlim([0 inf])
        hold on

        subplot(5,2,8)
        plot(time_plot(plot_start:plot_end), squeeze(Q_output(9,9,plot_start:plot_end)))
        ylabel('q_{99}')
        set(gca,'xticklabel',[])
        grid on
        xlim([0 inf])
        hold on

        subplot(5,2,9)
        plot(time_plot(plot_start:plot_end), squeeze(Q_output(5,5,plot_start:plot_end)))
        xlab = xlabel('Time (hrs)');
        ylabel('q_{55}')
        grid on
        xlim([0 inf])
        hold on

        subplot(5,2,10)
        plot(time_plot(plot_start:plot_end), squeeze(Q_output(10,10,plot_start:plot_end)))
        xlab = xlabel('Time (hrs)');
        ylabel('q_{10,10}')
        grid on
        xlim([0 inf])
        hold on
    end
    
%% ========================================================================
% Figure w/ Subfigures: R Covariance Plots

    if (R_cov_subplot)

        figure(fig_id); fig_id = fig_id + 1;
        set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
        subplot(4,2,1)
        %subplot(3,1,1)
        plot(time_plot(plot_start:plot_end), squeeze(R_output(1,1,plot_start:plot_end)))
        set(gca,'xticklabel',[])
        ylabel('R_{11}')
        grid on
        xlim([0 inf])
        hold on
        %plot(time(plot_start:plot_end), GPS_error(plot_start:plot_end,1).^2)
        %plot(time(plot_start:plot_end), 3*GPS_std_sample(plot_start:plot_end,1))
        %plot(time(plot_start:plot_end), -3*GPS_std_sample(plot_start:plot_end,1))

        subplot(4,2,2)
        plot(time_plot(plot_start:plot_end), squeeze(R_output(5,5,plot_start:plot_end)))
        ylabel('R_{55}')
        set(gca,'xticklabel',[])
        grid on
        xlim([0 inf])
        hold on
        %plot(time(plot_start:plot_end), 3*GPS_std_sample(plot_start:plot_end,4))
        %}

        subplot(4,2,3)
        %subplot(3,1,2)
        plot(time_plot(plot_start:plot_end), squeeze(R_output(2,2,plot_start:plot_end)))
        ylabel('R_{22}')
        set(gca,'xticklabel',[])
        grid on
        xlim([0 inf])
        hold on
        %plot(time(plot_start:plot_end), 3*GPS_std_sample(plot_start:plot_end,2))


        subplot(4,2,4)
        plot(time_plot(plot_start:plot_end), squeeze(R_output(6,6,plot_start:plot_end)))
        ylabel('R_{66}')
        set(gca,'xticklabel',[])
        grid on
        xlim([0 inf])
        hold on
        %plot(time(plot_start:plot_end), 3*GPS_std_sample(plot_start:plot_end,5))
        %}

        subplot(4,2,5)
        %subplot(3,1,3)
        plot(time_plot(plot_start:plot_end), squeeze(R_output(3,3,plot_start:plot_end)))
        ylabel('R_{33}')
        set(gca,'xticklabel',[])
        grid on
        xlim([0 inf])
        hold on
        %plot(time(plot_start:plot_end), 3*GPS_std_sample(plot_start:plot_end,3))

        subplot(4,2,6)
        plot(time_plot(plot_start:plot_end), squeeze(R_output(7,7,plot_start:plot_end)))
        ylabel('R_{77}')
        set(gca,'xticklabel',[])
        grid on
        xlim([0 inf])
        hold on
        %plot(time(plot_start:plot_end), 3*GPS_std_sample(plot_start:plot_end,6))

        subplot(4,2,7)
        plot(time_plot(plot_start:plot_end), squeeze(R_output(4,4,plot_start:plot_end)))
        ylabel('R_{44}')
        set(gca,'xticklabel',[])
        grid on
        xlim([0 inf])
        hold on
    end
%% ========================================================================
% Figure w/ Subfigures: Observed/Smoothed Residual Covariances
%{
figure(fig_id); fig_id = fig_id + 1;
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);

    %subplot(4,2,1)
    subplot(3,1,1)
    plot(time(plot_start:plot_end), Pv_obs(plot_start:plot_end, 1), ...
        '--k', 'LineWidth', 0.2)
    title('Observed and Smoothed Residual Covariances')
    hold on
    plot(time(plot_start:plot_end), Pv_smoothed(plot_start:plot_end, 1), ...
        '-k', 'LineWidth', 0.5)
    ylabel('x')
    set(gca,'xticklabel',[])
    grid on
    axis([time(plot_start) time(plot_end) -1.1*max(max(abs(Pv_obs(2:end)))) 1.1*max(max(abs(Pv_obs(2:end))))]);

    %{
    subplot(4,2,2)
    plot(time(plot_start:plot_end), Pv_obs(plot_start:plot_end, 5),...
        '--y', 'LineWidth', 0.05)
    hold on
    plot(time(plot_start:plot_end), Pv_smoothed(plot_start:plot_end, 5),...
        '-k', 'LineWidth', 0.5)
    ylabel('x_d_o_t')
    set(gca,'xticklabel',[])
    grid on
    %}
    
    %subplot(4,2,3)
    subplot(3,1,2)
    plot(time(plot_start:plot_end), Pv_obs(plot_start:plot_end, 2),...
        '--k', 'LineWidth', 0.2)
    hold on
    set(gca,'xticklabel',[])
    plot(time(plot_start:plot_end), Pv_smoothed(plot_start:plot_end, 2),...
        '-k', 'LineWidth', 0.5)
    ylabel('y')
    grid on
    axis([time(plot_start) time(plot_end) -1.1*max(max(abs(Pv_obs))) 1.1*max(max(abs(Pv_obs)))]);
    %{
    subplot(4,2,4)
    plot(time(plot_start:plot_end), Pv_obs(plot_start:plot_end, 6),...
        '--y', 'LineWidth', 0.05)
    hold on
    plot(time(plot_start:plot_end), Pv_smoothed(plot_start:plot_end, 6), ...
        '-k', 'LineWidth', 0.5)
    set(gca,'xticklabel',[])
    ylabel('y_d_o_t')
    grid on
    %}
    
    %subplot(4,2,5)
    subplot(3,1,3)
    plot(time(plot_start:plot_end), Pv_obs(plot_start:plot_end, 3),...
        '--k', 'LineWidth', 0.2)
    hold on
    plot(time(plot_start:plot_end), Pv_smoothed(plot_start:plot_end, 3), ...
        '-k', 'LineWidth', 0.5)
    set(gca,'xticklabel',[])
    ylabel('z')
    grid on
    axis([time(plot_start) time(plot_end) -1.1*max(max(abs(Pv_obs))) 1.1*max(max(abs(Pv_obs)))]);
    
    %{
    subplot(4,2,6)
    plot(time(plot_start:plot_end), Pv_obs(plot_start:plot_end, 7),...
        '--y', 'LineWidth', 0.05)
    hold on
    plot(time(plot_start:plot_end), Pv_smoothed(plot_start:plot_end, 7), ...
        '-k', 'LineWidth', 0.5)
    ylabel('z_d_o_t')
    set(gca,'xticklabel',[])
    grid on

    subplot(4,2,7)
    plot(time(plot_start:plot_end), Pv_obs(plot_start:plot_end, 4),...
        '--y', 'LineWidth', 0.05)
    hold on
    plot(time(plot_start:plot_end), Pv_smoothed(plot_start:plot_end, 4), ...
        '-k', 'LineWidth', 0.5)
   % xlab = xlabel('Time (hrs)');
   set(gca,'xticklabel',[])
    ylabel('\theta')
    grid on
    %}
    legend('Smoothed', 'Raw')
%}   

%% ========================================================================
% Figure w/ Subfigures: Theoretical Residual Covariances
%{
figure(fig_id); fig_id = fig_id + 1;
    subplot(4,2,1)
    plot(time(plot_start:plot_end), Pr_theo(plot_start:plot_end, 1))
    %title('Theoretical Residual Covariances')
    hold on
    ylabel('x')
    set(gca,'xticklabel',[])
    grid on
 
    subplot(4,2,2)
    plot(time(plot_start:plot_end), Pr_theo(plot_start:plot_end, 5))
    hold on
    ylabel('x_d_o_t')
    set(gca,'xticklabel',[])
    grid on
    
    subplot(4,2,3)
    plot(time(plot_start:plot_end), Pr_theo(plot_start:plot_end, 2))
    hold on
    set(gca,'xticklabel',[])
    ylabel('y')
    grid on
 
    subplot(4,2,4)
    plot(time(plot_start:plot_end), Pr_theo(plot_start:plot_end, 6))
    hold on
    set(gca,'xticklabel',[])
    ylabel('y_d_o_t')
    grid on
 
    subplot(4,2,5)
    plot(time(plot_start:plot_end), Pr_theo(plot_start:plot_end, 3))
    hold on
    set(gca,'xticklabel',[])
    ylabel('z')
    grid on
 
    subplot(4,2,6)
    plot(time(plot_start:plot_end), Pr_theo(plot_start:plot_end, 7))
    hold on
    ylabel('z_d_o_t')
    set(gca,'xticklabel',[])
    grid on

    subplot(4,2,7)
    plot(time(plot_start:plot_end), Pr_theo(plot_start:plot_end, 4))
    hold on
    %xlab = xlabel('Time (hrs)');
    ylabel('\theta')
    grid on
%}

%% ========================================================================
% Figure w/ Subfigures: Residual Covariance Error Signals (Theoretical - Observed)
if EKF_flag == 1
    DOM_subplot = 0;
end
if (DOM_subplot)

    figure(fig_id); fig_id = fig_id + 1;
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    subplot(4,2,1)
    %subplot(3,1,1)
    plot(time_plot(plot_start:plot_end), P_error(plot_start:plot_end, 1))
    hold on
    xlim([0 inf])
    ylabel('DOM x')
    set(gca,'xticklabel',[])
    grid on
    %axis([0 inf -100 100])
    
    subplot(4,2,2)
    plot(time_plot(plot_start:plot_end), P_error(plot_start:plot_end, 5))
    hold on
    ylabel('DOM x_d_o_t')
    set(gca,'xticklabel',[])
    grid on
    xlim([0 inf])
    %axis([0 inf -1 1])
    
    subplot(4,2,3)
    %subplot(3,1,2)
    plot(time_plot(plot_start:plot_end), P_error(plot_start:plot_end, 2))
    hold on
    set(gca,'xticklabel',[])
    ylabel('DOM y')
    grid on
    xlim([0 inf])
    %axis([0 inf -100 100])
    
    subplot(4,2,4)
    plot(time_plot(plot_start:plot_end), P_error(plot_start:plot_end, 6))
    hold on
    set(gca,'xticklabel',[])
    ylabel('DOM y_d_o_t')
    grid on
    xlim([0 inf])
    %axis([0 inf -1 1])
    
    subplot(4,2,5)
    %subplot(3,1,3)
    plot(time_plot(plot_start:plot_end), P_error(plot_start:plot_end, 3))
    hold on
    set(gca,'xticklabel',[])
    ylabel('DOM z')
    grid on
    xlim([0 inf])
    %axis([0 inf -100 100])
    
    subplot(4,2,6)
    plot(time_plot(plot_start:plot_end), P_error(plot_start:plot_end, 7))
    hold on
    ylabel('DOM z_d_o_t')
    set(gca,'xticklabel',[])
    grid on
    xlim([0 inf])
    %axis([0 inf -1 1])
    
    subplot(4,2,7)
    plot(time_plot(plot_start:plot_end), P_error(plot_start:plot_end, 4))
    hold on
    %xlab = xlabel('Time (hrs)');
    ylabel('DOM \theta')
    grid on
    xlim([0 inf])
    %axis([0 inf -1 1])
end

%% ========================================================================
% Figure: Normalized Innovation Squared Plot
if EKF_flag == 2
    if (NIS_plot)
        figure(fig_id); fig_id = fig_id + 1;
            set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
            plot(time_plot(plot_start:plot_end), NIS(plot_start:plot_end))
            ylabel('NIS - Divergence Ratio')
            xlab = xlabel('Time (hrs)');
            grid on
            hold on
            xlim([0 inf])
            %axis([0 inf -100 100])
    end
end

%% ========================================================================
% Figure: True Degree of Divergence Plot
%{
figure(fig_id); fig_id = fig_id + 1;
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    plot(time(plot_start:plot_end), DOD_true(plot_start:plot_end))
    ylabel('DOD True')
    xlab = xlabel('Time (hrs)');
    grid on
    hold on
    %}

%% ========================================================================
% Figure: Difference in Degree of Divergence Plot
if EKF_flag == 2
    if (deltaDOD_plot)
        figure(fig_id); fig_id = fig_id + 1;
            set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
            plot(time_plot(plot_start:plot_end), deltaDOD(plot_start:plot_end))
            ylabel('\Delta DOD')
            xlab = xlabel('Time (hrs)');
            grid on
            xlim([0 inf])
            %axis([0 inf -100 100])
            hold on
    end
end
% =========================================================================
% Figure: Trace of the Covariances
    %{
    for i = 1:plot_end
        trace_Q(i) = trace(Q_output(:,:,i));
        trace_R(i) = trace(R_output(:,:,i));
    end
    
    figure(fig_id); fig_id = fig_id + 1;
    plot(time(plot_start:plot_end), trace_Q(plot_start:plot_end))
    hold on
    plot([0 time(plot_end)],[trace(Q0) trace(Q0)],'--')
    title('Trace of Q')
    xlab = xlabel('Time (hrs)');
    ylabel('Tr(Q)')
    legend('Adapted Q', 'Initial Q0')
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    grid on
    
    % Figure: Trace of R Covariance
    figure(fig_id); fig_id = fig_id + 1;
    plot(time(plot_start:plot_end), trace_R(plot_start:plot_end))
    hold on
    plot([0 time(plot_end)],[trace(R0) trace(R0)],'--')
    title('Trace of R')
    xlab = xlabel('Time (hrs)');
    ylabel('Tr(R)')
    legend('Adapted R', 'Initial R0')
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    grid on
    
    % Figure: Trace of P Covariance (Posteriori)
    figure(fig_id); fig_id = fig_id + 1;
    plot(time, trace_P)
    title('Trace of P^+')
    xlab = xlabel('Time (hrs)');
    ylabel('Tr(P)')
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    grid on
    %}

end

%% ======================================================== 

%Sample Position Error Plot
%{
print_plots = 'off';
plot_start = 1;
step_size = 10;
plot_end = length(time_plot);

time_plot = time_plot*60; %Convert time to minutes

figure(fig_id); fig_id = fig_id+1;
subplot(3,1,1)  
    hold on
    grid on
    ylab = ylabel('$x$ (m)');
    set(ylab,'Interpreter','latex'); 
    set(gca,'XTickLabel',[]); 
    %axis([ 0 time(plot_end) -2.5 2.5])
    set(gcf, 'units','normalized','outerposition',[0 0 0.75 1]);
    plot(time_plot(plot_start:step_size:plot_end), r_clean(plot_start:step_size:plot_end,1), 'Color','k','LineWidth',2) 
    
    subplot(3,1,2)  
    hold on
    grid on
    set(gca,'XTickLabel',[]); 
    %axis([ 0 time(plot_end) -3 3])
    ylab = ylabel('$x_m$ (m)');
    set(ylab,'Interpreter','latex'); 
    set(gcf, 'units','normalized','outerposition',[0 0 0.75 1]);
    plot(time_plot(plot_start:step_size:plot_end), GPS_output(plot_start:step_size:plot_end,1), ':', 'Color','k', 'LineWidth',2) 
    
    subplot(3,1,3)  
    hold on
    grid on
    xlab = xlabel('Time (minutes)');  
    %axis([ 0 time(plot_end) -2 2])
    ylab = ylabel('$e_z$ (m)');
    set(ylab,'Interpreter','latex'); 
    set(gcf, 'units','normalized','outerposition',[0 0 0.75 1]);
    plot(time_plot(plot_start:step_size:plot_end), GPS_output(plot_start:step_size:plot_end,1)-r_clean(plot_start:step_size:plot_end,1), ':', 'Color','k', 'LineWidth',2) 

    if strcmp(print_plots, 'on')
    file_print = fullfile(pathname, 'Fig_PositionErrorsEx');
    print(file_print,  '-depsc');
    end
%}

%% =========================================================================  
% Figures: EKF State Error Histograms
nBins = 40; %number of bins to use in histograms

% Figure: EKF X-Position Error Histogram
if EKF_x_err_histo == 1
    figure(fig_id); fig_id = fig_id + 1;
        
        %Set data
        yA = EKF_error(n_start:n_end,1);
        yB = GPS_error(n_start:n_end,1);

        % Set extrema
        max_val = ceil(max(abs([yA; yB])));
        minY    = -max_val;
        maxY    = max_val;

        % Add in extrema placeholders to adjust bins to a common scale
        yA = [yA; minY; maxY];
        yB = [yB; minY; maxY];

        % Bin data
        [countsA, binsA] = hist(yA, nBins);
        [countsB, binsB] = hist(yB, nBins);

        % Remove extrema placeholders from counts
        histEnds = zeros(size(binsA));
        histEnds(1) = 1; %removes minimum placeholder
        histEnds(end) = 1; %removes maximum placeholder

        countsA2 = countsA - histEnds;
        countsB2 = countsB - histEnds;

        % Plot histograms
        hold on     
        histogram('BinEdges', [minY binsA], 'BinCounts', countsA2)
        histogram('BinEdges', [minY binsB], 'BinCounts', countsB2)
        set(gca, 'XLim', [minY maxY])        

        tlab = title('EKF-X Error Histogram');
        llab = legend('EKF', 'GPS');
        xlab = xlabel(' Error Distribution (m)');
        ylab = ylabel(' Frequency Count');
        set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
        grid on
        set(tlab,'Interpreter','latex');
        set(xlab,'Interpreter','latex');
        set(ylab,'Interpreter','latex');
         
        %if strcmp(print_flag, 'on')
        %    print -deps EKF_ErrorHistogram_X
        %end
end 

% Figure: EKF Y-Position Error Histogram
if EKF_y_err_histo == 1
    figure(fig_id); fig_id = fig_id + 1;
        %Set data
        yA = EKF_error(n_start:n_end,2);
        yB = GPS_error(n_start:n_end,2);

        % Set extrema
        max_val = ceil(max(abs([yA; yB])));
        minY    = -max_val;
        maxY    = max_val;

        % Add in extrema placeholders to adjust bins to a common scale
        yA = [yA; minY; maxY];
        yB = [yB; minY; maxY];

        % Bin data
        [countsA, binsA] = hist(yA, nBins);
        [countsB, binsB] = hist(yB, nBins);

        % Remove extrema placeholders from counts
        histEnds = zeros(size(binsA));
        histEnds(1) = 1; %removes minimum placeholder
        histEnds(end) = 1; %removes maximum placeholder

        countsA2 = countsA - histEnds;
        countsB2 = countsB - histEnds;

        % Plot histograms
        hold on     
        histogram('BinEdges', [minY binsA], 'BinCounts', countsA2)
        histogram('BinEdges', [minY binsB], 'BinCounts', countsB2)
        set(gca, 'XLim', [minY maxY])        

        tlab = title('EKF-Y Error Histogram');
        llab = legend('EKF', 'GPS');
        xlab = xlabel(' Error Distribution (m)');
        ylab = ylabel(' Frequency Count');
        set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
        grid on
        set(tlab,'Interpreter','latex');
        set(xlab,'Interpreter','latex');
        set(ylab,'Interpreter','latex');
      
        %if strcmp(print_flag, 'on')
        %    print -deps EKF_ErrorHistogram_Y
        %end
end 

% Figure: EKF Z-Position Error Histogram
if EKF_z_err_histo == 1
    figure(fig_id); fig_id = fig_id + 1;
        %Set data
        yA = EKF_error(n_start:n_end,3);
        yB = GPS_error(n_start:n_end,3);

        % Set extrema
        max_val = ceil(max(abs([yA; yB])));
        minY    = -max_val;
        maxY    = max_val;

        % Add in extrema placeholders to adjust bins to a common scale
        yA = [yA; minY; maxY];
        yB = [yB; minY; maxY];

        % Bin data
        [countsA, binsA] = hist(yA, nBins);
        [countsB, binsB] = hist(yB, nBins);

        % Remove extrema placeholders from counts
        histEnds = zeros(size(binsA));
        histEnds(1) = 1; %removes minimum placeholder
        histEnds(end) = 1; %removes maximum placeholder

        countsA2 = countsA - histEnds;
        countsB2 = countsB - histEnds;

        % Plot histograms
        hold on     
        histogram('BinEdges', [minY binsA], 'BinCounts', countsA2)
        histogram('BinEdges', [minY binsB], 'BinCounts', countsB2)
        set(gca, 'XLim', [minY maxY])        

        tlab = title('EKF-Z Error Histogram');
        llab = legend('EKF', 'GPS');
        xlab = xlabel(' Error Distribution (m)');
        ylab = ylabel(' Frequency Count');
        set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
        grid on
        set(tlab,'Interpreter','latex');
        set(xlab,'Interpreter','latex');
        set(ylab,'Interpreter','latex');
      
        %if strcmp(print_flag, 'on')
        %    print -deps EKF_ErrorHistogram_Z
        %end
end 

% Figure: EKF X-Velocity Error Histogram
if EKF_xd_err_histo == 1
    figure(fig_id); fig_id = fig_id + 1;
        
        %Set data
        yA = 1000*EKF_error(n_start:n_end,4);
        yB = 1000*GPS_error(n_start:n_end,4);

        % Set extrema
        max_val = ceil(max(abs([yA; yB])))+10;
        minY    = -max_val;
        maxY    = max_val;

        % Add in extrema placeholders to adjust bins to a common scale
        yA = [yA; minY; maxY];
        yB = [yB; minY; maxY];

        % Bin data
        [countsA, binsA] = hist(yA, nBins);
        [countsB, binsB] = hist(yB, nBins);

        % Remove extrema placeholders from counts
        histEnds = zeros(size(binsA));
        histEnds(1) = 1; %removes minimum placeholder
        histEnds(end) = 1; %removes maximum placeholder

        countsA2 = countsA - histEnds;
        countsB2 = countsB - histEnds;

        % Plot histograms
        hold on     
        histogram('BinEdges', [minY binsA], 'BinCounts', countsA2)
        histogram('BinEdges', [minY binsB], 'BinCounts', countsB2)
        set(gca, 'XLim', [minY maxY])        

        tlab = title('EKF-Xdot Error Histogram');
        llab = legend('EKF', 'GPS');
        xlab = xlabel(' Error Distribution (mm/s)');
        ylab = ylabel(' Frequency Count');
        set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
        grid on
        set(tlab,'Interpreter','latex');
        set(xlab,'Interpreter','latex');
        set(ylab,'Interpreter','latex');
         
        %if strcmp(print_flag, 'on')
        %    print -deps EKF_ErrorHistogram_Xdot
        %end
end 

% Figure: EKF Y-Velocity Error Histogram
if EKF_yd_err_histo == 1
    figure(fig_id); fig_id = fig_id + 1;
        %Set data
        yA = 1000*EKF_error(n_start:n_end,5);
        yB = 1000*GPS_error(n_start:n_end,5);

        % Set extrema
        max_val = ceil(max(abs([yA; yB])))+10;
        minY    = -max_val;
        maxY    = max_val;

        % Add in extrema placeholders to adjust bins to a common scale
        yA = [yA; minY; maxY];
        yB = [yB; minY; maxY];

        % Bin data
        [countsA, binsA] = hist(yA, nBins);
        [countsB, binsB] = hist(yB, nBins);

        % Remove extrema placeholders from counts
        histEnds = zeros(size(binsA));
        histEnds(1) = 1; %removes minimum placeholder
        histEnds(end) = 1; %removes maximum placeholder

        countsA2 = countsA - histEnds;
        countsB2 = countsB - histEnds;

        % Plot histograms
        hold on     
        histogram('BinEdges', [minY binsA], 'BinCounts', countsA2)
        histogram('BinEdges', [minY binsB], 'BinCounts', countsB2)
        set(gca, 'XLim', [minY maxY])        

        tlab = title('EKF-Ydot Error Histogram');
        llab = legend('EKF', 'GPS');
        xlab = xlabel(' Error Distribution (mm/s)');
        ylab = ylabel(' Frequency Count');
        set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
        grid on
        set(tlab,'Interpreter','latex');
        set(xlab,'Interpreter','latex');
        set(ylab,'Interpreter','latex');
      
        %if strcmp(print_flag, 'on')
        %    print -deps EKF_ErrorHistogram_Ydot
        %end
end 

% Figure: EKF Z-Velocity Error Histogram
if EKF_zd_err_histo == 1
    figure(fig_id); fig_id = fig_id + 1;
        %Set data
        yA = 1000*EKF_error(n_start:n_end,6);
        yB = 1000*GPS_error(n_start:n_end,6);

        % Set extrema
        max_val = ceil(max(abs([yA; yB])))+10;
        minY    = -max_val;
        maxY    = max_val;

        % Add in extrema placeholders to adjust bins to a common scale
        yA = [yA; minY; maxY];
        yB = [yB; minY; maxY];

        % Bin data
        [countsA, binsA] = hist(yA, nBins);
        [countsB, binsB] = hist(yB, nBins);

        % Remove extrema placeholders from counts
        histEnds = zeros(size(binsA));
        histEnds(1) = 1; %removes minimum placeholder
        histEnds(end) = 1; %removes maximum placeholder

        countsA2 = countsA - histEnds;
        countsB2 = countsB - histEnds;

        % Plot histograms
        hold on     
        histogram('BinEdges', [minY binsA], 'BinCounts', countsA2)
        histogram('BinEdges', [minY binsB], 'BinCounts', countsB2)
        set(gca, 'XLim', [minY maxY])        

        tlab = title('EKF-Zdot Error Histogram');
        llab = legend('EKF', 'GPS');
        xlab = xlabel(' Error Distribution (mm/s)');
        ylab = ylabel(' Frequency Count');
        set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
        grid on
        set(tlab,'Interpreter','latex');
        set(xlab,'Interpreter','latex');
        set(ylab,'Interpreter','latex');
      
        %if strcmp(print_flag, 'on')
        %    print -deps EKF_ErrorHistogram_Zdot
        %end
end

%% ========================================================    
print_plots = 'off';
plot_start = 1;
step_size = 10;
plot_end = length(time_plot);

%Sample True, Measured, Errors for a state
%{
figure(fig_id); fig_id = fig_id+1;
subplot(3,1,1)  
    hold on
    grid on
    ylab = ylabel('$x$ (m)');
    set(ylab,'Interpreter','latex'); 
    set(gca,'XTickLabel',[]); 
    %axis([ 0 time(plot_end) -2.5 2.5])
    set(gcf, 'units','normalized','outerposition',[0 0 0.75 1]);
    plot(time_plot(plot_start:step_size:plot_end), r_clean(plot_start:step_size:plot_end,1), 'Color','k','LineWidth',2) 
    
    subplot(3,1,2)  
    hold on
    grid on
    set(gca,'XTickLabel',[]); 
    %axis([ 0 time(plot_end) -3 3])
    ylab = ylabel('$x_m$ (m)');
    set(ylab,'Interpreter','latex'); 
    set(gcf, 'units','normalized','outerposition',[0 0 0.75 1]);
    plot(time_plot(plot_start:step_size:plot_end), GPS_output(plot_start:step_size:plot_end,1), ':', 'Color','k', 'LineWidth',2) 
    
    subplot(3,1,3)  
    hold on
    grid on
    xlab = xlabel('Time (minutes)');  
    %axis([ 0 time(plot_end) -2 2])
    ylab = ylabel('$e_z$ (m)');
    set(ylab,'Interpreter','latex'); 
    set(gcf, 'units','normalized','outerposition',[0 0 0.75 1]);
    plot(time_plot(plot_start:step_size:plot_end), GPS_output(plot_start:step_size:plot_end,1)-r_clean(plot_start:step_size:plot_end,1), ':', 'Color','k', 'LineWidth',2) 

    if strcmp(print_plots, 'on')
    file_print = fullfile(pathname, 'Fig_PositionErrorsEx');
    print(file_print,  '-depsc');
    end   
%==========================================================================    
    %}
%% Position Errors w/ Covariance (EKF Only) - Full Simulation
%{
    figure(fig_id); fig_id = fig_id+1;
    subplot(3,1,1)  
        hold on
        grid on
        if (Cov_plots)
                fill([time_plot(plot_start:plot_end)', fliplr(time_plot(plot_start:plot_end)')],...
                    [3*sig_x', fliplr(-3*sig_x')],...
                    [0.9 0.9 0.9], 'FaceAlpha', 0.25, 'LineStyle', '-.')
        end
        ylab = ylabel('$e_x$ (m)');
        set(ylab,'Interpreter','latex'); 
        set(gca,'XTickLabel',[]); 
        %axis([ 0 time(plot_end) -2.5 2.5])
        set(gcf, 'units','normalized','outerposition',[0 0 0.75 1]);
        plot(time_plot(plot_start:step_size:plot_end), errors_x(1:step_size:end), ':', 'Color','k','LineWidth',2) 
        if (Cov_plots)
                legend('EKF 3-\sigma Bounds', 'EKF Error')
        else
                legend('EKF Error')
        end 

        subplot(3,1,2)  
        hold on
        grid on
        if (Cov_plots)
                fill([time_plot(plot_start:plot_end)', fliplr(time_plot(plot_start:plot_end)')],...
                    [3*sig_y', fliplr(-3*sig_y')],...
                    [0.9 0.9 0.9], 'FaceAlpha', 0.25, 'LineStyle', '-.')
        end
        set(gca,'XTickLabel',[]); 
        %axis([ 0 time(plot_end) -3 3])
        ylab = ylabel('$e_y$ (m)');
        set(ylab,'Interpreter','latex'); 
        set(gcf, 'units','normalized','outerposition',[0 0 0.75 1]);
        plot(time_plot(plot_start:step_size:plot_end), (errors_y(1:step_size:end)), ':', 'Color','k', 'LineWidth',2) 

        subplot(3,1,3)  
        hold on
        grid on
        if (Cov_plots)
                fill([time_plot(plot_start:plot_end)', fliplr(time_plot(plot_start:plot_end)')],...
                    [3*sig_z', fliplr(-3*sig_z')],...
                    [0.9 0.9 0.9], 'FaceAlpha', 0.25, 'LineStyle', '-.')
        end
        xlab = xlabel('Time (minutes)');  
        %axis([ 0 time(plot_end) -2 2])
        ylab = ylabel('$e_z$ (m)');
        set(ylab,'Interpreter','latex'); 
        set(gcf, 'units','normalized','outerposition',[0 0 0.75 1]);
        plot(time_plot(plot_start:step_size:plot_end), (errors_z(1:step_size:end)), ':', 'Color','k', 'LineWidth',2) 
    if strcmp(print_plots, 'on')
        %file_print = fullfile(pathname, 'Fig_PositionErrors');
        %print(file_print,  '-depsc');
    end
%}

%% Position Errors w/ Covariance Side-by-Side (EKF Only) - Full Simulation
if EKF_StatesError_pos ==1

    figure(fig_id); fig_id = fig_id+1;
    set(gcf, 'units','normalized','outerposition',[0 0 0.5 1]);

        h(1) = subplot(3,2,1);
        plot(time_plot, EKF_output(:,1),'Color','k','LineWidth',2) 
        ylab = ylabel('$\hat{x}$ (m)');    
        set(ylab,'Interpreter','latex'); 
        set(gca,'XTickLabel',[]);
        grid on
        hold on
        xlim([0 inf])

        subplot(3,2,2)  
        hold on
        grid on
        if (Cov_plots)
                fill([time_plot(plot_start:plot_end)', fliplr(time_plot(plot_start:plot_end)')],...
                    [3*sig_x', fliplr(-3*sig_x')],...
                    [0.9 0.9 0.9], 'FaceAlpha', 0.6, 'LineStyle', 'none') 
                plot(time_plot, 3*sig_x,'-.k','LineWidth',1)
                plot(time_plot,-3*sig_x,'-.k','LineWidth',1) 
        end
        ylab = ylabel('$e_x$ (m)');
        set(ylab,'Interpreter','latex'); 
        set(gca,'XTickLabel',[]); 
        xlim([0 inf])
        plot(time_plot(plot_start:step_size:plot_end), errors_x(1:step_size:end), ':', 'Color','k','LineWidth',2) 

        h(2) = subplot(3,2,3);
        plot(time_plot, EKF_output(:,2),'k') 
        set(gca,'FontSize',AT_size);
        ylab = ylabel('$\hat{y}$ (m)');    
        set(ylab,'Interpreter','latex'); 
        set(gca,'XTickLabel',[]);
        grid on
        hold on
        %axis([0 time(length(time)) -200 200])
        xlim([0 inf])

        subplot(3,2,4)  
        hold on
        grid on
        if (Cov_plots)
                fill([time_plot(plot_start:plot_end)', fliplr(time_plot(plot_start:plot_end)')],...
                    [3*sig_y', fliplr(-3*sig_y')],...
                    [0.9 0.9 0.9], 'FaceAlpha', 0.6, 'LineStyle', 'none') 
                plot(time_plot, 3*sig_y,'-.k','LineWidth',1)
                plot(time_plot,-3*sig_y,'-.k','LineWidth',1) 
        end
        set(gca,'XTickLabel',[]); 
        xlim([0 inf])    
        %axis([ 0 time(plot_end) -3 3])
        ylab = ylabel('$e_y$ (m)');
        set(ylab,'Interpreter','latex'); 
        plot(time_plot(plot_start:step_size:plot_end), (errors_y(1:step_size:end)), ':', 'Color','k') 

        h(3) = subplot(3,2,5);
        plot(time_plot, EKF_output(:,3),'k') 
        set(gca,'FontSize',AT_size);
        xlab = xlabel('Time (minutes)') ;
        set(xlab,'Interpreter','latex'); 
        ylab = ylabel('$\hat{z}$ (m)');    
        set(ylab,'Interpreter','latex'); 
        grid on
        hold on
        xlim([0 inf])

        subplot(3,2,6)  
        hold on
        grid on
        if (Cov_plots)
                fill([time_plot(plot_start:plot_end)', fliplr(time_plot(plot_start:plot_end)')],...
                    [3*sig_z', fliplr(-3*sig_z')],...
                    [0.9 0.9 0.9], 'FaceAlpha', 0.6, 'LineStyle', 'none') 
                plot(time_plot, 3*sig_z,'-.k','LineWidth',1)
                plot(time_plot,-3*sig_z,'-.k','LineWidth',1) 
        end
        %axis([ 0 time(plot_end) -2 2])
        xlim([0 inf])    
        xlab = xlabel('Time (minutes)') ;
        set(xlab,'Interpreter','latex'); 
        ylab = ylabel('$e_z$ (m)');
        set(ylab,'Interpreter','latex'); 
        plot(time_plot(plot_start:step_size:plot_end), (errors_z(1:step_size:end)), ':', 'Color','k', 'LineWidth',2) 

    if strcmp(print_plots, 'on')
        %file_print = fullfile(pathname, 'Fig_PositionErrors');
        %print(file_print,  '-depsc');
    end
end


%% Velocity Errors w/ Covariance Side-by-Side (EKF Only) - Full Simulation
if EKF_StatesError_vel ==1

    figure(fig_id); fig_id = fig_id+1;
    set(gcf, 'units','normalized','outerposition',[0 0 0.5 1]);

        h(1) = subplot(3,2,1);
        plot(time_plot, EKF_output(:,6),'k','LineWidth',2) 
        ylab = ylabel('$\hat{\dot{x}}$ (m/s)');    
        set(ylab,'Interpreter','latex'); 
        set(gca,'XTickLabel',[]);
        grid on
        hold on
        xlim([0 inf])

        subplot(3,2,2)  
        hold on
        grid on
        if (Cov_plots)
                fill([time_plot(plot_start:plot_end)', fliplr(time_plot(plot_start:plot_end)')],...
                    [3*sig_xdot', fliplr(-3*sig_xdot')],...
                    [0.9 0.9 0.9], 'FaceAlpha', 0.6, 'LineStyle', 'none') 
                plot(time_plot, 3*sig_xdot,'-.k','LineWidth',1)
                plot(time_plot,-3*sig_xdot,'-.k','LineWidth',1) 
        end
        ylab = ylabel('$e_{\dot{x}}$ (m/s)');
        set(ylab,'Interpreter','latex'); 
        set(gca,'XTickLabel',[]); 
        xlim([0 inf])
        %axis([ 0 time(plot_end) -2.5 2.5])
        plot(time_plot(plot_start:step_size:plot_end), errors_xdot(1:step_size:end), ':', 'Color','k','LineWidth',2) 

        h(2) = subplot(3,2,3);
        plot(time_plot, EKF_output(:,7),'k') 
        ylab = ylabel('$\hat{\dot{y}}$ (m/s)');    
        set(ylab,'Interpreter','latex'); 
        set(gca,'XTickLabel',[]);
        grid on
        hold on
        xlim([0 inf])

        subplot(3,2,4)  
        hold on
        grid on
        if (Cov_plots)
                fill([time_plot(plot_start:plot_end)', fliplr(time_plot(plot_start:plot_end)')],...
                    [3*sig_ydot', fliplr(-3*sig_ydot')],...
                    [0.9 0.9 0.9], 'FaceAlpha', 0.6, 'LineStyle', 'none') 
                plot(time_plot, 3*sig_ydot,'-.k','LineWidth',1)
                plot(time_plot,-3*sig_ydot,'-.k','LineWidth',1) 
        end
        set(gca,'XTickLabel',[]); 
        xlim([0 inf])    
        ylab = ylabel('$e_{\dot{y}}$ (m/s)');
        set(ylab,'Interpreter','latex'); 
        plot(time_plot(plot_start:step_size:plot_end), (errors_ydot(1:step_size:end)), ':', 'Color','k') 

        h(3) = subplot(3,2,5);
        plot(time_plot, EKF_output(:,8),'k') 
        xlab = xlabel('Time (minutes)') ;
        set(xlab,'Interpreter','latex'); 
        ylab = ylabel('$\hat{\dot{z}}$ (m/s)');    
        set(ylab,'Interpreter','latex'); 
        grid on
        hold on
        xlim([0 inf])

        subplot(3,2,6)  
        hold on
        grid on
        if (Cov_plots)
                fill([time_plot(plot_start:plot_end)', fliplr(time_plot(plot_start:plot_end)')],...
                    [3*sig_zdot', fliplr(-3*sig_zdot')],...
                    [0.9 0.9 0.9], 'FaceAlpha', 0.6, 'LineStyle', 'none') 
                plot(time_plot, 3*sig_zdot,'-.k','LineWidth',1)
                plot(time_plot,-3*sig_zdot,'-.k','LineWidth',1) 
        end
        %axis([ 0 time(plot_end) -2 2])
        xlim([0 inf])    
        xlab = xlabel('Time (minutes)') ;
        set(xlab,'Interpreter','latex'); 
        ylab = ylabel('$e_{\dot{z}}$ (m/s)');
        set(ylab,'Interpreter','latex'); 
        plot(time_plot(plot_start:step_size:plot_end), (errors_zdot(1:step_size:end)), ':', 'Color','k', 'LineWidth',2) 

    if strcmp(print_plots, 'on')
        %file_print = fullfile(pathname, 'Fig_PositionErrors');
        %print(file_print,  '-depsc');
    end
end

% =========================================================================
end
% =========================================================================