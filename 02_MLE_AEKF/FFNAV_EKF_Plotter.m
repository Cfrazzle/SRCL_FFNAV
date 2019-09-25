function [] = FFNAV_EKF_Plotter(output_file)
%% FFNAV_EKF_plots=========================================================
% Description: This function contstructs plots for the FFNAV EKF
% simulation, based on a passed .mat file that can be created by the
% Run_EKF.m program.
%
% Inputs:
%   output_file - .mat file containing simulation results of the EKF
%   time_step - Defines size of the time step
%   time_end - Defines duration of the simulation
%   make_plots - Turns plot output on/off
%   print_plots - Turns .eps plot printing on/off
%   show_RMS - Turns RMS claculations on/off
%
% Outputs:
%   Result plots (optional)
%   
% Other Functions Called:
%
% Created by: Cory Fraser - July 2016
% Last Edited: Cory Fraser - August 12, 2017
% =========================================================================
close all
load(output_file)
    TF_size       = 24;               %Title font size
    AF_size       = 20;               %Axes label font size
    AT_size       = 20;               %Tick label font size

set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', 40)
set(0,'DefaultTextFontName', 'Times New Roman')
set(0,'DefaultTextFontSize', 40)
set(0,'defaultLineLineWidth', 2)
    
time = time/3600;   %Simulation time from seconds to hours, for plots

%% Plots for EKF ==========================================================

%Figure: 3D Position Plot for GPS
%
figure(1)
    plot3(GPS_output(:,1)*10^3,GPS_output(:,2)*10^3,GPS_output(:,3)*10^3);
    %title('GPS Relative Position Measurements')
    xlabel('x Position (m)')
    ylabel('y Position (m)')
    zlabel('z Position (m)')
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    %axis([-1200*max(GPS_output(:,1)) 1200*max(GPS_output(:,1)) -1200*max(GPS_output(:,2)) ...
    %    1200*max(GPS_output(:,2)) -1200*max(GPS_output(:,1)) 1200*max(GPS_output(:,1))]);
    axis([ -800 800 -800 800 -100 100])
    view(0,90)
    %grid on
    hold on
    Target   = scatter3(0, 0, 0, 'filled', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'y', 'MarkerEdgeColor', 'k', 'SizeData', 10^2);
    legend('Chaser Trajectory', 'Target Spacecraft', 'Location','NorthEast')
    %file_print = fullfile(pathname, 'Fig_NoisyGPS');
    %print(file_print,  '-deps');   
    
%Figure: 3D Position Plot NERMs
%{
figure(2)
    plot3(NERM_output(:,1)*10^3, NERM_output(:,2)*10^3, NERM_output(:,3)*10^3);
    title('NERM Dynamics Relative Position Estimates')
    xlabel('X Position (m)')
    ylabel('Y Position (m)')
    zlabel('Z Position (m)')
    %axis([-200 200 -200 200])
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    axis([-2000*max(GPS_output(:,1)) 2000*max(GPS_output(:,1)) -4000*max(GPS_output(:,1)) ...
        4000*max(GPS_output(:,1)) -11000*max(GPS_output(:,3)) 11000*max(GPS_output(:,3))]);
    view(-35,50)
    grid on
    hold on
    Target   = scatter3(0, 0, 0, 'filled', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'y', 'MarkerEdgeColor', 'k', 'SizeData', 5^2);
    legend('Chaser Trajectory', 'Target Spacecraft', 'Location','NorthEast')
    if strcmp(print_plots, 'on')
        print -deps NERM_Prop
    end
%}

%Figure: 3D Position Plot for EKF
figure(3)
    plot3(EKF_output(:,1)*10^3, EKF_output(:,2)*10^3, EKF_output(:,3)*10^3);
    title('EKF Relative Position Estimates')
    xlabel('X Position (m)')
    ylabel('Y Position (m)')
    zlabel('Z Position (m)')
    %axis([-200 200 -200 200]);
    grid on
    axis([-1200*max(GPS_output(:,1)) 1200*max(GPS_output(:,1)) -1200*max(GPS_output(:,2)) ...
        1200*max(GPS_output(:,2)) -11000*max(GPS_output(:,3)) 11000*max(GPS_output(:,3))]);
    view(-35,50)
    grid on
    hold on
    Target   = scatter3(0, 0, 0, 'filled', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'y', 'MarkerEdgeColor', 'k', 'SizeData', 5^2);
    legend('Chaser Trajectory', 'Target Spacecraft', 'Location','NorthEast')
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    if strcmp(print_plots, 'on')
        print -deps EKF_Output
    end
%}
%==========================================================================  
% Figures or Figure w/ Subfigures: EKF State Errors
%{
% Figure: EKF X-Position Error
figure(4)
    %subplot(3,2,1)
    plot(time, errors_x*1000) 
    hold on
    plot(time, 3*sig_x*1000) 
    plot(time, -3*sig_x*1000) 
    legend('EKF Error', '3-\sigma Bounds')
    title('EKF-X Error')
    xlabel('Time (hrs)')
    ylabel('X Position (m)')
    legend('EKF Error', '3-\sigma Bounds')
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    grid on
    %axis([0 time(length(time)) -4 4])

%==========================================================================   
% Figure: EKF Y-Position Error
%subplot(3,2,3)
figure(5)
    plot(time, errors_y*1000) 
    hold on
    plot(time, 3*sig_y*1000) 
    plot(time, -3*sig_y*1000) 
    set(gca,'FontSize',AT_size);
    title('EKF-Y Error')
    xlabel('Time (hrs)')
    ylabel('Y Position (m)')
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    grid on
    %axis([0 time(length(time)) -4 4])
    if strcmp(print_plots, 'on')
        print -deps EKF_Errors
    end

%==========================================================================  
% Figure: EKF Z-Position Error
figure(6)
%subplot(3,2,5)
    plot(time, errors_z*1000) 
    hold on
    plot(time, 3*sig_z*1000) 
    plot(time, -3*sig_z*1000) 
    legend('EKF Error', '3-\sigma Bounds')
    set(gca,'FontSize',AT_size);
    title('EKF-Z Error')
    xlabel('Time (hrs)')
    ylabel('Z Position (m)')
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    grid on
    %axis([0 time(length(time)) -1 1])

%==========================================================================  
% Figure: EKF Theta Error
figure(7)    
%subplot(3,2,6)
    plot(time, errors_theta) 
    hold on
    plot(time, 3*sig_theta) 
    plot(time, -3*sig_theta) 
    legend('EKF Error', '3-\sigma Bounds')
    set(gca,'FontSize',AT_size);
    title('EKF \theta Error')
    ylabel('\theta (rad)')
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    grid on

%==========================================================================  
% Figure: EKF rt Error
figure(8)    
%subplot(3,2,6)
    plot(time, errors_rt) 
    hold on
    plot(time, 3*sig_rt)
    plot(time, -3*sig_rt)
    legend('EKF Error', '3-\sigma Bounds')
    set(gca,'FontSize',AT_size);
    title('EKF r_t Error')
    ylabel('r_t (km)')
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    grid on
    
%==========================================================================  
% Figure: EKF X-Velocity Error
figure(9)    
%subplot(3,2,2)
    plot(time, errors_xdot*1000) 
    hold on
    plot(time, 3*sig_xdot*1000) 
    plot(time, -3*sig_xdot*1000) 
    legend('EKF Error', '3-\sigma Bounds')
    set(gca,'FontSize',AT_size);
    title('EKF X-Velocity Error')
    ylabel('X Velocity (m/s)')
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    grid on
    %axis([0 time(length(time)) -.2 .2])

%==========================================================================   
% Figure: EKF Y-Velocity Error
figure(10)    
%subplot(3,2,4)
    plot(time, errors_ydot*1000) 
    hold on
    plot(time, 3*sig_ydot*1000) 
    plot(time, -3*sig_ydot*1000) 
    legend('EKF Error', '3-\sigma Bounds')
    set(gca,'FontSize',AT_size);
    title('EKF Y-Velocity Error')
    ylabel('Y Velocity (m/s)')
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    grid on
    %axis([0 time(length(time)) -.2 .2])

%==========================================================================  
% Figure: EKF Z-Velocity Error
figure(11)    
%subplot(3,2,6)
    plot(time, errors_zdot*1000) 
    hold on
    plot(time, 3*sig_zdot*1000) 
    plot(time, -3*sig_zdot*1000) 
    legend('EKF Error', '3-\sigma Bounds')
    set(gca,'FontSize',AT_size);
    title('EKF Z-Velocity Error')
    ylabel('Z Velocity (m/s)')
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    grid on
    %axis([0 time(length(time)) -0.04 0.04])

%==========================================================================  
% Figure: EKF Theta_dot Error
figure(12)    
%subplot(3,2,6)
    plot(time, errors_theta_dot) 
    hold on
    plot(time, 3*sig_theta_dot) 
    plot(time, -3*sig_theta_dot) 
    legend('EKF Error', '3-\sigma Bounds')
    set(gca,'FontSize',AT_size);
    title('EKF \theta dot Error')
    ylabel('\theta dot(rad/s)')
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    grid on

%==========================================================================  
% Figure: EKF rt_dot Error
figure(13)    
%subplot(3,2,6)
    plot(time, errors_rt_dot) 
    hold on
    plot(time, 3*sig_rt_dot)%*10^-4) 
    plot(time, -3*sig_rt_dot)%*10^-4) 
    legend('EKF Error', '3-\sigma Bounds')
    set(gca,'FontSize',AT_size);
    title('EKF r_t dot Error')
    ylabel('r_t dot (km/s)')
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    grid on
%}
%==========================================================================  
% Figures: EKF State Outputs
%{  
    figure(14)
    clf
    plot(time, 1000*EKF_output(:,1), '-', time, 1000*NERM_output(:,1), '-.', time, 1000*GPS_output(:,1), '--')
    legend('EKF Output', 'Dynamics', 'GPS Measurements')
    set(gca,'FontSize',AT_size);
    title('EKF X Position vs Time')
    ylabel('X Position (m)')
    xlabel('Time (hr)')
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    grid on  
    
    figure(15)
    clf
    plot(time, 1000*EKF_output(:,2), '-', time, 1000*NERM_output(:,2), '-.', time, 1000*GPS_output(:,2), '--')
    legend('EKF Output', 'Dynamics', 'GPS Measurements')
    set(gca,'FontSize',AT_size);
    title('EKF Output Y vs Time')
    xlabel('Time (hr)')
    ylabel('Y Position (m)')
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    grid on  
    
    figure(16)
    clf
    plot(time, 1000*EKF_output(:,3), '-', time, 1000*NERM_output(:,3), '-.', time, 1000*GPS_output(:,3), '--')
    legend('EKF Output', 'Dynamics', 'GPS Measurements')
    set(gca,'FontSize',AT_size);
    title('EKF Output Z vs Time')
    xlabel('Time (hrs)')
    ylabel('Z Position (m)')
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    %axis([0 time(length(time)) -5 5]);
    grid on  

    figure(17)
    clf
    plot(time, EKF_output(:,4), '-', time, NERM_output(:,4),'-.', time, theta_target, '--')
    set(gca,'FontSize',AT_size);
    title('EKF Output \theta vs Time')
    xlabel('Time (hrs)')
    ylabel('\theta (rad)')
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    %axis([0 time(length(time)) -1 20]);
    legend('EKF Output', 'Dynamics', 'Virtual Measurement') 
    grid on  
    
    figure(18)
    clf
    plot(time, EKF_output(:,5), '-', time, NERM_output(:,5), '-.', time, rt_clean, '--')
    set(gca,'FontSize',AT_size);
    title('EKF Output R_t vs Time')
    xlabel('Time (hrs)')
    ylabel('R_t (km)')
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    %axis([0 time(length(time)) 0 9000]);
    legend('EKF Output', 'Dynamics', 'Virtual Measurement')
    grid on  
    
    figure(19)
    clf
    plot(time, EKF_output(:,6)*10^3, '-', time, NERM_output(:,6)*10^3, '-.', time, GPS_output(:,5)*10^3, '--')
    set(gca,'FontSize',AT_size);
    title('EKF Output X-dot vs Time')
    xlabel('Time (hrs)')
    ylabel('X_d_o_t (m/s)')
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    legend('EKF Output', 'Dynamics', 'GPS Measurement')
    grid on  
    
    figure(20)
    clf
    plot(time, EKF_output(:,7)*10^3, '-', time, 1000*NERM_output(:,7), '-.', time, GPS_output(:,6)*10^3, '--')
    set(gca,'FontSize',AT_size);
    title('EKF Output Y-dot vs Time')
    xlabel('Time (hrs)')
    ylabel('Y_d_o_t (m/s)')
    legend('EKF Output', 'Dynamics', 'GPS Measurement')
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    grid on  
    
    figure(21)
    clf
    plot(time, EKF_output(:,8)*10^3, '-', time, 1000*NERM_output(:,8), '-.', time, GPS_output(:,7)*10^3, '--')
    hold on
    set(gca,'FontSize',AT_size);
    title('EKF Output Z-dot vs Time')
    xlabel('Time (hrs)')
    ylabel('Z_d_o_t (m/s)')
    legend('EKF Output', 'Dynamics', 'GPS Measurement')
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    grid on  
        
    figure(22)
    clf
    plot(time, EKF_output(:,9), time, NERM_output(:,9), '-.')
    set(gca,'FontSize',AT_size);
    title('EKF Output \theta-dot vs Time')
    xlabel('Time (hrs)')
    ylabel('\theta_dot (rad/s)')
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    legend('EKF Output', 'Dynamics')
    grid on  
    
    figure(23)
    clf
    plot(time, EKF_output(:,10), time, NERM_output(:,10), '-.')
    set(gca,'FontSize',AT_size);
    title('EKF Output R_t-dot vs Time')
    xlabel('Time (hrs)')
    ylabel('R_t-dot (km/s)')
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    %axis([0 time(length(time)) -1 1])
    legend('EKF Output', 'Dynamics')
    grid on  
 %}   
%==========================================================================  
%Figure w/ Subfigures: EKF State Outputs
%{
figure(15)
%title('EKF Output States')    
    h(1) = subplot(5,2,1);
    plot(time, 10^3*EKF_output(:,1)) 
    set(gca,'FontSize',AT_size);
    title('EKF Output States')
    ylabel('X(m)')
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    set(gca,'XTick',[]);
    %set(gca,'YTick',[]);
    grid on
    %axis([0 time(length(time)) -100 100])
    
    h(2) = subplot(5,2,3);
    plot(time, 10^3*EKF_output(:,2)) 
    set(gca,'FontSize',AT_size);
    ylabel('Y(m)')
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    set(gca,'XTick',[]);
    %set(gca,'YTick',[]);
    grid on
    %axis([0 time(length(time)) -200 200])
    
    h(3) = subplot(5,2,5);
    plot(time, 10^3*EKF_output(:,3)) 
    set(gca,'FontSize',AT_size);
    ylabel('Z (m)')
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    set(gca,'XTick',[]);
    %set(gca,'YTick',[]);
    grid on
    %axis([0 time(length(time)) -0.1 0.1])
    
    h(4) = subplot(5,2,7);
    plot(time, EKF_output(:,4)) 
    set(gca,'FontSize',AT_size);
    ylabel('\theta (rad)')
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    set(gca,'XTick',[]);
    %set(gca,'YTick',[]);
    grid on
    %axis([0 time(length(time)) 0 300])    
    
    h(5) = subplot(5,2,9);
    plot(time, EKF_output(:,5)) 
    set(gca,'FontSize',AT_size);
    xlabel('Time (hrs)')
    ylabel('r_t (km)')
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    grid on
    %axis([0 time(length(time)) 7499 7501])
    
    h(6) = subplot(5,2,2);
    plot(time, 10^3*EKF_output(:,6)) 
    set(gca,'FontSize',AT_size);
    ylabel('X_d_o_t (m/s)')
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    set(gca,'XTick',[]);
    %set(gca,'YTick',[]);
    grid on
    %axis([0 time(length(time)) -0.1 0.1])
    
    h(7) = subplot(5,2,4);
    plot(time, 10^3*EKF_output(:,7)) 
    set(gca,'FontSize',AT_size);
    ylabel('Y_d_o_t (m/s)')
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    set(gca,'XTick',[]);
    %set(gca,'YTick',[]);
    grid on
    %axis([0 time(length(time)) -0.2 0.2])
    
    h(8) = subplot(5,2,6); 
    plot(time, 10^3*EKF_output(:,8)) 
    set(gca,'FontSize',AT_size);
    ylabel('Z_d_o_t (m/s)')
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    set(gca,'XTick',[]);
    %set(gca,'YTick',[]);
    grid on
    %axis([0 time(length(time)) -0.01 0.01])
    
    h(9) = subplot(5,2,8);
    plot(time, EKF_output(:,9)) 
    set(gca,'FontSize',AT_size);
    ylabel('\theta_d_o_t (rad/s)')
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    set(gca,'XTick',[]);
    %set(gca,'YTick',[]);
    grid on
    %axis([0 time(length(time)) 0 2*10^-3])
    
    h(10) = subplot(5,2,10);
    plot(time, EKF_output(:,10)) 
    set(gca,'FontSize',AT_size);
    ylabel('r_t_,_d_o_t')
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    %set(gca,'XTick',[]);
    %set(gca,'YTick',[]);
    grid on
    %axis([0 time(length(time)) -2*10^-5 2*10^-5])
    xlabel('Time (hrs)')

    for j = 1:10
        p(j,:) = get(h(j), 'pos');
        p(j,4) = p(j,4) + 0.015;
        set(h(j), 'pos', p(j,:));
        set(get(h(j),'YLabel'),'Rotation',0)
    end   
    
  %}
%==========================================================================  
% Figure: EKF State Covariances
%{
    figure(16)
    clf
    num = 11;
    plot(time, EKF_output(:, num), time, EKF_output(:, num+11), time, EKF_output(:, num+22),...
        time, EKF_output(:, num+33), time, EKF_output(:, num+44), time, EKF_output(:, num+55),...
        time, EKF_output(:, num+66), time, EKF_output(:, num+77),time, EKF_output(:, num+88), time, EKF_output(:, num+99))
    set(gca,'FontSize',AT_size);
    title('State Covariances vs Time')
    xlabel('Time (hr)')
    ylabel('State Covariance')
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    legend('x', 'y', 'z', '\theta', 'rt', 'x_d_o_t', 'y_d_o_t', 'z_d_o_t', '\theta_d_o_t', 'rt_d_o_t')
    grid on  
 
    figure(17)
    plot(time, EKF_output(:, 55))
    title('r_t Covariance vs Time')
    xlabel('Time (hr)')
    ylabel('r_t Covariance')
    %set(gca,'FontName', 'Times')
    grid on  
 %} 
%==========================================================================  
% Figure: Kalman Gains
%{
    figure(16)
    clf
    num = 119; %Index for first column of output data cxontaining gain data
    plot(time, EKF_output(:, num), time, EKF_output(:, num+1), time, EKF_output(:, num+2),...
        time, EKF_output(:, num+3), time, EKF_output(:, num+4), time, EKF_output(:, num+5),...
        time, EKF_output(:, num+6), time, EKF_output(:, num+7))
    set(gca,'FontSize',AT_size);
    title('X-Kalman Gains vs Time')
    xlabel('Time (hr)')
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
    xlabel('Time (hr)')
    ylabel('Kalman Gain')
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    legend('x', 'y', 'z', '\theta', 'rt', 'x_d_o_t', 'y_d_o_t', 'z_d_o_t')
    grid on  
    
    figure(18)
    clf
    num = 135;
    plot(time, EKF_output(:, num), time, EKF_output(:, num+1), time, EKF_output(:, num+2),...
        time, EKF_output(:, num+3), time, EKF_output(:, num+4), time, EKF_output(:, num+5),...
        time, EKF_output(:, num+6), time, EKF_output(:, num+7))
    set(gca,'FontSize',AT_size);
    title('Z-Kalman Gains vs Time')
    xlabel('Time (hr)')
    ylabel('Kalman Gain')
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    legend('x', 'y', 'z', '\theta', 'rt', 'x_d_o_t', 'y_d_o_t', 'z_d_o_t')
    grid on  
    
    figure(19)
    clf
    num = 143;
    plot(time, EKF_output(:, num), time, EKF_output(:, num+1), time, EKF_output(:, num+2),...
        time, EKF_output(:, num+3), time, EKF_output(:, num+4), time, EKF_output(:, num+5),...
        time, EKF_output(:, num+6), time, EKF_output(:, num+7))
    set(gca,'FontSize',AT_size);
    title('\theta-Kalman Gains vs Time')
    xlabel('Time (hr)')
    ylabel('Kalman Gain')
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    legend('x', 'y', 'z', '\theta', 'rt', 'x_d_o_t', 'y_d_o_t', 'z_d_o_t')
    grid on  
    
    figure(20)
    clf
    num = 151;
    plot(time, EKF_output(:, num), time, EKF_output(:, num+1), time, EKF_output(:, num+2),...
        time, EKF_output(:, num+3), time, EKF_output(:, num+4), time, EKF_output(:, num+5),...
        time, EKF_output(:, num+6), time, EKF_output(:, num+7))
    set(gca,'FontSize',AT_size);
    title('r_t-Kalman Gains vs Time')
    xlabel('Time (hrs)')
    ylabel('Kalman Gain')
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    legend('x', 'y', 'z', '\theta', 'rt', 'x_d_o_t', 'y_d_o_t', 'z_d_o_t')
    grid on  
    
    figure(21)
    clf
    num = 159;
    plot(time, EKF_output(:, num), time, EKF_output(:, num+1), time, EKF_output(:, num+2),...
        time, EKF_output(:, num+3), time, EKF_output(:, num+4), time, EKF_output(:, num+5),...
        time, EKF_output(:, num+6), time, EKF_output(:, num+7))
    set(gca,'FontSize',AT_size);
    title('x_d_o_t-Kalman Gains vs Time')
    xlabel('Time (hrs)')
    ylabel('Kalman Gain')
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    legend('x', 'y', 'z', '\theta', 'rt', 'x_d_o_t', 'y_d_o_t', 'z_d_o_t')
    grid on  
    
    figure(22)
    clf
    num = 167;
    plot(time, EKF_output(:, num), time, EKF_output(:, num+1), time, EKF_output(:, num+2),...
        time, EKF_output(:, num+3), time, EKF_output(:, num+4), time, EKF_output(:, num+5),...
        time, EKF_output(:, num+6), time, EKF_output(:, num+7))
    set(gca,'FontSize',AT_size);
    title('y_d_o_t-Kalman Gains vs Time')
    xlabel('Time (hrs)')
    ylabel('Kalman Gain')
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    legend('x', 'y', 'z', '\theta', 'rt', 'x_d_o_t', 'y_d_o_t', 'z_d_o_t')
    grid on  
    
    figure(23)
    clf
    num = 175;
    plot(time, EKF_output(:, num), time, EKF_output(:, num+1), time, EKF_output(:, num+2),...
        time, EKF_output(:, num+3), time, EKF_output(:, num+4), time, EKF_output(:, num+5),...
        time, EKF_output(:, num+6), time, EKF_output(:, num+7))
    set(gca,'FontSize',AT_size);
    title('z_d_o_t-Kalman Gains vs Time')
    xlabel('Time (hrs)')
    ylabel('Kalman Gain')
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    legend('x', 'y', 'z', '\theta', 'rt', 'x_d_o_t', 'y_d_o_t', 'z_d_o_t')
    grid on  
    
    figure(24)
    clf
    num = 183;
    plot(time, EKF_output(:, num), time, EKF_output(:, num+1), time, EKF_output(:, num+2),...
        time, EKF_output(:, num+3), time, EKF_output(:, num+4), time, EKF_output(:, num+5),...
        time, EKF_output(:, num+6), time, EKF_output(:, num+7))
    set(gca,'FontSize',AT_size);
    title('\theta_d_o_t-Kalman Gains vs Time')
    xlabel('Time (hrs)')
    ylabel('Kalman Gain')
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    legend('x', 'y', 'z', '\theta', 'rt', 'x_d_o_t', 'y_d_o_t', 'z_d_o_t')
    grid on  
    
    figure(25)
    clf
    num = 191;
    plot(time, EKF_output(:, num), time, EKF_output(:, num+1), time, EKF_output(:, num+2),...
        time, EKF_output(:, num+3), time, EKF_output(:, num+4), time, EKF_output(:, num+5),...
        time, EKF_output(:, num+6), time, EKF_output(:, num+7))
    set(gca,'FontSize',AT_size);
    title('r_t_d_o_t-Kalman Gains vs Time')
    xlabel('Time (hrs)')
    ylabel('Kalman Gain')
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    legend('x', 'y', 'z', '\theta', 'rt', 'x_d_o_t', 'y_d_o_t', 'z_d_o_t')
    grid on  
    
    figure(26)
    clf
    num = 199;
    plot(time, EKF_output(:, num), time, EKF_output(:, num+1), time, EKF_output(:, num+2),...
        time, EKF_output(:, num+3), time, EKF_output(:, num+4), time, EKF_output(:, num+5),...
        time, EKF_output(:, num+6), time, EKF_output(:, num+7))
    set(gca,'FontSize',AT_size);
    title('EKF Residuals vs Time')
    xlabel('Time (hrs)')
    ylabel('Residual')
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    legend('x', 'y', 'z', '\theta', 'rt', 'x_d_o_t', 'y_d_o_t', 'z_d_o_t')
    grid on  
    
    figure(27)
    clf
    num = 11;
    plot(time, EKF_output(:, num), time, EKF_output(:, num+11), time, EKF_output(:, num+22),...
        time, EKF_output(:, num+33), '--', time, EKF_output(:, num+44), time, EKF_output(:, num+55),...
        time, EKF_output(:, num+66), time, EKF_output(:, num+77), time, EKF_output(:, num+88), ...
        time, EKF_output(:, num+99))
    set(gca,'FontSize',AT_size);
    title('EKF a Posteriori State Covariance Diagonal Elements vs Time')
    xlabel('Time (hrs)')
    ylabel('Diagonal Covariance Element Value')
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    axis([0 time(length(time)) 0 0.1])
    legend('x', 'y', 'z', '\theta',' r_t', 'x_d_o_t', 'y_d_o_t', 'z_d_o_t', '\theta_d_o_t', 'rt_d_o_t')
    grid on  
    %}
%==========================================================================  
% Animation: Relative Positions
%{
figure(30)
title('Relative Positions', ... 
    'Fontweight', 'bold', 'FontSize', 12)
curve1 = animatedline('color','r', 'LineWidth', 0.5);
curve2 = animatedline('color','b', 'LineWidth', 0.5);
curve3 = animatedline('color','g', 'LineWidth', 0.5);
set(gca, 'XLim', [-200 200], 'YLim', ...
    [-200 200], 'ZLim', [-2*10^-7 2*10^-7]); 
xlabel(' x (m) ')
ylabel(' y (m) ')
zlabel(' z (m) ')
view(-40, 30);
set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
grid on
hold on;
EKF  = scatter3(1000*EKF_output(1,1), 1000*EKF_output(1,2), 1000*EKF_output(1,3), 'MarkerEdgeColor', 'r');
NERM  = scatter3(1000*NERM_output(1,1), 1000*NERM_output(1,2), 1000*NERM_output(1,3), 'MarkerEdgeColor', 'b');
GPS  = scatter3(1000*GPS_output(1,1), 1000*GPS_output(1,2), 1000*GPS_output(1,3), 'MarkerEdgeColor', 'b');
legend('EKF', 'Dynamics', 'GPS')
for i = 1:100:30000;
    addpoints(curve1, 1000*EKF_output(i,1), 1000*EKF_output(i,2), 1000*EKF_output(i,3));
    addpoints(curve2, 1000*NERM_output(i,1), 1000*NERM_output(i,2), 1000*NERM_output(i,3));
    addpoints(curve3, 1000*GPS_output(i,1), 1000*GPS_output(i,2), 1000*GPS_output(i,3));
    EKF  = scatter3(1000*EKF_output(i,1), 1000*EKF_output(i,2), 1000*EKF_output(i,3), 'MarkerEdgeColor', 'r');
    NERM  = scatter3(1000*NERM_output(i,1), 1000*NERM_output(i,2), 1000*NERM_output(i,3), 'MarkerEdgeColor', 'b');
    GPS  = scatter3(1000*GPS_output(i,1), 1000*GPS_output(i,2), 1000*GPS_output(i,3), 'MarkerEdgeColor', 'b');
    drawnow
    delete(EKF);
    delete(GPS);
    delete(NERM);
end
%}
%==========================================================================  
% Animation: Relative Velocities
%{
figure(31)
title('Relative Velocities Velocities', ... 
    'Fontweight', 'bold', 'FontSize', 12)
curve1 = animatedline('color','r', 'LineWidth', 0.5);
curve2 = animatedline('color','b', 'LineWidth', 0.5);
curve3 = animatedline('color','g', 'LineWidth', 0.5);
%set(gca, 'XLim', [-.2 .2], 'YLim', ...
%    [-.2 .2], 'ZLim', [-2*10^-7 2*10^-7]); 
xlabel(' x Velocity (m/s) ')
ylabel(' y Velocity (m/s) ')
zlabel(' z Velocity (m/s) ')
view(-40, 30);
set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
grid on
hold on;
EKF  = scatter3(1000*EKF_output(1,4), 1000*EKF_output(1,5), 1000*EKF_output(1,6), 'MarkerEdgeColor', 'r');
NERM  = scatter3(1000*NERM_output(1,4), 1000*NERM_output(1,5), 1000*NERM_output(1,6), 'MarkerEdgeColor', 'b');
GPS  = scatter3(1000*GPS_output(1,4), 1000*GPS_output(1,5), 1000*GPS_output(1,6), 'MarkerEdgeColor', 'b');
legend('EKF', 'Dynamics', 'Measurements')
for i = 1:100:30000;
    addpoints(curve1, 1000*EKF_output(i,4), 1000*EKF_output(i,5), 1000*EKF_output(i,6));
    addpoints(curve2, 1000*NERM_output(i,4), 1000*NERM_output(i,5), 1000*NERM_output(i,6));
    addpoints(curve3, 1000*GPS_output(i,4), 1000*GPS_output(i,5), 1000*GPS_output(i,6));
    EKF  = scatter3(1000*EKF_output(i,4), 1000*EKF_output(i,5), 1000*EKF_output(i,6), 'MarkerEdgeColor', 'r');
    NERM  = scatter3(1000*NERM_output(i,4), 1000*NERM_output(i,5), 1000*NERM_output(i,6), 'MarkerEdgeColor', 'b');
    GPS  = scatter3(1000*GPS_output(i,4), 1000*GPS_output(i,5), 1000*GPS_output(i,6), 'MarkerEdgeColor', 'b');
    drawnow
    delete(EKF);
    delete(GPS);
    delete(NERM);
end
  %}  

%==========================================================================  
% Bar Graphs: Relative Position and Velocity RMS
filter_names = {'GPS', 'NERM EKF'};
RMS_pos_values = [GPS_RMSE_Avg_pos EKF_RMSE_Avg_pos]; 
RMS_vel_values = [GPS_RMSE_Avg_vel EKF_RMSE_Avg_vel]; 

figure(40)
bar(RMS_pos_values*1000)
set(gca,'xticklabel',filter_names);
ylabel('Position Error RMS (m)')
%set(gcf, 'units','normalized','outerposition',[0 0 1 1]);

figure(41)
bar(RMS_vel_values*1000)
set(gca,'xticklabel',filter_names);
ylabel('Velocity Error RMS (m/s)')
%set(gcf, 'units','normalized','outerposition',[0 0 1 1]);

%% Plots for Adaptive EKF =================================================  
plot_start  = 1;
plot_end    = 300;%length(time);

for i = 1:plot_end
    trace_Q(i) = trace(Q_output(:,:,i));
    trace_R(i) = trace(R_output(:,:,i));
end

% Figure: Trace of Q Covariance 
figure(50)
    plot(time(plot_start:plot_end), trace_Q(plot_start:plot_end))
    hold on
    plot([0 time(plot_end)],[trace(Q0) trace(Q0)],'--')
    title('Trace of Q')
    xlabel('Time (hrs)')
    ylabel('Tr(Q)')
    legend('Adapted Q', 'Initial Q0')
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    grid on
    
% Figure: Trace of R Covariance 
figure(51)
    plot(time(plot_start:plot_end), trace_R(plot_start:plot_end))
    hold on
    plot([0 time(plot_end)],[trace(R0) trace(R0)],'--')
    title('Trace of R')
    xlabel('Time (hrs)')
    ylabel('Tr(R)')
    legend('Adapted R', 'Initial R0')
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    grid on

% Figure: Trace of P Covariance (Posteriori)
figure(52)
    plot(time, trace_P)
    title('Trace of P^+')
    xlabel('Time (hrs)')
    ylabel('Tr(P)')
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    grid on
    
% Figure: Likelihood Function
figure(53)
    plot(time(plot_start:plot_end), L_sum(plot_start:plot_end))
    title('Log-Likelihood Function')
    xlabel('Time (hrs)')
    ylabel('-log L(\theta)')
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    grid on    

% =========================================================================
% Figure w/ Subfigures: Q Covariance Plots
figure(60)
set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    subplot(5,2,1)
    plot(time(plot_start:plot_end), squeeze(Q_output(1,1,plot_start:plot_end))) 
    set(gca,'xticklabel',[])
    ylabel('q_{11}')
    grid on
 
    subplot(5,2,2)
    plot(time(plot_start:plot_end), squeeze(Q_output(6,6,plot_start:plot_end))) 
    ylabel('q_{66}')
    set(gca,'xticklabel',[])
    grid on
    
    subplot(5,2,3)
    plot(time(plot_start:plot_end), squeeze(Q_output(2,2,plot_start:plot_end))) 
    ylabel('q_{22}')
    set(gca,'xticklabel',[])
    grid on
 
    subplot(5,2,4)
    plot(time(plot_start:plot_end), squeeze(Q_output(7,7,plot_start:plot_end))) 
    ylabel('q_{77}')
    set(gca,'xticklabel',[])
    grid on
 
    subplot(5,2,5)
    plot(time(plot_start:plot_end), squeeze(Q_output(3,3,plot_start:plot_end))) 
    ylabel('q_{33}')
    set(gca,'xticklabel',[])
    grid on
 
    subplot(5,2,6)
    plot(time(plot_start:plot_end), squeeze(Q_output(8,8,plot_start:plot_end))) 
    ylabel('q_{88}')
    set(gca,'xticklabel',[])
    grid on
 
    subplot(5,2,7)
    plot(time(plot_start:plot_end), squeeze(Q_output(4,4,plot_start:plot_end))) 
    ylabel('q_{44}')
    set(gca,'xticklabel',[])
    grid on
 
    subplot(5,2,8)
    plot(time(plot_start:plot_end), squeeze(Q_output(9,9,plot_start:plot_end))) 
    ylabel('q_{99}')
    set(gca,'xticklabel',[])
    grid on
 
    subplot(5,2,9)
    plot(time(plot_start:plot_end), squeeze(Q_output(5,5,plot_start:plot_end))) 
    %xlabel('Time (hrs)')
    ylabel('q_{55}')
    grid on
 
    subplot(5,2,10)
    plot(time(plot_start:plot_end), squeeze(Q_output(10,10,plot_start:plot_end))) 
    %xlabel('Time (hrs)')
    ylabel('q_{10,10}')
    grid on
 
% =========================================================================
% Figure w/ Subfigures: R Covariance Plots
figure(61)
set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    subplot(4,2,1)
    plot(time(plot_start:plot_end), squeeze(R_output(1,1,plot_start:plot_end))) 
    set(gca,'xticklabel',[])
    ylabel('R_{11}')
    grid on
 
    subplot(4,2,2)
    plot(time(plot_start:plot_end), squeeze(R_output(5,5,plot_start:plot_end))) 
    ylabel('R_{55}')
    set(gca,'xticklabel',[])
    grid on
    
    subplot(4,2,3)
    plot(time(plot_start:plot_end), squeeze(R_output(2,2,plot_start:plot_end))) 
    ylabel('R_{22}')
    set(gca,'xticklabel',[])
    grid on
 
    subplot(4,2,4)
    plot(time(plot_start:plot_end), squeeze(R_output(6,6,plot_start:plot_end))) 
    ylabel('R_{66}')
    set(gca,'xticklabel',[])
    grid on
 
    subplot(4,2,5)
    plot(time(plot_start:plot_end), squeeze(R_output(3,3,plot_start:plot_end))) 
    ylabel('R_{33}')
    set(gca,'xticklabel',[])
    grid on
 
    subplot(4,2,6)
    plot(time(plot_start:plot_end), squeeze(R_output(7,7,plot_start:plot_end))) 
    ylabel('R_{77}')
    set(gca,'xticklabel',[])
    grid on
 
    subplot(4,2,7)
    plot(time(plot_start:plot_end), squeeze(R_output(4,4,plot_start:plot_end))) 
    ylabel('R_{44}')
    set(gca,'xticklabel',[])
    grid on

  
    
% =========================================================================
% Figure w/ Subfigures: Observed/Smoothed Residual Covariances
%{
figure(70)
    subplot(4,2,1)
    plot(time(plot_start:plot_end), Pv_obs(plot_start:plot_end, 1), ...
        '--y', 'LineWidth', 0.05) 
    title('Observed and Smoothed Residual Covariances')
    hold on
    plot(time(plot_start:plot_end), Pv_smoothed(plot_start:plot_end, 1), ...
        '-k', 'LineWidth', 0.5) 
    ylabel('x')
    set(gca,'xticklabel',[])
    grid on
 
    subplot(4,2,2)
    plot(time(plot_start:plot_end), Pv_obs(plot_start:plot_end, 5),...
        '--y', 'LineWidth', 0.05)
    hold on
    plot(time(plot_start:plot_end), Pv_smoothed(plot_start:plot_end, 5),...
        '-k', 'LineWidth', 0.5)  
    ylabel('x_d_o_t')
    set(gca,'xticklabel',[])
    grid on
    
    subplot(4,2,3)
    plot(time(plot_start:plot_end), Pv_obs(plot_start:plot_end, 2),...
        '--y', 'LineWidth', 0.05)
    hold on
    set(gca,'xticklabel',[])
    plot(time(plot_start:plot_end), Pv_smoothed(plot_start:plot_end, 2),...
        '-k', 'LineWidth', 0.5)  
    ylabel('y')
    grid on
 
    subplot(4,2,4)
    plot(time(plot_start:plot_end), Pv_obs(plot_start:plot_end, 2),...
        '--y', 'LineWidth', 0.05)
    hold on
    plot(time(plot_start:plot_end), Pv_smoothed(plot_start:plot_end, 6), ...
        '-k', 'LineWidth', 0.5)  
    set(gca,'xticklabel',[])
    ylabel('y_d_o_t')
    grid on
 
    subplot(4,2,5)
    plot(time(plot_start:plot_end), Pv_obs(plot_start:plot_end, 2),...
        '--y', 'LineWidth', 0.05) 
    hold on
    plot(time(plot_start:plot_end), Pv_smoothed(plot_start:plot_end, 3), ...
        '-k', 'LineWidth', 0.5)  
    set(gca,'xticklabel',[])
    ylabel('z')
    grid on
 
    subplot(4,2,6)
    plot(time(plot_start:plot_end), Pv_obs(plot_start:plot_end, 2),...
        '--y', 'LineWidth', 0.05)
    hold on
    plot(time(plot_start:plot_end), Pv_smoothed(plot_start:plot_end, 7), ...
        '-k', 'LineWidth', 0.5)      
    ylabel('z_d_o_t')
    set(gca,'xticklabel',[])
    grid on

    subplot(4,2,7)
    plot(time(plot_start:plot_end), Pv_obs(plot_start:plot_end, 2),...
        '--y', 'LineWidth', 0.05)
    hold on
    plot(time(plot_start:plot_end), Pv_smoothed(plot_start:plot_end, 4), ...
        '-k', 'LineWidth', 0.5)    
    xlabel('Time (hrs)')
    ylabel('\theta')
    legend('Smoothed', 'Raw')
    grid on

% Figure w/ Subfigures: Theoretical Residual Covariances    
figure(71)
    subplot(4,2,1)
    plot(time(plot_start:plot_end), Pr_theo(plot_start:plot_end, 1)) 
    title('Theoretical Residual Covariances')
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
    xlabel('Time (hrs)')
    ylabel('\theta')
    grid on
  
% Figure w/ Subfigures: Residual Covariance Error Signals (Theoretical - Observed)
figure(72)
    subplot(4,2,1)
    plot(time(plot_start:plot_end), P_error(plot_start:plot_end, 1)) 
    title('Errors from Residual Covariances')
    hold on
    ylabel('x')
    set(gca,'xticklabel',[])
    grid on
 
    subplot(4,2,2)
    plot(time(plot_start:plot_end), P_error(plot_start:plot_end, 5)) 
    hold on
    ylabel('x_d_o_t')
    set(gca,'xticklabel',[])
    grid on
    
    subplot(4,2,3)
    plot(time(plot_start:plot_end), P_error(plot_start:plot_end, 2)) 
    hold on
    set(gca,'xticklabel',[])
    ylabel('y')
    grid on
 
    subplot(4,2,4)
    plot(time(plot_start:plot_end), P_error(plot_start:plot_end, 6)) 
    hold on
    set(gca,'xticklabel',[])
    ylabel('y_d_o_t')
    grid on
 
    subplot(4,2,5)
    plot(time(plot_start:plot_end), P_error(plot_start:plot_end, 3)) 
    hold on
    set(gca,'xticklabel',[])
    ylabel('z')
    grid on
 
    subplot(4,2,6)
    plot(time(plot_start:plot_end), P_error(plot_start:plot_end, 7)) 
    hold on
    ylabel('z_d_o_t')
    set(gca,'xticklabel',[])
    grid on

    subplot(4,2,7)
    plot(time(plot_start:plot_end), P_error(plot_start:plot_end, 4)) 
    hold on
    xlabel('Time (hrs)')
    ylabel('\theta')
    grid on

% Figure: Residual Covariance Error Derivative Signal
%{
figure(73)
    subplot(4,2,1)
    plot(time(plot_start:plot_end), Pdot_error(plot_start:plot_end, 1)) 
    title('Derivative of Errors from Residual Covariances')
    hold on
    ylabel('x')
    set(gca,'xticklabel',[])
    grid on
 
    subplot(4,2,2)
    plot(time(plot_start:plot_end), Pdot_error(plot_start:plot_end, 5)) 
    hold on
    ylabel('x_d_o_t')
    set(gca,'xticklabel',[])
    grid on
    
    subplot(4,2,3)
    plot(time(plot_start:plot_end), Pdot_error(plot_start:plot_end, 2)) 
    hold on
    set(gca,'xticklabel',[])
    ylabel('y')
    grid on
 
    subplot(4,2,4)
    plot(time(plot_start:plot_end), Pdot_error(plot_start:plot_end, 6)) 
    hold on
    set(gca,'xticklabel',[])
    ylabel('y_d_o_t')
    grid on
 
    subplot(4,2,5)
    plot(time(plot_start:plot_end), Pdot_error(plot_start:plot_end, 3)) 
    hold on
    set(gca,'xticklabel',[])
    ylabel('z')
    grid on
 
    subplot(4,2,6)
    plot(time(plot_start:plot_end), Pdot_error(plot_start:plot_end, 7)) 
    hold on
    ylabel('z_d_o_t')
    set(gca,'xticklabel',[])
    grid on

    subplot(4,2,7)
    plot(time(plot_start:plot_end), Pdot_error(plot_start:plot_end, 4)) 
    hold on
    xlabel('Time (hrs)')
    ylabel('\theta')
    grid on    
%}
% =========================================================================
%Plots for the DMSAC Gains
%{
plot_start  = 1;
plot_end    = 1000;

for i = 1:plot_end
    trace_Kpe(i) = trace(K_pe(:,:,i));
    trace_KIe(i) = trace(K_Ie(:,:,i));
    trace_Ke(i) = trace(K_e(:,:,i));
end

% Figure: Trace of K_Ie
figure(80)
    plot(time(plot_start:plot_end), trace_KIe(plot_start:plot_end))
    hold on
    title('Trace of Integral Gain Matrix K_I_e')
    xlabel('Time (hrs)')
    ylabel('Tr(K_I_e)')
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    grid on 

% Figure: Trace of K_pe
figure(81)
    plot(time(plot_start:plot_end), trace_Kpe(plot_start:plot_end))
    hold on
    title('Trace of Proportional Gain Matrix K_p_e')
    xlabel('Time (hrs)')
    ylabel('Tr(K_p_e)')
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    grid on

% Figure: Trace of K_e
figure(82)
    plot(time(plot_start:plot_end), trace_Ke(plot_start:plot_end))
    hold on
    title('Trace of Gain Matrix K_e')
    xlabel('Time (hrs)')
    ylabel('Tr(K_e)')
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    grid on
%}    
% =========================================================================
%}