function [] = FFNAV_HCW_EKF_Plotter(output_file)
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
% Last Edited: Cory Fraser - July 11, 2017
%% ========================================================================
%close all
load(output_file)
    TF_size       = 24;               %Title font size
    AF_size       = 20;               %Axes label font size
    AT_size       = 20;               %Tick label font size

set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', 20)
set(0,'DefaultTextFontName', 'Times New Roman')
set(0,'DefaultTextFontSize', 24)
%set(0,'defaultLineLineWidth', 2)
    
time = time/3600;   %converting simulation time in seconds to hours for plots

%{
%Noisy Relative Position GPS Measurements
figure(1)
    plot3(GPS_output(:,1)*10^3,GPS_output(:,2)*10^3,GPS_output(:,3)*10^3);
    title('Noisy Relative Position GPS Measurements')
    xlabel('X Position (m)')
    ylabel('Y Position (m)')
    zlabel('Z Position (m)')
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    axis([-2000*max(GPS_output(:,1)) 2000*max(GPS_output(:,1)) -2000*max(GPS_output(:,1)) ...
        2000*max(GPS_output(:,1)) -11000*max(GPS_output(:,3)) 11000*max(GPS_output(:,3))]);
    view(-35,50)
    grid on
    hold on
    Target   = scatter3(0, 0, 0, 'filled', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'y', 'MarkerEdgeColor', 'k', 'SizeData', 5^2);
    legend('Chaser Trajectory', 'Target Spacecraft', 'Location','NorthEast')
    if strcmp(print_plots, 'on')
        print -deps NoisyGPS
    end

%Propagated Nonlinear Equations of Motion

figure(2)
    plot3(HCW_output(:,1)*10^3, HCW_output(:,2)*10^3, HCW_output(:,3)*10^3);
    title('Dynamics Relative Position Estimates')
    xlabel('X Position (m)')
    ylabel('Y Position (m)')
    zlabel('Z Position (m)')
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    axis([-2000*max(GPS_output(:,1)) 2000*max(GPS_output(:,1)) -2000*max(GPS_output(:,1)) ...
        2000*max(GPS_output(:,1)) -11000*max(GPS_output(:,3)) 11000*max(GPS_output(:,3))]);
    view(-35,50)
    grid on
    hold on
    Target   = scatter3(0, 0, 0, 'filled', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'y', 'MarkerEdgeColor', 'k', 'SizeData', 5^2);
    legend('Chaser Trajectory', 'Target Spacecraft', 'Location','NorthEast')
    if strcmp(print_plots, 'on')
        print -deps NERM_Prop
    end

%EKF Corrected Relative Position
figure(3)
    plot3(EKF_output(:,1)*10^3, EKF_output(:,2)*10^3, EKF_output(:,3)*10^3);
    title('Relative Position Estimates - A Posteriori')
    xlabel('X Position (m)')
    ylabel('Y Position (m)')
    zlabel('Z Position (m)')
    %axis([-200 200 -200 200]);
    grid on
    axis([-2000*max(GPS_output(:,1)) 2000*max(GPS_output(:,1)) -2000*max(GPS_output(:,1)) ...
        2000*max(GPS_output(:,1)) -11000*max(GPS_output(:,3)) 11000*max(GPS_output(:,3))]);
    view(-35,50)
    grid on
    hold on
    Target   = scatter3(0, 0, 0, 'filled', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'y', 'MarkerEdgeColor', 'k', 'SizeData', 5^2);
    legend('Chaser Trajectory', 'Target Spacecraft', 'Location','NorthEast')
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    if strcmp(print_plots, 'on')
        print -deps EKF_Output
    end

%==========================================================================  
%
% Figure: EKF X-Position Error
figure(4)
    %subplot(3,2,1)
    plot(time, errors_x) 
    hold on
    plot(time, 3*sig_x) 
    plot(time, -3*sig_x) 
    legend('EKF Error', '3-\sigma Bounds')
    title('EKF-X Error')
    xlabel('Time (hrs)')
    ylabel('X Position (km)')
    legend('EKF Error', '3-\sigma Bounds')
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    grid on
    %axis([0 time(length(time)) -4 4])

%==========================================================================   
% Figure: EKF Y-Position Error
%subplot(3,2,3)
figure(5)
    plot(time, errors_y) 
    hold on
    plot(time, 3*sig_y) 
    plot(time, -3*sig_y) 
    set(gca,'FontSize',AT_size);
    title('EKF-Y Error')
    xlabel('Time (hrs)')
    ylabel('Y Position (km)')
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
    plot(time, errors_z) 
    hold on
    plot(time, 3*sig_z) 
    plot(time, -3*sig_z) 
    legend('EKF Error', '3-\sigma Bounds')
    set(gca,'FontSize',AT_size);
    title('EKF-Z Error')
    xlabel('Time (hrs)')
    ylabel('Z Position (km)')
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    grid on
    %axis([0 time(length(time)) -1 1])

     
%==========================================================================  
% Figure: EKF X-Velocity Error
figure(7)    
%subplot(3,2,2)
    plot(time, errors_xdot) 
    hold on
    plot(time, 3*sig_xdot) 
    plot(time, -3*sig_xdot) 
    legend('EKF Error', '3-\sigma Bounds')
    set(gca,'FontSize',AT_size);
    title('EKF X-Velocity Error')
    ylabel('X Velocity (km/s)')
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    grid on
    %axis([0 time(length(time)) -.2 .2])

%==========================================================================   
% Figure: EKF Y-Velocity Error
figure(8)    
%subplot(3,2,4)
    plot(time, errors_ydot) 
    hold on
    plot(time, 3*sig_ydot) 
    plot(time, -3*sig_ydot) 
    legend('EKF Error', '3-\sigma Bounds')
    set(gca,'FontSize',AT_size);
    title('EKF Y-Velocity Error')
    ylabel('Y Velocity (km/s)')
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    grid on
    %axis([0 time(length(time)) -.2 .2])

%==========================================================================  
% Figure: EKF Z-Velocity Error
figure(9)    
%subplot(3,2,6)
    plot(time, errors_zdot) 
    hold on
    plot(time, 3*sig_zdot) 
    plot(time, -3*sig_zdot) 
    legend('EKF Error', '3-\sigma Bounds')
    set(gca,'FontSize',AT_size);
    title('EKF Z-Velocity Error')
    ylabel('Z Velocity (km/s)')
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    grid on
    %axis([0 time(length(time)) -0.04 0.04])

%==========================================================================  
% Figure: Remaining corrected state variables
%  
    figure(14)
    clf
    plot(time, 1000*EKF_output(:,1), '-', time, 1000*HCW_output(:,1), '-.', time, 1000*GPS_output(:,1), '--')
    legend('EKF Output', 'Dynamics', 'GPS Measurements')
    set(gca,'FontSize',AT_size);
    title('EKF X Position vs Time')
    ylabel('X Position (m)')
    xlabel('Time (hr)')
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    grid on  
    
    figure(15)
    clf
    plot(time, 1000*EKF_output(:,2), '-', time, 1000*HCW_output(:,2), '-.', time, 1000*GPS_output(:,2), '--')
    legend('EKF Output', 'Dynamics', 'GPS Measurements')
    set(gca,'FontSize',AT_size);
    title('EKF Output Y vs Time')
    xlabel('Time (hr)')
    ylabel('Y Position (m)')
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    grid on  
    
    figure(16)
    clf
    plot(time, 1000*EKF_output(:,3), '-', time, 1000*HCW_output(:,3), '-.', time, 1000*GPS_output(:,3), '--')
    legend('EKF Output', 'Dynamics', 'GPS Measurements')
    set(gca,'FontSize',AT_size);
    title('EKF Output Z vs Time')
    xlabel('Time (hrs)')
    ylabel('Z Position (m)')
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    %axis([0 time(length(time)) -5 5]);
    grid on  

   figure(19)
    clf
    plot(time, EKF_output(:,4)*10^3, '-', time, HCW_output(:,4)*10^3, '-.', time, GPS_output(:,4)*10^3, '--')
    set(gca,'FontSize',AT_size);
    title('EKF Output X-dot vs Time')
    xlabel('Time (hrs)')
    ylabel('X_d_o_t (m/s)')
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    legend('EKF Output', 'Dynamics', 'GPS Measurement')
    grid on  
    
    figure(20)
    clf
    plot(time, EKF_output(:,5)*10^3, '-', time, 1000*HCW_output(:,5), '-.', time, GPS_output(:,5)*10^3, '--')
    set(gca,'FontSize',AT_size);
    title('EKF Output Y-dot vs Time')
    xlabel('Time (hrs)')
    ylabel('Y_d_o_t (m/s)')
    legend('EKF Output', 'Dynamics', 'GPS Measurement')
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    grid on  
    
    figure(21)
    clf
    plot(time, EKF_output(:,6)*10^3, '-', time, 1000*HCW_output(:,6), '-.', time, GPS_output(:,6)*10^3, '--')
    hold on
    set(gca,'FontSize',AT_size);
    title('EKF Output Z-dot vs Time')
    xlabel('Time (hrs)')
    ylabel('Z_d_o_t (m/s)')
    legend('EKF Output', 'Dynamics', 'GPS Measurement')
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    grid on  
        
        
%==========================================================================  
%Subfigures: EKF Outputs
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
% Figure: Covariances
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
%
figure(30)
title('Relative Positions', ... 
    'Fontweight', 'bold', 'FontSize', 12)
curve1 = animatedline('color','r', 'LineWidth', 0.5);
curve2 = animatedline('color','b', 'LineWidth', 0.5);
curve3 = animatedline('color','g', 'LineWidth', 0.5);
set(gca, 'XLim', [-2000*max(GPS_output(:,1)) 2000*max(GPS_output(:,1))], 'YLim', ...
    [-2000*max(GPS_output(:,1)) 2000*max(GPS_output(:,1))], 'ZLim', [-11000*max(GPS_output(:,3)) 11000*max(GPS_output(:,3))]); 
  
xlabel(' x (m) ')
ylabel(' y (m) ')
zlabel(' z (m) ')
view(-40, 30);
set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
grid on
hold on;
EKF  = scatter3(1000*EKF_output(1,1), 1000*EKF_output(1,2), 1000*EKF_output(1,3), 'MarkerEdgeColor', 'r');
HCW  = scatter3(1000*HCW_output(1,1), 1000*HCW_output(1,2), 1000*HCW_output(1,3), 'MarkerEdgeColor', 'b');
GPS  = scatter3(1000*GPS_output(1,1), 1000*GPS_output(1,2), 1000*GPS_output(1,3), 'MarkerEdgeColor', 'b');
legend('EKF', 'Dynamics', 'GPS')
for i = 1:100:30000;
    delete(EKF);
    delete(GPS);
    delete(HCW);
    addpoints(curve1, 1000*EKF_output(i,1), 1000*EKF_output(i,2), 1000*EKF_output(i,3));
    addpoints(curve2, 1000*HCW_output(i,1), 1000*HCW_output(i,2), 1000*HCW_output(i,3));
    addpoints(curve3, 1000*GPS_output(i,1), 1000*GPS_output(i,2), 1000*GPS_output(i,3));
    EKF  = scatter3(1000*EKF_output(i,1), 1000*EKF_output(i,2), 1000*EKF_output(i,3), 'MarkerEdgeColor', 'r');
    HCW  = scatter3(1000*HCW_output(i,1), 1000*HCW_output(i,2), 1000*HCW_output(i,3), 'MarkerEdgeColor', 'b');
    GPS  = scatter3(1000*GPS_output(i,1), 1000*GPS_output(i,2), 1000*GPS_output(i,3), 'MarkerEdgeColor', 'b');
    drawnow
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
%}
% Bar Graph - Position RMS
filter_names = {'GPS', 'HCW EKF'};
RMS_pos_values = [GPS_RMSE_Avg_pos EKF_RMSE_Avg_pos]; 
RMS_vel_values = [GPS_RMSE_Avg_vel EKF_RMSE_Avg_vel]; 

figure(50)
bar(RMS_pos_values*1000)
set(gca,'xticklabel',filter_names);
ylabel('Position Error RMS (m)')
%set(gcf, 'units','normalized','outerposition',[0 0 1 1]);


figure(51)
bar(RMS_vel_values*1000)
set(gca,'xticklabel',filter_names);
ylabel('Velocity Error RMS (m/s)')
%set(gcf, 'units','normalized','outerposition',[0 0 1 1]);

