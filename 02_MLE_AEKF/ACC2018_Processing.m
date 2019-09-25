%% FFNAV_EKF_ACC_Processing ===============================================
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
% Created by: Cory Fraser - Sept 21, 2017
% Latest Edit: Cory Fraser - Sept 21, 2017
%% ========================================================================
%Initialization
clc
close all
fprintf('\n--------------------------------------------------------------\n')
fprintf('   FFNAV EKF SIMULATION POST-PROCESSING                       \n')
fprintf('--------------------------------------------------------------\n')
pathname = fileparts('C:\Users\Cory\Desktop\FFNAV Data\');
% ========================================================================
% No Adaptations
filename = 'FFNAV_PEO_NoAdaptResults.mat';
file_in = fullfile(pathname, filename);
load(file_in)
fprintf('\nParameter file %s has been loaded...\n', filename)

L_NoAdapt = L_sum;
errors_x_NoAdapt = errors_x;
errors_y_NoAdapt = errors_y;
errors_z_NoAdapt = errors_z;
errors_xdot_NoAdapt = errors_xdot;
errors_ydot_NoAdapt = errors_ydot;
errors_zdot_NoAdapt = errors_zdot;
Q_NoAdapt = Q_output;
R_NoAdapt = R_output;

% ========================================================================
% Q Adaptations
filename = 'FFNAV_PEO_QAdaptResults.mat';
file_in = fullfile(pathname, filename);
load(file_in)
fprintf('\nParameter file %s has been loaded...\n', filename)

L_QAdapt = L_sum;
errors_x_QAdapt = errors_x;
errors_y_QAdapt = errors_y;
errors_z_QAdapt = errors_z;
errors_xdot_QAdapt = errors_xdot;
errors_ydot_QAdapt = errors_ydot;
errors_zdot_QAdapt = errors_zdot;
Q_QAdapt = Q_output;
R_QAdapt = R_output;

% ========================================================================
% R Adaptations
filename = 'FFNAV_PEO_RAdaptResults.mat';
file_in = fullfile(pathname, filename);
load(file_in)
fprintf('\nParameter file %s has been loaded...\n', filename)

L_RAdapt = L_sum;
errors_x_RAdapt = errors_x;
errors_y_RAdapt = errors_y;
errors_z_RAdapt = errors_z;
errors_xdot_RAdapt = errors_xdot;
errors_ydot_RAdapt = errors_ydot;
errors_zdot_RAdapt = errors_zdot;
Q_RAdapt = Q_output;
R_RAdapt = R_output;

% ========================================================================
% QR Adaptations
filename = 'FFNAV_PEO_QRAdaptResults.mat';
file_in = fullfile(pathname, filename);
load(file_in)
fprintf('\nParameter file %s has been loaded...\n', filename)

L_QRAdapt = L_sum;
errors_x_QRAdapt = errors_x;
errors_y_QRAdapt = errors_y;
errors_z_QRAdapt = errors_z;
errors_xdot_QRAdapt = errors_xdot;
errors_ydot_QRAdapt = errors_ydot;
errors_zdot_QRAdapt = errors_zdot;
Q_QRAdapt = Q_output;
R_QRAdapt = R_output;


%% ========================================================================
%Plots
print_plots = 'off';
plot_start = 1;
plot_end = 12000;

% Figure: Likelihood Function
figure(1)
    %title('Log-Likelihood Function')
    xlabel('Time (s)')
    y = ylabel('$-\ln L(\theta)$');
    set(y,'Interpreter','latex');
    axis([ 0 time(plot_end) -400 0])
    set(gcf, 'units','normalized','outerposition',[0 0 0.75 1]);
    %grid on    
    hold on
    plot(time(plot_start:plot_end), L_NoAdapt(plot_start:plot_end), ':', 'Color','k','LineWidth',2)
    plot(time(plot_start:plot_end), L_QAdapt(plot_start:plot_end), '-.', 'Color', 'b')
    plot(time(plot_start:plot_end), L_RAdapt(plot_start:plot_end), ':', 'Color', 'r')
    plot(time(plot_start:plot_end), L_QRAdapt(plot_start:plot_end), '-', 'Color', 'k')
    legend('EKF', 'Q-AEKF', 'R-AEKF', 'QR-AEKF', 'Location', 'NorthWest')
if strcmp(print_plots, 'on')
    file_print = fullfile(pathname, 'Fig_Likelihood');
    print(file_print,  '-depsc');
end
    
set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', 32)
set(0,'DefaultTextFontName', 'Times New Roman')
set(0,'DefaultTextFontSize', 32)
set(0,'defaultLineLineWidth', 1.25)

figure(2)
    subplot(3,1,1)  
    hold on
    y = ylabel('$e_x$ (m)');
    set(y,'Interpreter','latex');
    set(gca,'XTick',[]);
    axis([ 0 time(plot_end) -2.5 2.5])
    set(gcf, 'units','normalized','outerposition',[0 0 0.75 1]);
    plot(time(plot_start:plot_end), errors_x_NoAdapt(plot_start:plot_end)*1000, ':', 'Color','k','LineWidth',2) 
    plot(time(plot_start:plot_end), errors_x_QAdapt(plot_start:plot_end)*1000, '-.', 'Color','b') 
    plot(time(plot_start:plot_end), errors_x_RAdapt(plot_start:plot_end)*1000, ':', 'Color','r') 
    plot(time(plot_start:plot_end), errors_x_QRAdapt(plot_start:plot_end)*1000, '-', 'Color','k') 
    %legend('EKF', 'Q-AEKF', 'R-AEKF', 'QR-AEKF', 'Location', 'NorthEast')
    legend({'EKF', 'Q-AEKF', 'R-AEKF', 'QR-AEKF'}, 'Position', [0.75 0.78 0.1 0.1]) 
    
    subplot(3,1,2)  
    hold on
    set(gca,'XTick',[]);
    axis([ 0 time(plot_end) -3 3])
    y = ylabel('$e_y$ (m)');
    set(y,'Interpreter','latex');
    set(gcf, 'units','normalized','outerposition',[0 0 0.75 1]);
    plot(time(plot_start:plot_end), (errors_y_NoAdapt(plot_start:plot_end))*1000, ':', 'Color','k', 'LineWidth',2) 
    plot(time(plot_start:plot_end), (errors_y_QAdapt(plot_start:plot_end))*1000, '-.', 'Color','b') 
    plot(time(plot_start:plot_end), (errors_y_RAdapt(plot_start:plot_end))*1000, ':', 'Color','r') 
    plot(time(plot_start:plot_end), (errors_y_QRAdapt(plot_start:plot_end))*1000, '-', 'Color','k') 
    
    subplot(3,1,3)  
    hold on
    xlabel('Time (s)')
    axis([ 0 time(plot_end) -2 2])
    y = ylabel('$e_z$ (m)');
    set(y,'Interpreter','latex');
    set(gcf, 'units','normalized','outerposition',[0 0 0.75 1]);
    plot(time(plot_start:plot_end), (errors_z_NoAdapt(plot_start:plot_end))*1000, ':', 'Color','k', 'LineWidth',2) 
    plot(time(plot_start:plot_end), (errors_z_QAdapt(plot_start:plot_end))*1000, '-.', 'Color','b') 
    plot(time(plot_start:plot_end), (errors_z_RAdapt(plot_start:plot_end))*1000, ':', 'Color','r') 
    plot(time(plot_start:plot_end), (errors_z_QRAdapt(plot_start:plot_end))*1000, '-', 'Color','k')
if strcmp(print_plots, 'on')
    file_print = fullfile(pathname, 'Fig_PositionErrors');
    print(file_print,  '-depsc');
end   
    
figure(3)
    subplot(3,1,1)  
    hold on
    y = ylabel('$e_{\dot{x}}$ (m/s)');
    set(y,'Interpreter','latex');
    set(gca,'XTick',[]);
    axis([ 0 time(plot_end) -.125 .125])
    set(gcf, 'units','normalized','outerposition',[0 0 0.75 1]);
    plot(time(plot_start:plot_end), errors_xdot_NoAdapt(plot_start:plot_end)*1000, ':', 'Color','k','LineWidth',2) 
    plot(time(plot_start:plot_end), errors_xdot_QAdapt(plot_start:plot_end)*1000, '-.', 'Color','b') 
    plot(time(plot_start:plot_end), errors_xdot_RAdapt(plot_start:plot_end)*1000, ':', 'Color','r') 
    plot(time(plot_start:plot_end), errors_xdot_QRAdapt(plot_start:plot_end)*1000, '-', 'Color','k') 
    %legend('EKF', 'Q-AEKF', 'R-AEKF', 'QR-AEKF', 'Location', 'NorthEast')
    legend({'EKF', 'Q-AEKF', 'R-AEKF', 'QR-AEKF'}, 'Position', [0.75 0.78 0.1 0.1]) 

    
    subplot(3,1,2)  
    hold on
    set(gca,'XTick',[]);
    axis([ 0 time(plot_end) -.125 .125])
    y = ylabel('$e_{\dot{y}}$ (m/s)');
    set(y,'Interpreter','latex');
    set(gcf, 'units','normalized','outerposition',[0 0 0.75 1]);
    plot(time(plot_start:plot_end), (errors_ydot_NoAdapt(plot_start:plot_end))*1000, ':', 'Color','k','LineWidth',2) 
    plot(time(plot_start:plot_end), (errors_ydot_QAdapt(plot_start:plot_end))*1000, '-.', 'Color','b') 
    plot(time(plot_start:plot_end), (errors_ydot_RAdapt(plot_start:plot_end))*1000, ':', 'Color','r') 
    plot(time(plot_start:plot_end), (errors_ydot_QRAdapt(plot_start:plot_end))*1000, '-', 'Color','k') 
    
    subplot(3,1,3)  
    hold on
    xlabel('Time (s)')
    axis([ 0 time(plot_end) -.125 .125])
    y = ylabel('$e_{\dot{z}}$ (m/s)');
    set(y,'Interpreter','latex');
    set(gcf, 'units','normalized','outerposition',[0 0 0.75 1]);
    plot(time(plot_start:plot_end), (errors_zdot_NoAdapt(plot_start:plot_end))*1000, ':', 'Color','k','LineWidth',2) 
    plot(time(plot_start:plot_end), (errors_zdot_QAdapt(plot_start:plot_end))*1000, '-.', 'Color','b') 
    plot(time(plot_start:plot_end), (errors_zdot_RAdapt(plot_start:plot_end))*1000, ':', 'Color','r') 
    plot(time(plot_start:plot_end), (errors_zdot_QRAdapt(plot_start:plot_end))*1000, '-', 'Color','k')

if strcmp(print_plots, 'on')
    file_print = fullfile(pathname, 'Fig_VelocityErrors');
    print(file_print,  '-depsc');
end

%====================================
%Plots
plot_start = 1;
plot_end = floor(length(time)/2);

set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', 32)
set(0,'DefaultTextFontName', 'Times New Roman')
set(0,'DefaultTextFontSize', 32)
set(0,'defaultLineLineWidth', 1.25)

figure(5)
    subplot(3,1,1)  
    hold on
    y = ylabel('$e_x$ (m)');
    set(y,'Interpreter','latex');
    set(gca,'XTick',[]);
    axis([ 0 time(plot_end) -2 2])
    set(gcf, 'units','normalized','outerposition',[0 0 0.75 1]);
    %plot(time(plot_start:plot_end), errors_x_NoAdapt(plot_start:plot_end)*1000, ':', 'Color','k','LineWidth',2) 
    plot(time(plot_start:plot_end), errors_x_QAdapt(plot_start:plot_end)*1000, '-.', 'Color','b') 
   % plot(time(plot_start:plot_end), errors_x_RAdapt(plot_start:plot_end)*1000, ':', 'Color','r') 
    plot(time(plot_start:plot_end), errors_x_QRAdapt(plot_start:plot_end)*1000, '-', 'Color','k') 
    %legend('Q-AEKF', 'QR-AEKF', 'Location', 'NorthEast')
    legend({'Q-AEKF', 'QR-AEKF'}, 'Position', [0.75 0.875 0.1 0.05])  

    %zoomPlot(time(plot_start:plot_end),errors_x_QAdapt(plot_start:plot_end)*1000,[ 0 500],[0.25 0.75 0.2 0.2]) ;
    %zoomPlot(time(plot_start:plot_end),errors_x_QRAdapt(plot_start:plot_end)*1000,[ 0 500],[0.25 0.75 0.2 0.2]) ;
    
    subplot(3,1,2)  
    hold on
    set(gca,'XTick',[]);
    axis([ 0 time(plot_end) -2 2])
    y = ylabel('$e_y$ (m)');
    set(y,'Interpreter','latex');
    set(gcf, 'units','normalized','outerposition',[0 0 0.75 1]);
    %plot(time(plot_start:plot_end), (errors_y_NoAdapt(plot_start:plot_end))*1000, ':', 'Color','k', 'LineWidth',2) 
    plot(time(plot_start:plot_end), (errors_y_QAdapt(plot_start:plot_end))*1000, '-.', 'Color','b') 
    %plot(time(plot_start:plot_end), (errors_y_RAdapt(plot_start:plot_end))*1000, ':', 'Color','r') 
    plot(time(plot_start:plot_end), (errors_y_QRAdapt(plot_start:plot_end))*1000, '-', 'Color','k') 
    
    subplot(3,1,3)  
    hold on
    xlabel('Time (s)')
    axis([ 0 time(plot_end) -2 2])
    y = ylabel('$e_z$ (m)');
    set(y,'Interpreter','latex');
    set(gcf, 'units','normalized','outerposition',[0 0 0.75 1]);
    %plot(time(plot_start:plot_end), (errors_z_NoAdapt(plot_start:plot_end))*1000, ':', 'Color','k', 'LineWidth',2) 
    plot(time(plot_start:plot_end), (errors_z_QAdapt(plot_start:plot_end))*1000, '-.', 'Color','b') 
    %plot(time(plot_start:plot_end), (errors_z_RAdapt(plot_start:plot_end))*1000, ':', 'Color','r') 
    plot(time(plot_start:plot_end), (errors_z_QRAdapt(plot_start:plot_end))*1000, '-', 'Color','k')
if strcmp(print_plots, 'on')
    file_print = fullfile(pathname, 'Fig_PositionErrorsFull');
    print(file_print,  '-depsc');
end   
    
figure(6)
    subplot(3,1,1)  
    hold on
    y = ylabel('$e_{\dot{x}}$ (m/s)');
    set(y,'Interpreter','latex');
    set(gca,'XTick',[]);
    axis([ 0 time(plot_end) -.125 .125])
    set(gcf, 'units','normalized','outerposition',[0 0 0.75 1]);
    %plot(time(plot_start:plot_end), errors_xdot_NoAdapt(plot_start:plot_end)*1000, ':', 'Color','k','LineWidth',2) 
    plot(time(plot_start:plot_end), errors_xdot_QAdapt(plot_start:plot_end)*1000, '-.', 'Color','b') 
    %plot(time(plot_start:plot_end), errors_xdot_RAdapt(plot_start:plot_end)*1000, ':', 'Color','r') 
    plot(time(plot_start:plot_end), errors_xdot_QRAdapt(plot_start:plot_end)*1000, '-', 'Color','k') 
    %legend('Q-AEKF', 'QR-AEKF', 'Location', 'NorthEast')
    legend({'Q-AEKF', 'QR-AEKF'}, 'Position', [0.75 0.875 0.1 0.05])  
    
    subplot(3,1,2)  
    hold on
    set(gca,'XTick',[]);
    axis([ 0 time(plot_end) -.125 .125])
    y = ylabel('$e_{\dot{y}}$ (m/s)');
    set(y,'Interpreter','latex');
    set(gcf, 'units','normalized','outerposition',[0 0 0.75 1]);
    %plot(time(plot_start:plot_end), (errors_ydot_NoAdapt(plot_start:plot_end))*1000, ':', 'Color','k','LineWidth',2) 
    plot(time(plot_start:plot_end), (errors_ydot_QAdapt(plot_start:plot_end))*1000, '-.', 'Color','b') 
   % plot(time(plot_start:plot_end), (errors_ydot_RAdapt(plot_start:plot_end))*1000, ':', 'Color','r') 
    plot(time(plot_start:plot_end), (errors_ydot_QRAdapt(plot_start:plot_end))*1000, '-', 'Color','k') 
    
    subplot(3,1,3)  
    hold on
    xlabel('Time (s)')
    axis([ 0 time(plot_end) -.125 .125])
    y = ylabel('$e_{\dot{z}}$ (m/s)');
    set(y,'Interpreter','latex');
    set(gcf, 'units','normalized','outerposition',[0 0 0.75 1]);
    %plot(time(plot_start:plot_end), (errors_zdot_NoAdapt(plot_start:plot_end))*1000, ':', 'Color','k','LineWidth',2) 
    plot(time(plot_start:plot_end), (errors_zdot_QAdapt(plot_start:plot_end))*1000, '-.', 'Color','b') 
   % plot(time(plot_start:plot_end), (errors_zdot_RAdapt(plot_start:plot_end))*1000, ':', 'Color','r') 
    plot(time(plot_start:plot_end), (errors_zdot_QRAdapt(plot_start:plot_end))*1000, '-', 'Color','k')

if strcmp(print_plots, 'on')
    file_print = fullfile(pathname, 'Fig_VelocityErrorsFull');
    print(file_print,  '-depsc');
end

%Q Covariance Plot
figure(60)
set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    subplot(5,2,1)
    hold on
    plot(time(plot_start:plot_end), squeeze(Q_NoAdapt(1,1,plot_start:plot_end)), ':', 'Color','k', 'LineWidth',2) 
    plot(time(plot_start:plot_end), squeeze(Q_QAdapt(1,1,plot_start:plot_end)), '-.', 'Color','b') 
    plot(time(plot_start:plot_end), squeeze(Q_RAdapt(1,1,plot_start:plot_end)), ':', 'Color','r') 
    plot(time(plot_start:plot_end), squeeze(Q_QRAdapt(1,1,plot_start:plot_end)), '-', 'Color','k') 
    set(gca,'xticklabel',[])
    %axis([ 0 time(plot_end) -.125 .125])
    y = ylabel('$q_{11}$');
    set(y,'Interpreter','latex');

    
    subplot(5,2,2)
    hold on
    plot(time(plot_start:plot_end), squeeze(Q_NoAdapt(6,6,plot_start:plot_end)), ':', 'Color','k', 'LineWidth',2) 
    plot(time(plot_start:plot_end), squeeze(Q_QAdapt(6,6,plot_start:plot_end)), '-.', 'Color','b') 
    plot(time(plot_start:plot_end), squeeze(Q_RAdapt(6,6,plot_start:plot_end)), ':', 'Color','r') 
    plot(time(plot_start:plot_end), squeeze(Q_QRAdapt(6,6,plot_start:plot_end)), '-', 'Color','k') 
    set(gca,'xticklabel',[])
    y = ylabel('$q_{66}$');
    set(y,'Interpreter','latex');

    subplot(5,2,3)
    hold on
    plot(time(plot_start:plot_end), squeeze(Q_NoAdapt(2,2,plot_start:plot_end)), ':', 'Color','k', 'LineWidth',2)
    plot(time(plot_start:plot_end), squeeze(Q_QAdapt(2,2,plot_start:plot_end)), '-.', 'Color','b') 
    plot(time(plot_start:plot_end), squeeze(Q_RAdapt(2,2,plot_start:plot_end)), ':', 'Color','r') 
    plot(time(plot_start:plot_end), squeeze(Q_QRAdapt(2,2,plot_start:plot_end)), '-', 'Color','k')
    set(gca,'xticklabel',[])
    y = ylabel('$q_{22}$');
    set(y,'Interpreter','latex');

    subplot(5,2,4)
    hold on
    plot(time(plot_start:plot_end), squeeze(Q_NoAdapt(7,7,plot_start:plot_end)), ':', 'Color','k', 'LineWidth',2)
    plot(time(plot_start:plot_end), squeeze(Q_QAdapt(7,7,plot_start:plot_end)), '-.', 'Color','b') 
    plot(time(plot_start:plot_end), squeeze(Q_RAdapt(7,7,plot_start:plot_end)), ':', 'Color','r') 
    plot(time(plot_start:plot_end), squeeze(Q_QRAdapt(7,7,plot_start:plot_end)), '-', 'Color','k')
    set(gca,'xticklabel',[])
    y = ylabel('$q_{77}$');
    set(y,'Interpreter','latex');
 
    subplot(5,2,5)
    hold on
    plot(time(plot_start:plot_end), squeeze(Q_NoAdapt(3,3,plot_start:plot_end)), ':', 'Color','k', 'LineWidth',2)
    plot(time(plot_start:plot_end), squeeze(Q_QAdapt(3,3,plot_start:plot_end)), '-.', 'Color','b') 
    plot(time(plot_start:plot_end), squeeze(Q_RAdapt(3,3,plot_start:plot_end)), ':', 'Color','r') 
    plot(time(plot_start:plot_end), squeeze(Q_QRAdapt(3,3,plot_start:plot_end)), '-', 'Color','k')
    set(gca,'xticklabel',[])
    y = ylabel('$q_{33}$');
    set(y,'Interpreter','latex');

    subplot(5,2,6)
    hold on
    plot(time(plot_start:plot_end), squeeze(Q_NoAdapt(8,8,plot_start:plot_end)), ':', 'Color','k', 'LineWidth',2)
    plot(time(plot_start:plot_end), squeeze(Q_QAdapt(8,8,plot_start:plot_end)), '-.', 'Color','b') 
    plot(time(plot_start:plot_end), squeeze(Q_RAdapt(8,8,plot_start:plot_end)), ':', 'Color','r') 
    plot(time(plot_start:plot_end), squeeze(Q_QRAdapt(8,8,plot_start:plot_end)), '-', 'Color','k')
    set(gca,'xticklabel',[])
    y = ylabel('$q_{88}$');
    set(y,'Interpreter','latex');
 
    subplot(5,2,7)
    hold on
    plot(time(plot_start:plot_end), squeeze(Q_NoAdapt(4,4,plot_start:plot_end)), ':', 'Color','k', 'LineWidth',2)
    plot(time(plot_start:plot_end), squeeze(Q_output(4,4,plot_start:plot_end)), '-.', 'Color','b') 
    plot(time(plot_start:plot_end), squeeze(Q_RAdapt(4,4,plot_start:plot_end)), ':', 'Color','r') 
    plot(time(plot_start:plot_end), squeeze(Q_QRAdapt(4,4,plot_start:plot_end)), '-', 'Color','k') 
    set(gca,'xticklabel',[])
    y = ylabel('$q_{44}$');
    set(y,'Interpreter','latex');

    subplot(5,2,8)
    hold on
    plot(time(plot_start:plot_end), squeeze(Q_NoAdapt(9,9,plot_start:plot_end)), ':', 'Color','k', 'LineWidth',2)
    plot(time(plot_start:plot_end), squeeze(Q_QAdapt(9,9,plot_start:plot_end)), '-.', 'Color','b') 
    plot(time(plot_start:plot_end), squeeze(Q_RAdapt(9,9,plot_start:plot_end)), ':', 'Color','r') 
    plot(time(plot_start:plot_end), squeeze(Q_QRAdapt(9,9,plot_start:plot_end)), '-', 'Color','k')
    set(gca,'xticklabel',[])
    y = ylabel('$q_{99}$');
    set(y,'Interpreter','latex');
 
    subplot(5,2,9)
    hold on
    plot(time(plot_start:plot_end), squeeze(Q_NoAdapt(5,5,plot_start:plot_end)), ':', 'Color','k', 'LineWidth',2)
    plot(time(plot_start:plot_end), squeeze(Q_QAdapt(5,5,plot_start:plot_end)), '-.', 'Color','b') 
    plot(time(plot_start:plot_end), squeeze(Q_RAdapt(5,5,plot_start:plot_end)), ':', 'Color','r') 
    plot(time(plot_start:plot_end), squeeze(Q_QRAdapt(5,5,plot_start:plot_end)), '-', 'Color','k')
    %xlabel('Time (hrs)')
    y = ylabel('$q_{55}$');
    set(y,'Interpreter','latex');

    subplot(5,2,10)
    hold on
    plot(time(plot_start:plot_end), squeeze(Q_NoAdapt(10,10,plot_start:plot_end)), ':', 'Color','k', 'LineWidth',2)
    plot(time(plot_start:plot_end), squeeze(Q_QAdapt(10,10,plot_start:plot_end)), '-.', 'Color','b') 
    plot(time(plot_start:plot_end), squeeze(Q_RAdapt(10,10,plot_start:plot_end)), ':', 'Color','r') 
    plot(time(plot_start:plot_end), squeeze(Q_QRAdapt(10,10,plot_start:plot_end)), '-', 'Color','k')
    %xlabel('Time (hrs)')
    y = ylabel('$q_{10,10}$');
    set(y,'Interpreter','latex');
    legend({'EKF', 'Q-AEKF', 'R-AEKF', 'QR-AEKF'}, 'Position', [0.75 0.78 0.1 0.1]) 

    if strcmp(print_plots, 'on')
        file_print = fullfile(pathname, 'Fig_Qcovariances');
        print(file_print,  '-depsc');
    end

figure(61)
set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    subplot(4,2,1)
    hold on
    plot(time(plot_start:plot_end), squeeze(R_NoAdapt(1,1,plot_start:plot_end)), ':', 'Color','k', 'LineWidth',2)
    plot(time(plot_start:plot_end), squeeze(R_QAdapt(1,1,plot_start:plot_end)), '-.', 'Color','b') 
    plot(time(plot_start:plot_end), squeeze(R_RAdapt(1,1,plot_start:plot_end)), ':', 'Color','r') 
    plot(time(plot_start:plot_end), squeeze(R_QRAdapt(1,1,plot_start:plot_end)), '-', 'Color','k') 
    set(gca,'xticklabel',[])
    y = ylabel('$r_{11}$');
    set(y,'Interpreter','latex');

    subplot(4,2,2)
    hold on
    plot(time(plot_start:plot_end), squeeze(R_NoAdapt(5,5,plot_start:plot_end)), ':', 'Color','k', 'LineWidth',2) 
    plot(time(plot_start:plot_end), squeeze(R_QAdapt(5,5,plot_start:plot_end)), '-.', 'Color','b') 
    plot(time(plot_start:plot_end), squeeze(R_RAdapt(5,5,plot_start:plot_end)), ':', 'Color','r') 
    plot(time(plot_start:plot_end), squeeze(R_QRAdapt(5,5,plot_start:plot_end)), '-', 'Color','k') 
    set(gca,'xticklabel',[])
    y = ylabel('$r_{55}$');
    set(y,'Interpreter','latex');
    
    
    subplot(4,2,3)
    hold on
    plot(time(plot_start:plot_end), squeeze(R_NoAdapt(2,2,plot_start:plot_end)), ':', 'Color','k', 'LineWidth',2) 
    plot(time(plot_start:plot_end), squeeze(R_QAdapt(2,2,plot_start:plot_end)), '-.', 'Color','b') 
    plot(time(plot_start:plot_end), squeeze(R_RAdapt(2,2,plot_start:plot_end)), ':', 'Color','r') 
    plot(time(plot_start:plot_end), squeeze(R_QRAdapt(2,2,plot_start:plot_end)), '-', 'Color','k') 
    set(gca,'xticklabel',[])
    y = ylabel('$r_{22}$');
    set(y,'Interpreter','latex');
    
    subplot(4,2,4)
    hold on
    plot(time(plot_start:plot_end), squeeze(R_NoAdapt(6,6,plot_start:plot_end)), ':', 'Color','k', 'LineWidth',2) 
    plot(time(plot_start:plot_end), squeeze(R_QAdapt(6,6,plot_start:plot_end)), '-.', 'Color','b') 
    plot(time(plot_start:plot_end), squeeze(R_RAdapt(6,6,plot_start:plot_end)), ':', 'Color','r') 
    plot(time(plot_start:plot_end), squeeze(R_QRAdapt(6,6,plot_start:plot_end)), '-', 'Color','k') 
    set(gca,'xticklabel',[])
    y = ylabel('$r_{66}$');
    set(y,'Interpreter','latex');
    
    subplot(4,2,5)
    hold on
    plot(time(plot_start:plot_end), squeeze(R_NoAdapt(3,3,plot_start:plot_end)), ':', 'Color','k', 'LineWidth',2) 
    plot(time(plot_start:plot_end), squeeze(R_QAdapt(3,3,plot_start:plot_end)), '-.', 'Color','b') 
    plot(time(plot_start:plot_end), squeeze(R_RAdapt(3,3,plot_start:plot_end)), ':', 'Color','r') 
    plot(time(plot_start:plot_end), squeeze(R_QRAdapt(3,3,plot_start:plot_end)), '-', 'Color','k') 
    set(gca,'xticklabel',[])
    y = ylabel('$r_{33}$');
    set(y,'Interpreter','latex');
    
    subplot(4,2,6)
    hold on
    plot(time(plot_start:plot_end), squeeze(R_NoAdapt(7,7,plot_start:plot_end)), ':', 'Color','k', 'LineWidth',2) 
    plot(time(plot_start:plot_end), squeeze(R_QAdapt(7,7,plot_start:plot_end)), '-.', 'Color','b') 
    plot(time(plot_start:plot_end), squeeze(R_RAdapt(7,7,plot_start:plot_end)), ':', 'Color','r') 
    plot(time(plot_start:plot_end), squeeze(R_QRAdapt(7,7,plot_start:plot_end)), '-', 'Color','k') 
    set(gca,'xticklabel',[])
    y = ylabel('$r_{77}$');
    set(y,'Interpreter','latex');
    
    subplot(4,2,7)
    hold on
    plot(time(plot_start:plot_end), squeeze(R_NoAdapt(4,4,plot_start:plot_end)), ':', 'Color','k', 'LineWidth',2) 
    plot(time(plot_start:plot_end), squeeze(R_QAdapt(4,4,plot_start:plot_end)),'-.', 'Color','b') 
    plot(time(plot_start:plot_end), squeeze(R_RAdapt(4,4,plot_start:plot_end)), ':', 'Color','r') 
    plot(time(plot_start:plot_end), squeeze(R_QRAdapt(4,4,plot_start:plot_end)), '-', 'Color','k') 
    set(gca,'xticklabel',[])
    y = ylabel('$r_{44}$');
    set(y,'Interpreter','latex');
    legend({'EKF', 'Q-AEKF', 'R-AEKF', 'QR-AEKF'}, 'Position', [0.7 0.125 0.1 0.1]) 

    if strcmp(print_plots, 'on')
        file_print = fullfile(pathname, 'Fig_Rcovariances');
        print(file_print,  '-depsc');
    end