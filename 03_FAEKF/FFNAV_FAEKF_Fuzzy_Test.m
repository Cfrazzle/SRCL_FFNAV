% FFNAV Fuzzy Adaptation Testbed ==========================================
% Description: This script is used to test the various fuzzy adaptation
% functions. A plot of the SISO-FLS control curve is created, and plots of
% the MISO-FLS control surface and contour are produced.
%
% Inputs:
%   None
%
% Outputs:
%   SISO Control Curve Plot
%   MISO Control Surface Plot
%   MISO Control Contour Plot

% Other Function Calls:
%   FFNAV_EKF_Fuzzy_Type1_SISO.m
%          - Fuzzy inference for the single-input single-output system
%   FFNAV_EKF_Fuzzy_Type1_MISO.m
%          - Fuzzy inference for the multi-input single-output system
%
% Created by:  Cory Fraser - JAN 09, 2018
% Latest Edit: Cory Fraser - MAY 13, 2018
% Copyright(c) 2018 by Cory Fraser
% =========================================================================

%% Fuzzy Adaptations - Testing =====================================
clc; clear; close all;

%Create membership functions and fuzzy properties
FFNAV_FAEKF_FuzzyInitialize(1); % 1 = show membership function plots
load('FAEKF_FuzzyProps.mat', 'FUZZY_props')

% Implement SISO Fuzzy Inference System
test_range = 1.1;
u_SISO = [-test_range:0.1:test_range];

for i = 1:length(u_SISO)
    [lambda_SISO(i)] = FFNAV_FAEKF_Fuzzy_Type1_SISO(u_SISO(i), FUZZY_props);
end

% Implement MISO Fuzzy Inference System
[u1_MISO,u2_MISO] = meshgrid(-test_range:0.1:test_range); 
for i = 1:1:length(u1_MISO)
    for j = 1:1:length(u1_MISO)
        lambda_MISO(i,j) = FFNAV_FAEKF_Fuzzy_Type1_MISO([u1_MISO(i,j) u2_MISO(i,j)], FUZZY_props);
    end
end


    TF_size       = 28;  %Title font size
    AF_size       = 28;  %Axes label font size
    AT_size       = 28;  %Tick label font size
    Line_width    = 2;   %Line width

    set(0,'DefaultAxesFontSize', AF_size)
    set(0,'DefaultTextFontSize', TF_size)
    set(0,'defaultLineLineWidth', Line_width)
    
%Figure - SISO Control Curve
fig = figure(10);
% Set the paper size
% fig.Resize = 'off';
fig.PaperUnits = 'inches';
fig.Units = 'inches';
fig.PaperPositionMode = 'manual';
fig.PaperPosition = [0, 0, 10, 5];
fig.PaperSize = [10, 5];
fig.Position = [0, 0, 9.9, 4.9];
fig.Color = [255,255,255]/255; % Background color
fig.InvertHardcopy = 'off'; % Prevent the background color from chaning on save
hold on
set(gcf, 'units','normalized','outerposition',[0.1 0.1 0.75 0.8]);
plot(u_SISO, lambda_SISO, '-k')
ax = gca; % Set the axes properties
ax.FontName = 'LaTeX';
ax.TickLabelInterpreter = 'LaTeX';
ax.YLim = [-1.25, 1.25];
%ax.YTick = 0:1;
ax.XLim = [-1.25, 1.25];
ax.XColor = 'black';
ax.YColor = 'black';
ax.Box = 'off';
ax.Color = [255, 255, 255]/255;
ax.YGrid = 'off';
ax.XGrid = 'off';
%t=title('Fuzzy Control Curve (Normalized)');
%t.Color = 'black';
%t.Interpreter = 'LaTeX';
%t.FontSize = 24;
ylab = ylabel(('Output $\lambda_i$'));
%set(get(gca,'YLabel'),'Rotation',0)
ylab.Color = 'black';
ylab.Interpreter = 'LaTeX';
ylab.FontSize = 36;
%set(ylab, 'Units', 'Normalized', 'Position', [-0.05, 0.175, 0]);
xlab = xlabel(('Input $u_i$'));
xlab.Color = 'black';
xlab.Interpreter = 'LaTeX';
xlab.FontSize = 36;
%axis square
%set(xlab, 'Units', 'Normalized', 'Position', [-0.05, 0.175, 0]);


%MISO Control Surface
figure(13)
surf(u1_MISO,u2_MISO,lambda_MISO)
%surf(X,Y,Z,'FaceColor','red','EdgeColor','none')
title('Fuzzy Control Surface (Normalized)')
xlabel('Input u_1')
ylabel('Input u_2')
%colormap hsv
h = colorbar;
title(h,'Output \lambda')
caxis([0 1])
camlight left; 
lighting phong

figure(14)
contourf(u1_MISO,u2_MISO,lambda_MISO)
title('Fuzzy Control Contour (Normalized)')
xlabel('Input u_1')
ylabel('Input u_2')
h = colorbar;
title(h,'Output \lambda')
caxis([0 1])
%colormap hsv

%==========================================================================