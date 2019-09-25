function [FUZZY_props] = FFNAV_FAEKF_FuzzyInitialize(Plot_flag)
% FFNAV Fuzzy Adaptation Initialization ===================================
% Description: This script defines the linguistic variables, membership 
% functions, universes of discourse, and implication methods for the FLS.
% Plots of the I/O membership functions can be created based on user input. 
%
% Inputs: 
%   Plot_flag - Flag for membership function plots (1 = make plots)
%
% Outputs:
%   FUZZY_props - structure containing various properties
%
% Other Function Calls:
%   None
%
% Created by: Cory Fraser  - NOV 15, 2017
% Latest Edit: Cory Fraser - JUL 04, 2018
% Copyright(c) 2018 by Cory Fraser
%==========================================================================

%% Linguistic Values

%Inputs NEG HI, NEG LO, ZERO, POS LO, POS HI
LingValsIn = {'NH', 'NL','ZERO', 'PL', 'PH'};
n_LingValsIn = length(LingValsIn);

%SISO Outputs ZERO, MIN, MAX
%LingValsOut_SISO = {'ZERO', 'MIN', 'MAX'};
%n_LingValsOut_SISO = length(LingValsOut_SISO);

%SISO Outputs 
LingValsOut_SISO = {'NMAX', 'NMIN','ZERO', 'PMIN', 'PMAX'};
n_LingValsOut_SISO = length(LingValsOut_SISO);

%MISO Outputs ZERO, PL, MED, PH, MAX
LingValsOut_MISO = {'ZERO', 'LOW', 'MED', 'HIGH', 'MAX'};
n_LingValsOut_MISO = length(LingValsOut_MISO);

%==========================================================================
%% Rule Bases

%SISO System (1 = NEG HI, 2 = NEG LO, 3 = ZERO, 4 = POS LO, 5 = POS HI)
rules_SISO    = [5 4 3 2 1];
n_Rules_SISO  = length(rules_SISO);

%MISO System (1 = ZERO, 2 = LOW, 3 = MED, 4 = HIGH, 5 = MAX)
rules_MISO   = [ 5 4 3 4 5
                 4 3 2 3 4
                 3 2 1 2 3
                 4 3 2 3 4
                 5 4 3 4 5 ];
n_Rules_MISO  = numel(rules_MISO);        

%==========================================================================
%% Centers of the Membership Functions

%        [  NH     NL   ZERO   PL    PH ]
%c_IN =  [  -1   -0.55   0    0.55   1  ];    %Spread - Outwards
 c_IN =  [  -1   -0.5    0    0.5    1  ];    %Spread - Equal
%c_IN =  [  -1   -0.25   0    0.25   1  ];    %Spread - Inwards

%          [  NH       NL   ZERO   PL     PH   ]
%c_SISO =  [  -1     -0.75    0   0.75     1   ];   %Spread - Outwards
%c_SISO =  [  -1     -0.5     0   0.5      1   ];   %Spread - Equal
c_SISO =  [  -1     -0.25    0   0.25     1   ];   %Spread - Inwards

% Calculating the spread (inwards)
c_new = sign(c_SISO).*c_SISO.^2;

%          [  ZERO    PL      M       PH    MAX  ]
%c_MISO =  [   0     0.5    0.75     0.9    1    ];   %Spread - Higher
 c_MISO =  [   0     0.25   0.5     0.75    1    ];   %Spread - Equal
%c_MISO =  [   0     0.1    0.25     0.5    1    ];   %Spread - Lower

%==========================================================================
%% Membership Function Definitions

n_pts      = 101;                   %Number of points
X_pts      = linspace(-1,1,n_pts);   %Universe of discourse - Inputs
Y_pts      = linspace(-1,1,n_pts);   %Universe of discourse - Outputs

sigma_IN   = 1/12;                   %Gaussian Curvature - Inputs
sigma_SISO = 1/12;                   %Gaussian Curvature - Output
sigma_MISO = 1/12;                   %Gaussian Curvature - Output

sigmoid_a   = 25;                   %Sigmoidal Curvature
sigmoid_c   = 0.75;                 %Sigmoidal Center

%==========================================================================
% Input Membership Functions
in_mf = zeros(n_LingValsIn,n_pts);

% Outer Functions - Sigmoidal for Saturation
in_mf(1,1:n_pts) = 1./(1+exp(sigmoid_a*(X_pts+sigmoid_c)));;
in_mf(n_LingValsIn,1:n_pts) = 1./(1+exp(-sigmoid_a*(X_pts-sigmoid_c)));

% Inner Functions - Gaussian 
for i = 2:n_LingValsIn-1
    in_mf(i,1:n_pts) = exp(-((X_pts-c_IN(i)).^2)/(2*sigma_IN^2));
end

%==========================================================================
% SISO System Output Membership Functions
out_mf_SISO = zeros(n_LingValsOut_SISO,n_pts);

% Gaussian MSFs
for i = 1:n_LingValsOut_SISO
    out_mf_SISO(i,1:n_pts) = exp(-((Y_pts-c_SISO(i)).^2)/(2*sigma_SISO^2));
end

%==========================================================================
% MISO System Output Membership Function
out_mf_MISO = zeros(n_LingValsOut_MISO,n_pts);

% Gaussian MSFs
for i = 1:n_LingValsOut_MISO
    out_mf_MISO(i,1:n_pts) = exp(-((Y_pts-c_MISO(i)).^2)/(2*sigma_MISO^2));
end

%==========================================================================
% Plot Membership Functions
if (Plot_flag)
    
    TF_size       = 28;  %Title font size
    AF_size       = 28;  %Axes label font size
    AT_size       = 28;  %Tick label font size
    Line_width    = 2;   %Line width

    set(0,'DefaultAxesFontSize', AF_size)
    set(0,'DefaultTextFontSize', TF_size)
    set(0,'defaultLineLineWidth', Line_width)
    
    %fprintf(' - Plotting Membership Functions \n')

    %Figure - Subplots of membership function
    fig = figure(90);
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
    
    set(gcf, 'units','normalized','outerposition',[0.1 0.1 0.75 0.8]);
    subplot(2,1,1)
    plot(X_pts, in_mf, '-k')
    set(gca,'xtick',c_IN);
    t1 = text(c_IN-[0.00, 0.04, 0.08, 0.04, 0.09],1.15*ones(1,n_LingValsIn), LingValsIn);
    set(t1,'Interpreter','latex');
    t2 = text([-1.1 1.0],[1 1],{'1.0'});
    set(t2,'Interpreter','latex');
    axis([ -1.2 1.2 0 1.2])
    ax = gca; % Set the axes properties
    ax.FontName = 'LaTeX';
    ax.TickLabelInterpreter = 'LaTeX';
    ax.YLim = [0, 1.25];
    ax.YTick = 0:1;
    ax.XLim = [-1, 1];
    ax.XColor = 'black';
    ax.YColor = 'white';
    ax.Box = 'off';
    ax.Color = [255, 255, 255]/255;
    ax.YGrid = 'off';
    ax.XGrid = 'off';
    ylab = ylabel(sprintf('Input \n Membership \n Functions'));
    set(get(gca,'YLabel'),'Rotation',0)
    ylab.Color = 'black';
    ylab.Interpreter = 'LaTeX';
    ylab.FontSize = 28;
    set(ylab, 'Units', 'Normalized', 'Position', [-0.05, 0.175, 0]);

    subplot(2,1,2)
    plot(Y_pts, out_mf_SISO, '-k')
    set(gca,'xtick',c_SISO);
    t1 = text(c_SISO-[0.0, 0.1, 0.09, 0.09, 0.21],1.15*ones(1,n_LingValsOut_SISO),LingValsOut_SISO);
    set(t1,'Interpreter','latex');
    t2 = text([-1.1 1.0],[1 1],{'1.0'});
    set(t2,'Interpreter','latex');
    axis([ Y_pts(1)-0.2 Y_pts(end)+0.2 0 1.2])
    ax = gca; % Set the axes properties
    ax.FontName = 'LaTeX';
    ax.TickLabelInterpreter = 'LaTeX';
    ax.YLim = [0, 1.25];
    ax.YTick = 0:1;
    ax.XLim = [-1, 1];
    ax.XColor = 'black';
    ax.YColor = 'white';
    ax.Box = 'off';
    ax.Color = [255, 255, 255]/255;
    ax.YGrid = 'off';
    ax.XGrid = 'off';
    ylab = ylabel(sprintf('Output \n Membership \n Functions'));
    set(get(gca,'YLabel'),'Rotation',0)
    ylab.Color = 'black';
    ylab.Interpreter = 'LaTeX';
    ylab.FontSize = 28;
    set(ylab, 'Units', 'Normalized', 'Position', [-0.05, 0.175, 0]);
    %}
    
    %Figure - Input Membership Functions
    fig = figure(91);
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
    
    set(gcf, 'units','normalized','outerposition',[0.1 0.1 0.75 0.8]);
    plot(X_pts, in_mf, '-k')
    set(gca,'xtick',c_IN);
    t1 = text(c_IN-[0.00, 0.04, 0.08, 0.04, 0.09],1.05*ones(1,n_LingValsIn), LingValsIn);
    set(t1,'Interpreter','latex');
    t2 = text([-1.1 1.0],[1 1],{'1.0'});
    set(t2,'Interpreter','latex');
    axis([ -1.2 1.2 0 1.2])
    ax = gca; % Set the axes properties
    ax.FontName = 'LaTeX';
    ax.TickLabelInterpreter = 'LaTeX';
    ax.YLim = [0, 1.25];
    ax.YTick = 0:1;
    ax.XLim = [-1, 1];
    ax.XColor = 'black';
    ax.YColor = 'white';
    ax.Box = 'off';
    ax.Color = [255, 255, 255]/255;
    ax.YGrid = 'off';
    ax.XGrid = 'off';
    %t=title('Input Membership Functions');
    %t.Color = 'black';
    %t.Interpreter = 'LaTeX';
    %t.FontSize = 24;
    ylab = ylabel(sprintf('Input \n Membership \n Functions'));
    set(get(gca,'YLabel'),'Rotation',0)
    ylab.Color = 'black';
    ylab.Interpreter = 'LaTeX';
    ylab.FontSize = 28;
    set(ylab, 'Units', 'Normalized', 'Position', [-0.05, 0.3, 0]);

    
%Figure - Output Membership Functions
    fig = figure(92);
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
    set(gcf, 'units','normalized','outerposition',[0.1 0.1 0.75 0.8]);
    
    plot(Y_pts, out_mf_SISO, '-k')
    set(gca,'xtick',c_SISO);
    t1 = text(c_SISO-[0.0, 0.1, 0.09, 0.09, 0.21],1.05*ones(1,n_LingValsOut_SISO),LingValsOut_SISO);
    set(t1,'Interpreter','latex');
    t2 = text([-1.1 1.0],[1 1],{'1.0'});
    set(t2,'Interpreter','latex');
    axis([ Y_pts(1)-0.2 Y_pts(end)+0.2 0 1.2])
    ax = gca; % Set the axes properties
    ax.FontName = 'LaTeX';
    ax.TickLabelInterpreter = 'LaTeX';
    ax.YLim = [0, 1.25];
    ax.YTick = 0:1;
    ax.XLim = [-1, 1];
    ax.XColor = 'black';
    ax.YColor = 'white';
    ax.Box = 'off';
    ax.Color = [255, 255, 255]/255;
    ax.YGrid = 'off';
    ax.XGrid = 'off';
    %t=title('Output Membership Functions');
    %t.Color = 'black';
    %t.Interpreter = 'LaTeX';
    %t.FontSize = 24;
    ylab = ylabel(sprintf('Output \n Membership \n Functions'));
    set(get(gca,'YLabel'),'Rotation',0)
    ylab.Color = 'black';
    ylab.Interpreter = 'LaTeX';
    ylab.FontSize = 28;
    set(ylab, 'Units', 'Normalized', 'Position', [-0.05, 0.3, 0]);
    %}
    
    
    %{
    figure(93)
    plot(Y_pts, out_mf_MISO, '-k')
    title('Output Membership Functions (MISO)')
    set(gca,'xtick',c_MISO);
    text(c_MISO-0.05,1.1*ones(1,n_LingValsOut_MISO), LingValsOut_MISO)
    axis([ -0.2 1.2 0 1.2])
    hold on
    %}
    
end

%==========================================================================
%% Initialize Output Fuzzy Sets

% 1 = Implied Fuzzy Sets (Faster)
% 0 = Overall (Aggregated) Fuzzy Set

implied_flag = 1;

if implied_flag == 1
    implied_mf_MISO = zeros(1,n_pts);
    implied_mf_SISO = implied_mf_MISO;
    aggregate_mf = 0;
    fprintf(' Implied Fuzzy Sets \n')
else
    implied_mf_MISO = zeros(n_LingValsIn^2,n_pts);
    implied_mf_SISO = zeros(n_LingValsIn,n_pts);
    aggregate_mf = zeros(1,n_pts);
    fprintf(' Aggregated Fuzzy Sets \n')
end

%==========================================================================
% Form a structure for Fuzzy Properties
FUZZY_props.rules_MISO    = rules_MISO;
FUZZY_props.rules_SISO    = rules_SISO;
FUZZY_props.n_pts         = n_pts;
FUZZY_props.X_pts         = X_pts;
FUZZY_props.Y_pts         = Y_pts;

% Membership Function Definitions
FUZZY_props.c_IN           = c_IN;
FUZZY_props.c_SISO         = c_SISO;
FUZZY_props.c_MISO         = c_MISO;
FUZZY_props.sigmoid_a      = sigmoid_a;
FUZZY_props.sigmoid_c      = sigmoid_c;

FUZZY_props.implied_flag   = implied_flag;
FUZZY_props.sigma_IN       = sigma_IN;
FUZZY_props.sigma_SISO     = sigma_SISO;
FUZZY_props.sigma_MISO     = sigma_MISO;

% Membership Functions
FUZZY_props.out_mf_SISO     = out_mf_SISO;
FUZZY_props.out_mf_MISO     = out_mf_MISO;
FUZZY_props.implied_mf_SISO = implied_mf_SISO;
FUZZY_props.implied_mf_MISO = implied_mf_MISO;
FUZZY_props.aggregate_mf    = aggregate_mf;

% Number of rules and lingustic variables
FUZZY_props.n_LingValsIn        = n_LingValsIn;
FUZZY_props.n_LingValsOut_SISO  = n_LingValsOut_SISO;
FUZZY_props.n_LingValsOut_MISO  = n_LingValsOut_MISO;

save('FAEKF_FuzzyProps.mat', 'FUZZY_props');
end
%==========================================================================