function [Q_update, R_update, NIS,deltaDOD,DOD_true,epsilon_Q,epsilon_R] = ...
    FFNAV_FAEKF_FuzzyAdapt(state_pre,Pr_theo, Pr_smooth,residuals,...
    time_step, mu, Q0, R0, Q_pre, R_pre, time, FUZZY_props,Q_adapt_flag, R_adapt_flag)
% FFNAV Fuzzy Adaptations ============================================
% Description: This function adapts the process noise covariance Q and
% measurement noise covariance R matrices using a Fuzzy Logic System. 
%
% Inputs:
%   state_pre    - Previous state estimate
%   Pr_theo      - Theoretical covariance of the innovations
%   Pr_smooth    - Smoothed sampled covariance of the innovations
%   residuals    - Original innovations from the current time step
%   time_step    - Time step of the simulation
%   mu           - Earth's gravitational parameter
%   Q0           - Initial process noise covariance matrix
%   R0           - Initial measurement noise covariance matrix
%   Q_pre        - Previous time step process noise covariance matrix
%   R_pre        - Previous time step measurement noise covariance matrix
%   time         - Current simulation time
%   Fuzzy_props. - Structure containing data for the Fuzzy System
%   Q_adapt_flag - Flag to adapt Q (0 = no adaptation) 
%   R_adapt_flag - Flag to adapt R (0 = no adaptation)
%
% Outputs:
%   Q_update    - Updated process noise covariance matrix
%   R_update    - Updated measurement noise covariance matrix
%   NIS         - Normalized innovation squared
%   deltaDOD    - Difference in Degree of Divergence
%   DOD_true    - True degree of Divergence
%   epsilonQ    - Scaling parameter for Q matrix
%   epsilonR    - Scaling parameters for R matrix (7x7 matrix)
%
% Other Functions Called:
%   FFNAV_FAEKF_Fuzzy_Type1_SISO - Computes fuzzy output from single input
%   FFNAV_FAEKF_Fuzzy_Type1_MISO - Computes fuzzy output from multiple inputs
%   findPSDS - Finds nearest positive semidefinite symmetric matrix
%   
% References:
%   -
%   -
%
% Notes:
%   1) 
%
% Created by:  Cory Fraser - JAN 13, 2018
% Latest Edit: Cory Fraser - SEP 11, 2018
% Copyright(c) 2018 by Cory Fraser

%% EKF PERFORMANCE METRICS ================================================

%Innovations Covariances
Pr_theo = diag(Pr_theo);            %Theoretical Covariance of Residuals
Pr_smooth = diag(Pr_smooth);        %Smoothed Covariance of Residuals
DOM = Pr_theo - Pr_smooth;          %Degree of Mismatch

%DOD_true = residuals'*residuals;   %Degree of Divergence - True
DOD_true = trace(Pr_smooth);        %Degree of Divergence - Smoothed True
DOD_theo = trace(Pr_theo);          %Degree of Divergence - Theoretical
deltaDOD  = DOD_true - DOD_theo;    %Difference in DOD

%Normalized Innovation Squared (NIS) - Divergence Ratio
NIS_threshold = 0;                          %NIS Threshold for fault
%NIS = residuals'*inv(Pr_theo)*residuals;   %NIS w/ current residuals
Pr_inv = Pr_theo\eye(7); 
NIS = residuals'*Pr_inv*residuals;    
NIS = DOD_true/DOD_theo;                    %NIS w/ current residuals

%% FUZZY ADAPTATIONS ======================================================

if (Q_adapt_flag)
    
% Scalar Q-Adaptations, using Degree of Divergence

    %PRISMA Gains
    %{
    g_DOD  = -5e-3; % Input Scaling Gain (PRISMA)
    h_DOD  = 1e-2;  % Output Scaling Gain (PRISMA)
    %}
    
    %PROBA-3 Gains
    g_DOD  = -3e-3; % Input Scaling Gain (PRISMA)
    h_DOD  = 1e-2;  % Output Scaling Gain (PRISMA)
    %g_DOD  = -0.33e-3; % Input Scaling Gain (PROBA-3)
    %h_DOD  = 1e-3;  % Output Scaling Gain (PROBA-3)
 
    %FLS Inputs
    u_DOD       = g_DOD*deltaDOD;   % FIS Inputs

    % Single-Input Single-Output Fuzzy Inference System
    lambda_DOD  = FFNAV_FAEKF_Fuzzy_Type1_SISO(u_DOD, FUZZY_props);   

    epsilon_Q   = 1+lambda_DOD*h_DOD; % Q-Scale Factor (make sure > 0)
    Q_update    = epsilon_Q*Q_pre;    % Updated Q

else
    Q_update    = Q_pre;
    epsilon_Q   = 1;
end
%=========================================================================
if (R_adapt_flag)

    % Additive R-Adaptations, using Degree of Mismatch
    %{
    % Input Scaling Gains (Original)
    g_r     = 5*10^-2;     %Gain for Position
    g_theta = 1*10^1;      %Gain for True Anomaly
    g_v     = 1*10^1;      %Gain for Velocity

    % Input Scaling Gains
    %g_r     = 5*10^-2;     %Gain for Position
    %g_theta = 1*10^1;      %Gain for True Anomaly
    %g_v     = 1*10^0;      %Gain for Velocity

    G = [g_r g_r g_r g_theta g_v g_v g_v]';

    % FIS Inputs
    u_DOM = G.*diag(DOM);

    % Single-Input Single-Output Fuzzy Inference System
    lambda = zeros(7,1);
    for i = 1:length(u_DOM)
        lambda(i) = FFNAV_FAEKF_Fuzzy_Type1_SISO(u_DOM(i), FUZZY_props);
    end

    %Output Scaling Gains (Original)
    h_r     = 10^-2;       %Gain for Position
    h_theta = 10^-4;       %Gain for True Anomaly
    h_v     = 5*10^-5;     %Gain for Velocity

    %Output Scaling Gains
    %h_r     = 10^0;       %Gain for Position
    %h_theta = 10^-4;       %Gain for True Anomaly
    %h_v     = 5*10^-5;     %Gain for Velocity

    H = [h_r h_r h_r h_theta h_v h_v h_v]';

    % R Additive Factor
    deltaR = diag(H.*lambda);
    R_update = R_pre + deltaR;
    epsilon_R = deltaR./R_pre;
    %}
    
    %=========================================================================
    % Scalar R-Adaptations, using Degree of Mismatch

    %Gains for PRISMA
    %{
    % Input Scaling Gains (Original - PRISMA)
    g_r     = 5e-2;     %Gain for Position
    g_theta = 1e1;      %Gain for True Anomaly
    g_v     = 1e0;      %Gain for Velocity
    
    %Output Scaling Gains (Original - PRISMA)
    h_r     = 1e-4;       %Gain for Position
    h_theta = 1e-4;       %Gain for True Anomaly
    h_v     = 1e-3;       %Gain for Velocity
    %}    
 
    %Gains for PROBA-3
    %{
    % Input Scaling Gains (Original - PROBA-3)
    g_r     = 5e-2;     %Gain for Position
    g_theta = 1e0;      %Gain for True Anomaly
    g_v     = 1e0;      %Gain for Velocity
    
    %Output Scaling Gains (Original - PROBA-3)
    h_r     = 1e-2;       %Gain for Position
    h_theta = 1e-3;       %Gain for True Anomaly
    h_v     = 1e-2;       %Gain for Velocity
    %}  
    %
    g_r     = 1e-2;     %Gain for Position
    g_theta = 0.5e3;      %Gain for True Anomaly
    g_v     = 0.5e0;      %Gain for Velocity
    
    h_r     = 0.5e-4;       %Gain for Position
    h_theta = 0.5e-4;       %Gain for True Anomaly
    h_v     = 0.5e-3;       %Gain for Velocity
    %}
    G = [g_r g_r g_r g_theta g_v g_v g_v]';
    H = [h_r h_r h_r h_theta h_v h_v h_v]';

    % FLS Inputs
    u_DOM = G.*diag(DOM);
    %if u_DOM(1) < 0
    %    temp = 2;
    %end        
        
    % Single-Input Single-Output Fuzzy Inference System
    lambda = zeros(7,1);
    for i = 1:length(u_DOM)
        lambda(i) = FFNAV_FAEKF_Fuzzy_Type1_SISO(u_DOM(i), FUZZY_props);
    end
   
    %Scale Factor and Update
    epsilon_R = eye(7)+diag(H.*lambda); % R Scaling Factor
    R_update = epsilon_R.*R_pre;   % Updated R
    
else
    R_update = R_pre;
    epsilon_R = eye(7,7);
end
%=======================================
% Check and correct for Positive Semidefinite Symmetric (PSDS)

%Q_update = findPSDS(Q_update);
%R_update = findPSDS(R_update);
        
% =========================================================================       

% Multi-Input Single-Output Fuzzy Inference System
%lambda2 = FFNAV_FAEKF_Fuzzy_Type1_MISO([-1 -1], FUZZY_props);

%==========================================================================
%% Weighted Kalman Filter (Fixed Weights)
%{
k = time/time_step;

alpha = 1.01;
beta = 1.0005;

%Q_adapt = Q0*alpha^(-2*(k+1)); Don't like this, k+1 means Q0 is not Q0
Q_adapt = Q0;%*alpha^(-2*(k));
R_adapt = R0*beta^(-2*(k));
%}

%==========================================================================
end
