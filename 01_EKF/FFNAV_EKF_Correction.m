function [output] = FFNAV_EKF_Correction(sensors, input_priori, R)
% FFNAV Extended Kalman Filter ============================================
% Description: This function completes the state and covariance correction 
% phase of the EKF algorithm, using measurments provided by GPS
%
% Inputs:
% A priori state vector
%    x_priori, y_priori, z_priori, theta_priori, rt_priori
%    x_dot_priori, y_dot_priori, z_dot_priori, theta_dot_priori, rt_dot_priori
%
% Measurements - Z vector
%       x_m, y_m, z_m, x_dot_m, y_dot_m, z_dot_m
%   
% Outputs:
%   a Posteriori State Estimates
%      x_post, y_post, z_post, theta_post, rt_post
%      x_dot_post, y_dot_post, z_dot_post, theta_dot_post, rt_dot_post
%
%   a Posteriori State Error Covariance Estimate
%      P_post
%
% Created by: Cory Fraser  - October 2015
% Latest Edit: Cory Fraser - July 03, 2017
%% ========================================================================
% Assign input vectors from the propagated states
    % States
    x_priori            = input_priori(1);
    y_priori            = input_priori(2);
    z_priori            = input_priori(3);
    theta_priori        = input_priori(4);
    rt_priori           = input_priori(5);
    x_dot_priori        = input_priori(6);
    y_dot_priori        = input_priori(7);
    z_dot_priori        = input_priori(8);
    theta_dot_priori    = input_priori(9);
    rt_dot_priori       = input_priori(10);
    
    state_priori        = zeros(10,1);
    state_priori =  [   x_priori; y_priori; z_priori; theta_priori; rt_priori; ...
                        x_dot_priori; y_dot_priori; z_dot_priori; ...
                        theta_dot_priori; rt_dot_priori ]; %(10 x 1)

% Unpropagated State Error Covariance (10x10)
    P_priori = [    input_priori(11:20)'
                    input_priori(21:30)'
                    input_priori(31:40)'
                    input_priori(41:50)'
                    input_priori(51:60)'
                    input_priori(61:70)'
                    input_priori(71:80)'
                    input_priori(81:90)'
                    input_priori(91:100)'
                    input_priori(101:110)' ];
          
%% ========================================================================
% EKF Correction Step: Correction of the a priori estimate using the measured data

% Measured position and velocity parameters
    x_m             = sensors(1);
    y_m             = sensors(2);
    z_m             = sensors(3);
    theta_m         = sensors(4);
    x_dot_m         = sensors(5);
    y_dot_m         = sensors(6);
    z_dot_m         = sensors(7);
    
%Assembling the measured outputs (8 x 1 matrix)
    Z_m = [     x_m
                y_m
                z_m
                theta_m
                x_dot_m
                y_dot_m
                z_dot_m];

%Defining the measurement model (7 x 10 matrix) %8x10
    H = [   1 0 0 0 0 0 0 0 0 0
            0 1 0 0 0 0 0 0 0 0
            0 0 1 0 0 0 0 0 0 0
            0 0 0 1 0 0 0 0 0 0 %Providing theta
            %0 0 0 0 1 0 0 0 0 0 %Providing rt
            0 0 0 0 0 1 0 0 0 0
            0 0 0 0 0 0 1 0 0 0
            0 0 0 0 0 0 0 1 0 0];

%Manual Covariance Entry
%{     
%Covariance of the measurement noise (depends on the sensors)
    sig2_rm           = 2*10^-4; %Good - July 7
    %sig2_rm           = 0.01; %(0.048)
    %sig2_rm           = 0.006; %(0.082)
    %sig2_rm           = 0.005; %(0.086)
    %sig2_rm           = 0.004; %(0.086)
    %sig2_rm           = 0.003; %(0.0819)
    %sig2_rm           = 2*10^-3; %Seems to smooth
    %sig2_rm           = 0.00005; %
    %sig2_rm           = 0.002; %From GPS error data stats (0.07)
    %sig2_rm           = 0.02; %(0.017)
    %sig2_rm           = 0.2; % (-0.84)
    %sig2_rm           = 2; % (-2.26, but positions look smoother)
    %sig2_rm           = 10; % (-3.85, but positions look smoother)
    %sig2_rm           = 20; % (Divergence)
    %sig2_rm           = 1*10^2;  %Old value - June 25 (0.04 EKF Avg, but now diverges)
    
    sig2_thetam       = 1*10^-9; %Good - July 7
    %sig2_thetam       = 1*10^5; %Old value - June 25
    %sig2_thetam       = 1*10^-7; %theta std stablizes at 1x10^-3
    
    sig2_vm           = 5*10^-9; %Good - July 6
    %sig2_vm           = 1*10^0; %Old value - June 25
    %sig2_vm           = 5*10^-5; %From GPS error data 
    %sig2_vm           = 5*10^-7; %Error bounds reduced
    %sig2_vm           = 5*10^-10; %Error bounds too low
    
    
    %Current Values
    sig2_rm           = 2*10^-4; %Good - July 7
    sig2_thetam       = 1*10^-9; %Good - July 7
    sig2_vm           = 5*10^-9; %Good - July 6
    
    R_k             = [ sig2_rm                           0 0 0 0 0 0 
                        0   sig2_rm                         0 0 0 0 0 
                        0 0     sig2_rm                       0 0 0 0 
                        0 0 0       sig2_thetam                 0 0 0 
                        0 0 0 0            sig2_vm                0 0 
                        0 0 0 0 0              sig2_vm             0   
                        0 0 0 0 0 0                sig2_vm             ];            
%}


%Calulating the Estimated Measurements (7x1 matrix)
    Z_est = H*state_priori;
  
%Calculating the residuals (7x1 matrix)
    residuals   = (Z_m - Z_est); 

%Calculating the theoretical Residual Covariance (7x7 Matrix)
    Pr_theo = H*P_priori*H'+R; 

%Calculating the Kalman Gain (10x7 matrix)    
    K = (P_priori*H')*inv(Pr_theo);
    %K = (P_priori*H')*inv(H*P_priori*H'+R_k);

%Correcting the State Estimate (10x1 matrix)
    correction  = K*residuals;                  
    state_post  = state_priori + correction;
    
%% ========================================================================
% EKF Step 4: Correct the state error covariance matrix to get P_post (10 x 10 matrix)

%Alternate Joseph's Form (Bierman, in Vigneron's)
    P_post = (eye(10)-K*H)*P_priori*(eye(10)-K*H)' + K*R*K'; 
    %ADD NOTE OF THIS IN WRITTEN EQUATIONS
%{
%Standards (Bierman, Aitken, not good)
    %P_post = P_priori - K*(P_priori*H')'; 
    %P_post = (eye(10)- K*H)*P_priori; 

%Generic (Ulrich, Sasidek, et al.)
    %P_post = P_priori - K*(H*P_priori*H'+R_k)*K';
    %P_post = P_priori - K*(Pr_theo)*K';

%Joseph's stabilized (Bierman) - Comparable to original method
    %P_bar   = P_priori - K*(P_priori*H')';
    %P_post = P_bar - P_bar*H'*K' + K*R_k*K';

% D'Amicos Method
%    alph    = P_priori*H';
%    beta    = inv(H*alph + R_k);
%    K       = alph*beta;
%    alph_t  = H*P_priori;
%    P_post  = P_priori - K*alph_t;

%Symmetrize P_post - no accuracy benefit
    %P_post = (P_post+P_post')/2; 
%}    
%% ========================================================================
% Converting data into output vector format
P_post = [ P_post(1,:) P_post(2,:) P_post(3,:) P_post(4,:) P_post(5,:) P_post(6,:)...
           P_post(7,:) P_post(8,:) P_post(9,:) P_post(10,:) ]; %(1 x 100)
       
Pr_theo = [ Pr_theo(1,:) Pr_theo(2,:) Pr_theo(3,:) Pr_theo(4,:) Pr_theo(5,:) Pr_theo(6,:)...
           Pr_theo(7,:)]; %(1 x 49)
       
K = [ K(1,:) K(2,:) K(3,:) K(4,:) K(5,:) K(6,:) K(7,:) K(8,:) K(9,:) K(10,:)]; %(1 x 70)
       
output = [ state_post' P_post Z_est' Pr_theo residuals' K]'; %(1x243)
end