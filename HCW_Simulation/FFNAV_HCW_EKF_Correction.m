function output = FFNAV_HCW_EKF_Correction(sensors, state_pre, R)
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
% Latest Edit: Cory Fraser - July 10, 2017
%% ========================================================================
% Assign input vectors from the propagated states
    % States
    x_priori            = state_pre(1);
    y_priori            = state_pre(2);
    z_priori            = state_pre(3);
    x_dot_priori        = state_pre(4);
    y_dot_priori        = state_pre(5);
    z_dot_priori        = state_pre(6);
    
    state_priori        = zeros(6,1);
    state_priori        = [ x_priori 
                            y_priori 
                            z_priori
                            x_dot_priori
                            y_dot_priori
                            z_dot_priori ]; 

    % Unpropagated State Error Covariance (6x6)
    P_priori            = [ state_pre(7:12)'
                            state_pre(13:18)'
                            state_pre(19:24)'
                            state_pre(25:30)'
                            state_pre(31:36)'
                            state_pre(37:42)' ];
        
    % Measured position and velocity parameters
    x_m             = sensors(1);
    y_m             = sensors(2);
    z_m             = sensors(3);
    x_dot_m         = sensors(4);
    y_dot_m         = sensors(5);
    z_dot_m         = sensors(6);
    
%% ========================================================================
% EKF Correction Step: Correction of the a priori estimate using the measured data

%Assembling the measured outputs (6 x 1 matrix)
    Z_m = [     x_m
                y_m
                z_m
                x_dot_m
                y_dot_m
                z_dot_m];

%%Defining the measurement model (6 x 6 matrix) 
    H = [   1 0 0 0 0 0
            0 1 0 0 0 0
            0 0 1 0 0 0
            0 0 0 1 0 0
            0 0 0 0 1 0
            0 0 0 0 0 1  ];
     
%Covariance of the measurement noise (depends on the sensors)
%{
    %Current Values
    sig2_rm           = 2*10^-4; %Good - July 7
    sig2_vm           = 5*10^-9; %Good - July 6
    
    R_k             = [ sig2_rm                 0 0 0 0 0  
                        0   sig2_rm               0 0 0 0  
                        0 0     sig2_rm             0 0 0  
                        0 0 0       sig2_vm           0 0 
                        0 0 0 0         sig2_vm         0   
                        0 0 0 0 0           sig2_vm      ];            
%}

R_k = R;
%Calculating the Kalman Gain (6 x 6 matrix)
%Standard Implementation 
K = (P_priori*H')*inv(H*P_priori*H'+R_k);

%Previous state vector is a 6 x 1 matrix

%Calulating the estimated outputs (6 x 1 matrix)
    Z_est = H*state_priori;

%Correcting the state 
    residuals   = (Z_m - Z_est);                %(6 x 1 matrix)
    correction  = K*residuals;                  %(6 x 1 matrix)
    state_post  = state_priori + correction;    %(6 x 1 matrix)
       
%% ========================================================================
% EKF Step 4: Correct the state error covariance matrix to get P_post (6 x 6 matrix)
P_post = P_priori - K*(H*P_priori*H'+R_k)*K';

%% ========================================================================
% Converting data into output vector format
P_post  = [ P_post(1,:) P_post(2,:) P_post(3,:) P_post(4,:) P_post(5,:) P_post(6,:)]; %(1 x 36)
K       = [ K(1,:) K(2,:) K(3,:) K(4,:) K(5,:) K(6,:)]; %(1 x 36)
output  = [ state_post' P_post Z_est' K residuals']'; %(1x90)


end