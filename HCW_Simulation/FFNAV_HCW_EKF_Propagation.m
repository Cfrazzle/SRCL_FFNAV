function output = FFNAV_HCW_EKF_Propagation(time_step, state_pre, n, Q)
% FFNAV Extended Kalman Filter ============================================
% Description: This function completes the state and covariance propagation
% (time update) of the EKF algorithm, and passes the a priori data to the
% correction phase
%
% Inputs:
%   Initial Covariance Estimate
%      P_k

%   A priori state vector
%    x_priori, y_priori, z_priori
%    x_dot_priori, y_dot_priori, z_dot_priori
%
% Measurements - Z vector
%       x_m, y_m, z_m, x_dot_m, y_dot_m, z_dot_m
%   
% Outputs:
%   a Posteriori State Estimates
%      x_post, y_post, z_post
%      x_dot_post, y_dot_post, z_dot_post
%
%   a Posteriori State Error Covariance Estimate
%      P_post
%
% Created by: Cory Fraser  - October 2015
% Latest Edit: Cory Fraser - July 10, 2017
%% ========================================================================
% Previous State Estimates
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
    P_k                 = [ state_pre(7:12)'
                            state_pre(13:18)'
                            state_pre(19:24)'
                            state_pre(25:30)'
                            state_pre(31:36)'
                            state_pre(37:42)' ]; 
          
%% ========================================================================
% Initialize Parameters
    T               = time_step;              %Time step for integration
    
%Covariance of the process noise (Depends on the system)
%{
    sig2_r           = 2*10^-6; %Good - July 6   
    sig2_v           = 5*10^-9; % Better, velocity error std < 0.03
       
     Q_k             = [ sig2_r          0 0 0 0 0 
                        0   sig2_r         0 0 0 0
                        0 0     sig2_r       0 0 0
                        0 0 0       sig2_v     0 0 
                        0 0 0 0         sig2_v   0 
                        0 0 0 0 0           sig2_v ];
%}
Q_k = Q;

%EKF Step 2: Propagation of the covariance matrix to get P_priori
% Propagation is completed using a discrete-time linear model of the system

%Partial derviations of the 3 linear equations to establish F

    %Derivatives of x-acceleration equation
    dxddot_dx           = 3*n^2;
    dxddot_dy           = 0;
    dxddot_dz           = 0;
    dxddot_dx_dot       = 0; 
    dxddot_dy_dot       = 2*n;
    dxddot_dz_dot       = 0;
    
    %Derivatives of y-acceleration equation
    dyddot_dx           = 0;
    dyddot_dy           = 0;
    dyddot_dz           = 0;
    dyddot_dx_dot       = -2*n; 
    dyddot_dy_dot       = 0;
    dyddot_dz_dot       = 0;
    
    %Derivatives of z-acceleration equation
    dzddot_dx           = 0;
    dzddot_dy           = 0;
    dzddot_dz           = -n*n;
    dzddot_dx_dot       = 0; 
    dzddot_dy_dot       = 0;
    dzddot_dz_dot       = 0;
    
    %Defining the state matrix Fk (6x6)
    F11 = zeros(3,3);
    F12 = eye(3,3);

    F21 = [dxddot_dx dxddot_dy dxddot_dz]; 
    F22 = [dxddot_dx_dot dxddot_dy_dot dxddot_dz_dot];

    F31 = [dyddot_dx dyddot_dy dyddot_dz ];
    F32 = [dyddot_dx_dot dyddot_dy_dot dyddot_dz_dot];

    F41 = [dzddot_dx dzddot_dy dzddot_dz];
    F42 = [dzddot_dx_dot dzddot_dy_dot dzddot_dz_dot ];

%Assembling the discrete state transition matrix Fk (6 x 6)
    F_k = [ F11 F12
            F21 F22
            F31 F32
            F41 F42 ];
   
%F_k = expm(A_k*T); %Use the below Numerical Approximation;
    phi_k = (F_k*T)*(((F_k*T)/2)*(((F_k*T)/3)*((F_k*T)/4 + eye(6,6))+eye(6,6))...
          +eye(6,6))+eye(6,6);
   
%Dynamics Propagation using RK Method        
    k1 = FFNAV_HCW_EKFdot([state_priori], n);
    k2 = FFNAV_HCW_EKFdot([state_priori + (1/3)*k1*T], n);
    k3 = FFNAV_HCW_EKFdot([state_priori + (1/6)*k1*T + (1/6)*k2*T], n); 
    k4 = FFNAV_HCW_EKFdot([state_priori + (1/8)*k1*T + (3/8)*k3*T], n);
    k5 = FFNAV_HCW_EKFdot([state_priori + (1/2)*k1*T - (3/2)*k3*T + 2*k4*T], n);
 
    state_priori = state_priori + (1/6)*(k1 + 4*k4 + k5)*T; %(6 x 1)   

%Calculating the a priori state error covariance (6 x 6)
    P_priori = phi_k*P_k*phi_k' + Q_k;
    
%% ========================================================================
% Converting data into output vector format
P_priori = [ P_priori(1,:) P_priori(2,:) P_priori(3,:) P_priori(4,:) P_priori(5,:) P_priori(6,:)]; %(1 x 36)
              
output = [ state_priori' P_priori]'; %(1x42)
end