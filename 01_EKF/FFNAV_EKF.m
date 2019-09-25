function output = FFNAV_EKF(sensors, time_step, mu, state_pre)
% FFNAV Extended Kalman Filter ============================================
% Description: This function completes the covariance proagation and state
% correction phases of the EKF algorithm.
%
% Inputs:
%   Initial Covariance Estimate
%      P_k

%   A priori state vector
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
% Latest Edit: Cory Fraser - April 17, 2017
%% ========================================================================
% Assign input vectors from the propagated states
    % States
    x_priori            = state_pre(1);
    y_priori            = state_pre(2);
    z_priori            = state_pre(3);
    theta_priori        = state_pre(4);
    rt_priori           = state_pre(5);
    x_dot_priori        = state_pre(6);
    y_dot_priori        = state_pre(7);
    z_dot_priori        = state_pre(8);
    theta_dot_priori    = state_pre(9);
    rt_dot_priori       = state_pre(10);
    
    rc  = sqrt((rt_priori + x_priori)^2 + y_priori^2 + z_priori^2);  
    
    state_priori = zeros(10,1);
    state_priori =  [   x_priori; y_priori; z_priori; theta_priori; rt_priori; ...
                        x_dot_priori; y_dot_priori; z_dot_priori; ...
                        theta_dot_priori; rt_dot_priori ]; %(10 x 1)

    % Unpropagated State Error Covariance (10x10)
    P_k = [ state_pre(11:20)'
            state_pre(21:30)'
            state_pre(31:40)'
            state_pre(41:50)'
            state_pre(51:60)'
            state_pre(61:70)'
            state_pre(71:80)'
            state_pre(81:90)'
            state_pre(91:100)'
            state_pre(101:110)' ];
       
    % Measured position and velocity parameters
    x_m             = sensors(1);
    y_m             = sensors(2);
    z_m             = sensors(3);
    theta_m         = sensors(4);
    x_dot_m         = sensors(5);
    y_dot_m         = sensors(6);
    z_dot_m         = sensors(7);
    

%% ========================================================================
% Initialize Parameters
    T               = time_step;              %Time step for integration
    
%{
    %Values as of June 12, 2017
%Covariance of the process noise (Depends on the system)
    sig2_r           = 2*10^1;  % 2*10^1, diverges below 1*10^1
    sig2_theta       = 5*10^-8; 
    sig2_rt          = 5*10^-8; 
    sig2_v           = 5*10^-3; 
    sig2_thetadot    = 5*10^-10;
    sig2_rtdot       = 5*10^-10;
   %}
    
    sig2_r           = 2*10^1; 
    sig2_theta       = 5*10^-8; 
    sig2_rt          = 5*10^-8; 
    sig2_v           = 5*10^-3; 
    sig2_thetadot    = 5*10^-10;
    sig2_rtdot       = 5*10^-10;
    Q_k             = [ sig2_r                           0 0 0 0 0 0 0 0 0 
                        0   sig2_r                         0 0 0 0 0 0 0 0 
                        0 0     sig2_r                       0 0 0 0 0 0 0 
                        0 0 0       sig2_theta                 0 0 0 0 0 0 
                        0 0 0 0         sig2_rt                  0 0 0 0 0 
                        0 0 0 0 0           sig2_v                 0 0 0 0 
                        0 0 0 0 0 0             sig2_v               0 0 0 
                        0 0 0 0 0 0 0               sig2_v             0 0 
                        0 0 0 0 0 0 0 0                 sig2_thetadot    0 
                        0 0 0 0 0 0 0 0 0                   sig2_rtdot    ];
    
%EKF Step 2: Propagation of the covariance matrix to get P_priori
% Propagation is completed using a discrete-time linear model of the system

%Partial derviations of the 5 non-linear equations to establish F

    %Derivatives of x-acceleration equation
    dxddot_dx           = (mu/rc^3)*(3*((rt_priori+x_priori)/rc)^2 - 1) + theta_dot_priori^2;
    dxddot_dy           = 3*mu*((rt_priori+x_priori)/(rc^5))*y_priori - ...
                          2*theta_dot_priori*rt_dot_priori/rt_priori;
    dxddot_dz           = 3*mu*((rt_priori+x_priori)/rc^5)*z_priori;
    dxddot_dtheta       = 0;
    dxddot_drt          = 2*theta_dot_priori*rt_dot_priori*y_priori/rt_priori^2 - ...
                          2*mu/rt_priori^3 + (mu/rc^3)*(3*((rt_priori+x_priori)/rc)^2 - 1);
    dxddot_dx_dot       = 0; 
    dxddot_dy_dot       = 2*theta_dot_priori;
    dxddot_dz_dot       = 0;
    dxddot_dtheta_dot   = 2*(theta_dot_priori*x_priori + y_dot_priori - ...
                          (rt_dot_priori/rt_priori)*y_priori);
    dxddot_drt_dot      = -2*theta_dot_priori*y_priori/rt_priori;

    %Derivatives of y-acceleration equation
    dyddot_dx           = 2*theta_dot_priori*rt_dot_priori/rt_priori + ...
                          3*mu*((rt_priori+x_priori)/rc^5)*y_priori;
    dyddot_dy           = (mu/rc^3)*(3*(y_priori/rc)^2-1) + theta_dot_priori^2;
    dyddot_dz           = 3*mu*(y_priori/rc^5)*z_priori;
    dyddot_dtheta       = 0;
    dyddot_drt          = 3*mu*((rt_priori+x_priori)/rc^5)*y_priori - ...
                          2*theta_dot_priori*rt_dot_priori*x_priori/rt_priori^2;
    dyddot_dx_dot       = -2*theta_dot_priori; 
    dyddot_dy_dot       = 0;
    dyddot_dz_dot       = 0;
    dyddot_dtheta_dot   = 2*(theta_dot_priori*y_priori - x_dot_priori + ...
                          rt_dot_priori*x_priori/rt_priori);
    dyddot_drt_dot      = 2*theta_dot_priori*x_priori/rt_priori;

    %Derivatives of z-acceleration equation
    dzddot_dx           = 3*mu*((rt_priori+x_priori)/rc^5)*z_priori;
    dzddot_dy           = 3*mu*y_priori*z_priori/rc^5;
    dzddot_dz           = (mu/rc^3)*(3*(z_priori/rc)^2-1);
    dzddot_dtheta       = 0;
    dzddot_drt          = 3*mu*((rt_priori+x_priori)/rc^5)*z_priori;
    dzddot_dx_dot       = 0; 
    dzddot_dy_dot       = 0;
    dzddot_dz_dot       = 0;
    dzddot_dtheta_dot   = 0;
    dzddot_drt_dot      = 0;

    %Derivatives of theta-acceleration equation
    dtheta_ddot_dx          = 0;
    dtheta_ddot_dy          = 0;
    dtheta_ddot_dz          = 0;
    dtheta_ddot_dtheta      = 0;
    dtheta_ddot_drt         = 2*rt_dot_priori*theta_dot_priori/rt_priori^2;
    dtheta_ddot_dx_dot      = 0; 
    dtheta_ddot_dy_dot      = 0;
    dtheta_ddot_dz_dot      = 0;
    dtheta_ddot_dtheta_dot  = -2*rt_dot_priori/rt_priori;
    dtheta_ddot_drt_dot     = -2*theta_dot_priori/rt_priori;
    
    %Derivatives of rt-acceleration equation
    drt_ddot_dx             = 0;
    drt_ddot_dy             = 0;
    drt_ddot_dz             = 0;
    drt_ddot_dtheta         = 0;
    drt_ddot_drt            = theta_dot_priori^2 + 2*mu/rt_priori^3;
    drt_ddot_dx_dot         = 0; 
    drt_ddot_dy_dot         = 0;
    drt_ddot_dz_dot         = 0;
    drt_ddot_dtheta_dot     = 2*theta_dot_priori*rt_priori;
    drt_ddot_drt_dot        = 0;

    %Defining the state matrix Ak
    F11 = zeros(5,5);
    F12 = eye(5,5);

    F21 = [dxddot_dx dxddot_dy dxddot_dz dxddot_dtheta dxddot_drt]; 
    F22 = [dxddot_dx_dot dxddot_dy_dot dxddot_dz_dot dxddot_dtheta_dot dxddot_drt_dot];

    F31 = [dyddot_dx dyddot_dy dyddot_dz dyddot_dtheta dyddot_drt];
    F32 = [dyddot_dx_dot dyddot_dy_dot dyddot_dz_dot dyddot_dtheta_dot dyddot_drt_dot];

    F41 = [dzddot_dx dzddot_dy dzddot_dz dzddot_dtheta dzddot_drt];
    F42 = [dzddot_dx_dot dzddot_dy_dot dzddot_dz_dot dzddot_dtheta_dot dzddot_drt_dot];

    F51 = [dtheta_ddot_dx dtheta_ddot_dy dtheta_ddot_dz dtheta_ddot_dtheta dtheta_ddot_drt];
    F52 = [dtheta_ddot_dx_dot dtheta_ddot_dy_dot dtheta_ddot_dz_dot ...
           dtheta_ddot_dtheta_dot dtheta_ddot_drt_dot];
    
    F61 = [drt_ddot_dx drt_ddot_dy drt_ddot_dz drt_ddot_dtheta drt_ddot_drt];
    F62 = [drt_ddot_dx_dot drt_ddot_dy_dot drt_ddot_dz_dot drt_ddot_dtheta_dot ...
           drt_ddot_drt_dot];

%Assembling the discrete state transition matrix Fk (10 x 10)
    F_k = [ F11 F12
            F21 F22
            F31 F32
            F41 F42
            F51 F52
            F61 F62 ];
   
    %F_k = expm(A_k*T); %Use the below Numerical Approximation;
    phi_k = (F_k*T)*(((F_k*T)/2)*(((F_k*T)/3)*((F_k*T)/4 + eye(10,10))+eye(10,10))...
          +eye(10,10))+eye(10,10);
         
%Dynamics Propagation using RK Method        
    k1 = FFNAV_EKFdot([state_priori]);
    k2 = FFNAV_EKFdot([state_priori + (1/3)*k1*T]);
    k3 = FFNAV_EKFdot([state_priori + (1/6)*k1*T + (1/6)*k2*T]); 
    k4 = FFNAV_EKFdot([state_priori + (1/8)*k1*T + (3/8)*k3*T]);
    k5 = FFNAV_EKFdot([state_priori + (1/2)*k1*T - (3/2)*k3*T + 2*k4*T]);
 
    state_priori = state_priori + (1/6)*(k1 + 4*k4 + k5)*T; %(10 x 1)   

%Calculating the a priori state error covariance (10 x 10)
    P_priori = phi_k*P_k*phi_k' + Q_k;
%% ========================================================================
% EKF Correction Step: Correction of the a priori estimate using the measured data

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
            0 0 0 1 0 0 0 0 0 0
            0 0 0 0 0 1 0 0 0 0
            0 0 0 0 0 0 1 0 0 0
            0 0 0 0 0 0 0 1 0 0];
     
%Covariance of the measurement noise (depends on the sensors)
    sig2_rm           = 1*10^2;                     
    sig2_thetam       = 1*10^-5;
    sig2_vm           = 1*10^0; 
    
    R_k             = [ sig2_rm                           0 0 0 0 0 0 
                        0   sig2_rm                         0 0 0 0 0 
                        0 0     sig2_rm                       0 0 0 0 
                        0 0 0       sig2_thetam                 0 0 0 
                        0 0 0 0            sig2_vm                0 0 
                        0 0 0 0 0              sig2_vm             0   
                        0 0 0 0 0 0                sig2_vm             ];            
        
%Calculating the Kalman Gain (10 x 7 matrix)
%Standard Implementation 
K = (P_priori*H')*inv(H*P_priori*H'+R_k);

%Previous state vector is a 10 x 1 matrix

%Calulating the estimated outputs (7 x 1 matrix)
    Z_est = H*state_priori;

%Correcting the state 
    residuals   = (Z_m - Z_est);                %(7 x 1 matrix)
    correction  = K*residuals;                  %(10 x 1 matrix)
    state_post  = state_priori + correction;    %(10 x 1 matrix)
       
%% ========================================================================
% EKF Step 4: Correct the state error covariance matrix to get P_post (10 x 10 matrix)
P_post = P_priori - K*(H*P_priori*H'+R_k)*K';

%Alternate methods of calculating P
%P_post = (P_post+P_post')/2; %Method to symmetrize P_post - no accuracy benefit

%Joseph's stabilized, from Bierman - Works, comparable to original method
    %P_bar   = P_priori - K*(P_priori*H')';
    %P_post = P_bar - P_bar*H'*K' + K*R_k*K';
    
%Alternate Joseph's form from Bierman - Works, similar as above, used in Vigneron's
%P_post = (eye(10)-K*H)*P_priori*(eye(10)-K*H)' + K*R_k*K'; 

%Standard used in Bierman, not good
    %P_post = P_priori - K*(P_priori*H')'; 

% D'Amicos Method
%{
    alph    = P_priori*H';
    beta    = inv(H*alph + R_k);
    K       = alph*beta;

    alph_t  = H*P_priori;
    P_post  = P_priori - K*alph_t;
%}
    
%% ========================================================================
% Converting data into output vector format
P_post = [ P_post(1,:) P_post(2,:) P_post(3,:) P_post(4,:) P_post(5,:) P_post(6,:)...
           P_post(7,:) P_post(8,:) P_post(9,:) P_post(10,:) ]; %(1 x 100)
       
K = [ K(1,:) K(2,:) K(3,:) K(4,:) K(5,:) K(6,:) K(7,:) K(8,:) K(9,:) K(10,:)]; %(1 x 70)
       
output = [ state_post' P_post Z_est' K residuals']'; %(1x194)
end