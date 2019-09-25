function output = FFNAV_EKF_Propagation(time_step, mu, state_pre, Q)
% FFNAV Extended Kalman Filter ============================================
% Description: This function completes the state and covariance propagation
% (time update) of the EKF algorithm, and passes the a priori data to the
% correction phase
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
% Latest Edit: Cory Fraser - June 25, 2017
%% ========================================================================
% Previous State Estimates
    x            = state_pre(1);
    y            = state_pre(2);
    z            = state_pre(3);
    theta        = state_pre(4);
    rt           = state_pre(5);
    x_dot        = state_pre(6);
    y_dot        = state_pre(7);
    z_dot        = state_pre(8);
    theta_dot    = state_pre(9);
    rt_dot       = state_pre(10);
    
    rc  = sqrt((rt + x)^2 + y^2 + z^2);  
    
    state_priori        = zeros(10,1);
    state_priori =  [   x; y; z; theta; rt; ...
                        x_dot; y_dot; z_dot; ...
                        theta_dot; rt_dot ]; %(10 x 1)

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
          
%% ========================================================================
% Initialize Parameters
    T               = time_step;              %Time step for integration

%Manual Covariance 
%{    
%Covariance of the process noise (Depends on the system)
    sig2_r           = 2*10^-6; %Good - July 6 
    %sig2_r           = 2*10^2;  %Old value - June 25   
    %sig2_r           = 2*10^-2; %Better (0.07 EKF improvement)
    %sig2_r           = 2*10^-3; %Better (0.018 EKF improvement)
    %sig2_r           = 2*10^-5; 
    %sig2_r           = 2*10^-7; 
    
    
    sig2_theta       = 5*10^-10; %Good - July 6 
    %sig2_theta       = 5*10^-8; %Old value - June 25
    %sig2_theta       = 5*10^-5; %Worse, theta std increasing past 3
    
    
    sig2_rt          = 5*10^-10; %Good - July 6 
    %sig2_rt          = 5*10^-3; 
    
    sig2_v           = 5*10^-9; % Better, velocity error std < 0.03
    %sig2_v           = 5*10^-3; %Original - velocity error std < 0.8
    %sig2_v           = 5*10^-5; % Better, velocity error std < 0.25
    %sig2_v           = 5*10^-7; % Better, velocity error std < 0.08
    
    
    sig2_thetadot    = 5*10^-15;   %Good - July 6 
    %sig2_thetadot    = 5*10^-12;   % Better covariance convergence for rt and theta
    %sig2_thetadot    = 5*10^-10;   %Original - theta_dot std < 0.025
    %sig2_thetadot    = 5*10^-7;    % Poor covariance convergence for rt and theta
    %sig2_thetadot    = 5*10^-5;    % Poor covariance convergence for rt and theta
    
    sig2_rtdot       = 5*10^-10; %Original - rt_dot std < 0.01
    %sig2_rtdot       = 5*10^-5; %Larger SS covariance for rt and theta
    %sig2_rtdot       = 5*10^-1; %Larger SS covariance for rt and theta
    %sig2_rtdot       = 5*10^0; %Larger SS covariance for rt and theta
    %sig2_rtdot       = 5*10^1; %Larger SS covariance for rt and theta
      
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
 %}   
   
Q_k = Q;
% =========================================================================
%EKF Step 2: Propagation of the covariance matrix to get P_priori
% Propagation is completed using a discrete-time linear model of the system

%Partial derviations of the 5 non-linear equations to establish F

    %Derivatives of x-acceleration equation
    dxddot_dx           = (mu/rc^3)*(3*((rt+x)/rc)^2 - 1) + theta_dot^2;
    dxddot_dy           = 3*mu*((rt+x)/(rc^5))*y - ...
                          2*theta_dot*rt_dot/rt;
    dxddot_dz           = 3*mu*((rt+x)/rc^5)*z;
    dxddot_dtheta       = 0;
    dxddot_drt          = 2*theta_dot*rt_dot*y/rt^2 - ...
                          2*mu/rt^3 + (mu/rc^3)*(3*((rt+x)/rc)^2 - 1);
    dxddot_dx_dot       = 0; 
    dxddot_dy_dot       = 2*theta_dot;
    dxddot_dz_dot       = 0;
    dxddot_dtheta_dot   = 2*(theta_dot*x + y_dot - ...
                          (rt_dot/rt)*y);
    dxddot_drt_dot      = -2*theta_dot*y/rt;

    %Derivatives of y-acceleration equation
    dyddot_dx           = 2*theta_dot*rt_dot/rt + ...
                          3*mu*((rt+x)/rc^5)*y;
    dyddot_dy           = (mu/rc^3)*(3*(y/rc)^2-1) + theta_dot^2;
    dyddot_dz           = 3*mu*(y/rc^5)*z;
    dyddot_dtheta       = 0;
    dyddot_drt          = 3*mu*((rt+x)/rc^5)*y - ...
                          2*theta_dot*rt_dot*x/rt^2;
    dyddot_dx_dot       = -2*theta_dot; 
    dyddot_dy_dot       = 0;
    dyddot_dz_dot       = 0;
    dyddot_dtheta_dot   = 2*(theta_dot*y - x_dot + ...
                          rt_dot*x/rt);
    dyddot_drt_dot      = 2*theta_dot*x/rt;

    %Derivatives of z-acceleration equation
    dzddot_dx           = 3*mu*((rt+x)/rc^5)*z;
    dzddot_dy           = 3*mu*y*z/rc^5;
    dzddot_dz           = (mu/rc^3)*(3*(z/rc)^2-1);
    dzddot_dtheta       = 0;
    dzddot_drt          = 3*mu*((rt+x)/rc^5)*z;
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
    dtheta_ddot_drt         = 2*rt_dot*theta_dot/rt^2;
    dtheta_ddot_dx_dot      = 0; 
    dtheta_ddot_dy_dot      = 0;
    dtheta_ddot_dz_dot      = 0;
    dtheta_ddot_dtheta_dot  = -2*rt_dot/rt;
    dtheta_ddot_drt_dot     = -2*theta_dot/rt;
    
    %Derivatives of rt-acceleration equation
    drt_ddot_dx             = 0;
    drt_ddot_dy             = 0;
    drt_ddot_dz             = 0;
    drt_ddot_dtheta         = 0;
    drt_ddot_drt            = theta_dot^2 + 2*mu/rt^3;
    drt_ddot_dx_dot         = 0; 
    drt_ddot_dy_dot         = 0;
    drt_ddot_dz_dot         = 0;
    drt_ddot_dtheta_dot     = 2*theta_dot*rt;
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
% Converting data into output vector format
P_priori = [ P_priori(1,:) P_priori(2,:) P_priori(3,:) P_priori(4,:) P_priori(5,:) P_priori(6,:)...
           P_priori(7,:) P_priori(8,:) P_priori(9,:) P_priori(10,:) ]; %(1 x 100)
              
output = [ state_priori' P_priori]'; %(1x110)
end