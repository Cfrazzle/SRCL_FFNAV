function [output] = FFNAV_FAEKF_Propagation(input_pre, time_step, mu, Q)
% FFNAV Extended Kalman Filter Propagation ================================
% Description: This function completes the state and covariance propagation
% (time update) of the EKF algorithm. The resulting a priori data is
% then passed to the correction phase of the EKF algorithm.
%
% Inputs:
%   input_pre   - Vector of previous states and covariances
%   time_step   - Time step of the simulation
%   mu          - Earth's gravitational parameter
%   Q           - Process noise covariance matrix
%
% Outputs:
%   output      - Vector of a priori states and covariances 
%
% Other Functions Called:
%   FFNAV_EKFdot.m  - Calculates the state vector derivative components
%   FFNAV_EKF_STM.m - Calculates the state transition matrix
%
% Created by:  Cory Fraser - OCT 31, 2017
% Latest Edit: Cory Fraser - JUL 04, 2018
% Copyright(c) 2018 by Cory Fraser
% ========================================================================

%% Initialize and assign data
    T = time_step;              %Time step for integration 

% Previous State Vector (10x1 Matrix)
    x            = input_pre(1);
    y            = input_pre(2);
    z            = input_pre(3);
    theta        = input_pre(4);
    rt           = input_pre(5);
    x_dot        = input_pre(6);
    y_dot        = input_pre(7);
    z_dot        = input_pre(8);
    theta_dot    = input_pre(9);
    rt_dot       = input_pre(10);
       
    %Check theta for wrapping
    %if theta > 2*pi
    %    theta = AngleRangeRad(theta)
    %end
    
    state_pre =  [   x; y; z; theta; rt; ...
                        x_dot; y_dot; z_dot; ...
                        theta_dot; rt_dot ];

% Previous State Error Covariance (10x10 Matrix)
    P_pre = [ input_pre(11:20)'
            input_pre(21:30)'
            input_pre(31:40)'
            input_pre(41:50)'
            input_pre(51:60)'
            input_pre(61:70)'
            input_pre(71:80)'
            input_pre(81:90)'
            input_pre(91:100)'
            input_pre(101:110)' ];        
% =========================================================================
%% EKF Propagation Step

% Runge-Kutta-Merson integration to propagate state vector (10x1 Matrix)
    k1 = FFNAV_EKFdot([state_pre], mu);
    k2 = FFNAV_EKFdot([state_pre + (1/3)*k1*T], mu);
    k3 = FFNAV_EKFdot([state_pre + (1/6)*k1*T + (1/6)*k2*T], mu); 
    k4 = FFNAV_EKFdot([state_pre + (1/8)*k1*T + (3/8)*k3*T], mu);
    k5 = FFNAV_EKFdot([state_pre + (1/2)*k1*T - (3/2)*k3*T + 2*k4*T], mu);
 
    state_priori = state_pre + (1/6)*(k1 + 4*k4 + k5)*T;

% Discrete-time state transition matrix to propagate the covariance (10x10 Matrix)
    phi = FFNAV_EKF_STM(input_pre, T, mu);
    P_priori = phi*P_pre*phi' + Q;
% =========================================================================
%% Convert data into output vector format (1x110 Matrix)

P_priori = [ P_priori(1,:) P_priori(2,:) P_priori(3,:) P_priori(4,:)...
    P_priori(5,:) P_priori(6,:) P_priori(7,:) P_priori(8,:)...
    P_priori(9,:) P_priori(10,:) ];
              
output = [ state_priori' P_priori]'; 

end
% =========================================================================