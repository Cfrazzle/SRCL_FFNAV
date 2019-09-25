function [output] = FFNAV_FAEKF_Correction(input_priori, sensors, R)
% FFNAV Extended Kalman Filter ============================================
% Description: This function completes the state and covariance correction
% (measurement update) of the EKF algorithm. The resulting a posteriori 
% data is used as the state estimate for the next time step.
%
% Inputs:
%   input_priori - Vector of a priori states and covariances
%   sensors      - Vector of measurements 
%   R            - Measurement noise covariance matrix
%   
% Outputs:
%   output      - Vector of a priori states and covariances 
%
% Created by: Cory Fraser - OCT 31, 2017
% Last Edit : Cory Fraser - JUL 04, 2018
% Copyright(c) 2018 by Cory Fraser
% =========================================================================

%% Initialize and assign data
    coder.extrinsic('chol');

% Propagated State Vector (10x1 Matrix)
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

state_priori =  [   x_priori; y_priori; z_priori; theta_priori; rt_priori; ...
                    x_dot_priori; y_dot_priori; z_dot_priori; ...
                    theta_dot_priori; rt_dot_priori ]; %(10 x 1)

% Propagated State Error Covariance (10x10 Matrix)
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

% Measurement Vector (7x1 Matrix)
x_m       = sensors(1);
y_m       = sensors(2);
z_m       = sensors(3);
theta_m   = sensors(4);
%rt_m   = sensors(4);
x_dot_m   = sensors(5);
y_dot_m   = sensors(6);
z_dot_m   = sensors(7);

Z_m = [     x_m
            y_m
            z_m
            theta_m
            %rt_m
            x_dot_m
            y_dot_m
            z_dot_m ];
            %rho_SDPR    ];

%% Defining the measurement model (7 x 10 matrix) 
%Relative position, velocity, and theta (7 measurements)
%
H = [   1 0 0 0 0 0 0 0 0 0
        0 1 0 0 0 0 0 0 0 0
        0 0 1 0 0 0 0 0 0 0
        0 0 0 1 0 0 0 0 0 0 
        0 0 0 0 0 1 0 0 0 0
        0 0 0 0 0 0 1 0 0 0
        0 0 0 0 0 0 0 1 0 0 ];
%}
%Relative position, velocity, and target radius (7 measurements)    
%{
    H = [   1 0 0 0 0 0 0 0 0 0
            0 1 0 0 0 0 0 0 0 0
            0 0 1 0 0 0 0 0 0 0
            0 0 0 0 1 0 0 0 0 0
            0 0 0 0 0 1 0 0 0 0
            0 0 0 0 0 0 1 0 0 0
            0 0 0 0 0 0 0 1 0 0 ];
%}

%% EKF Correction Step

% Calulating the estimated measurements (7x1 matrix)
Z_est = H*state_priori;

%Calculating the innovations (7x1 matrix)
innovations   = Z_m - Z_est;

%Calculating the theoretical covariance of the innovations (7x7 Matrix)
Pr_theo = H*P_priori*H'+R;

%Pr_inv  = inv(Pr_theo);        %This inversion gives numerical issues
Pr_inv = Pr_theo\eye(7);        %Seems sufficient for EKF, FAEKF

%Inverting the innovations covariance using Cholesky Decomposition
%{
    Pr_inv = zeros(7);
    PSDS_flag = 0;    %Flag for Positive Semi-Definite Symmetric
    [M_upper, PSDS_flag] = chol(Pr_theo);
       
    if PSDS_flag~=0
        fprintf('Pr_theo is not PDS \n')
    end
    M_inv = M_upper \ eye(7);
    Pr_inv = M_inv * M_inv';
%}

%Calculating the Kalman Gain (10x7 Matrix)
%K = (P_priori*H')*inv(H*P_priori*H'+R_k);
K = (P_priori*H')*Pr_inv;

%Correcting the state estimate (10x1 Matrix)
correction  = K*innovations;
state_post  = state_priori + correction;

%Correct the state error covariance matrix - Joseph's Form (10x10 Matrix)
P_post = (eye(10)-K*H)*P_priori*(eye(10)-K*H)' + K*R*K';

%% ========================================================================
% Converting data into output vector format
P_post = [ P_post(1,:) P_post(2,:) P_post(3,:) P_post(4,:) P_post(5,:) P_post(6,:)...
           P_post(7,:) P_post(8,:) P_post(9,:) P_post(10,:) ]; %(1x100)
       
Pr_theo = [ Pr_theo(1,:) Pr_theo(2,:) Pr_theo(3,:) Pr_theo(4,:) Pr_theo(5,:) Pr_theo(6,:)...
           Pr_theo(7,:)]; %(1x49)
       
K = [ K(1,:) K(2,:) K(3,:) K(4,:) K(5,:) K(6,:) K(7,:) K(8,:) K(9,:) K(10,:)]; %(1x70)
       
output = [ state_post' P_post Z_est' Pr_theo innovations' K]'; %(1x243)
end