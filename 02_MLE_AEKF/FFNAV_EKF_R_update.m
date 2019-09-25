function R_adapt = FFNAV_EKF_R_update(state_priori, Pv, time_step, mu)
% FFNAV Extended Kalman Filter ============================================
% Description: 
%
%
% Created by: Cory Fraser  - August 19, 2018
% Latest Edit: Cory Fraser - August 19, 2018
%% ========================================================================
%Maximum Likelihood Estimation Technique (MLE)

% Previous State Estimates
%{
    x_priori            = state_priori(1);
    y_priori            = state_priori(2);
    z_priori            = state_priori(3);
    theta_priori        = state_priori(4);
    rt_priori           = state_priori(5);
    x_dot_priori        = state_priori(6);
    y_dot_priori        = state_priori(7);
    z_dot_priori        = state_priori(8);
    theta_dot_priori    = state_priori(9);
    rt_dot_priori       = state_priori(10);
    
    rc  = sqrt((rt_priori + x_priori)^2 + y_priori^2 + z_priori^2);  
%}    
% A Priori State Error Covariance (10x10)
    P_priori = [ state_priori(11:20)'
            state_priori(21:30)'
            state_priori(31:40)'
            state_priori(41:50)'
            state_priori(51:60)'
            state_priori(61:70)'
            state_priori(71:80)'
            state_priori(81:90)'
            state_priori(91:100)'
            state_priori(101:110)' ];
%{          
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
   
    %F_k = expm(A_k*time_step); %Use the below Numerical Approximation;
    Phi = (F_k*time_step)*(((F_k*time_step)/2)*(((F_k*time_step)/3)*((F_k*time_step)/4 + eye(10,10))+eye(10,10))...
          +eye(10,10))+eye(10,10);
%% ===================================================
%Controller Updates
    sig2_rm           = 2*10^-4;
    sig2_thetam       = 1*10^-9;
    sig2_vm           = 5*10^-9;

R   = [ sig2_rm                           0 0 0 0 0 0 
        0   sig2_rm                         0 0 0 0 0 
        0 0     sig2_rm                       0 0 0 0 
        0 0 0       sig2_thetam                 0 0 0 
        0 0 0 0            sig2_vm                0 0 
        0 0 0 0 0              sig2_vm              0   
        0 0 0 0 0 0                sig2_vm             ];  
%}
H = [       1 0 0 0 0 0 0 0 0 0
            0 1 0 0 0 0 0 0 0 0
            0 0 1 0 0 0 0 0 0 0
            0 0 0 1 0 0 0 0 0 0
            0 0 0 0 0 1 0 0 0 0
            0 0 0 0 0 0 1 0 0 0
            0 0 0 0 0 0 0 1 0 0];


C_res = diag(Pv);

%Adaptations by Mohamed/Schwarz
R_adapt = C_res + H*P_priori*H';

% =========================================================================