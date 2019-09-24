%% FFNAV_Prop ====================================================
% Description: This script propogates the NERMs using an RK5
% integration technqique
%
% Created by: Cory Fraser - Summer 2016
% Latest Edit: Cory Fraser - Winter 2016 
%% ========================================================================
%Initialization
fprintf('\n--------------------------------------------------------------\n')
fprintf('         NERM RELATIVE DYNAMICS PROPAGATION          \n')
fprintf('--------------------------------------------------------------\n')

filein = 'FFNAV_LargePEO_Params.mat';
load(filein)
fprintf('\nParameter file %s has been loaded...\n', filein)

 %% ======================================================================
% EKF Simulation Setup
T       = 1;        % Time step for the simulation
time_end        = 30000;    % Duration of the simulation
mu              = 398600.4418;          %Earth's Gravitational Constant (km^3/s^2)
time =  [0:T:time_end]';

% Assign input vectors from the propagated states
    % States
   
    x_priori            = state0(1);
    y_priori            = state0(2);
    z_priori            = state0(3);
    theta_priori        = state0(4);
    rt_priori           = state0(5);
    x_dot_priori        = state0(6);
    y_dot_priori        = state0(7);
    z_dot_priori        = state0(8);
    theta_dot_priori    = state0(9);
    rt_dot_priori       = state0(10);
   
    state        = zeros(10,1);
    state        = [        x_priori 
                            y_priori 
                            z_priori
                            theta_priori
                            rt_priori
                            x_dot_priori
                            y_dot_priori
                            z_dot_priori
                            theta_dot_priori
                            rt_dot_priori    ]; 

tic
for i = 1:T:time_end
   
    k1 = FFNAV_EKFdot([state(:,i)]);
    k2 = FFNAV_EKFdot([state(:,i) + (1/3)*k1*T]);
    k3 = FFNAV_EKFdot([state(:,i) + (1/6)*k1*T + (1/6)*k2*T]); 
    k4 = FFNAV_EKFdot([state(:,i) + (1/8)*k1*T + (3/8)*k3*T]);
    k5 = FFNAV_EKFdot([state(:,i) + (1/2)*k1*T - (3/2)*k3*T + 2*k4*T]);
 
    state(:,i+1) = state(:,i) + (1/6)*(k1 + 4*k4 + k5)*T; %(10 x 1)
end
toc


fprintf('Simulation complete.\n')
fprintf('--------------------------------------------------------------\n')
TF_size       = 20;               %Title font size
    AF_size       = 16;               %Axes label font size
    AT_size       = 12;               %Tick label font size

    plot3(state(1,:)*10^3, state(2,:)*10^3, state(3,:)*10^3);
    set(gca,'FontSize',20);
    title('Propagated Nonlinear Equations of Motion','FontSize', 20)
    xlabel('X Position (m)', 'FontSize', 20)
    ylabel('Y Position (m)', 'FontSize', 20)
    zlabel('Z Position (m)', 'FontSize', 20)
    set(gca,'FontName', 'Times')
    grid on
