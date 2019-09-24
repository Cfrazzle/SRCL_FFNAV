%% FFNAV EKF TRoubleshooting ==============================================
% Description: This 
%
% Function Calls:
%   COEtoRV     - Evaluates the state vector given orbital elements
%   COEfromRV   - Determines the orbital elements given the ECI state vector
%   
% Created by: Cory Fraser  - Apr.16, 2017
% Latest Edit: Cory Fraser - Apr.16, 2017
%% ========================================================================
clc; clear;
addpath('C:\Users\Cory\OneDrive\University\MECH5105 - Orbital Mechanics and Control\Orbital Mechanics Code')
global mu
mu              = 398600.4418;          %Earth's Gravitational Constant (km^3/s^2)

%Target Initial Conditions
COE(1).ecc      = 0;                   %True Anomaly (rad)
COE(1).SMA      = 7500;               %Semi-major Axis (km)
COE(1).inc      = 30;                  %Eccentricity
COE(1).RAAN     = 0;                 %Inclination (rad)
COE(1).AoP      = 0;                  %Right-ascension of the Ascending Node (rad)
COE(1).TA       = 0;                  %Argument of Perigee (rad)

%Chaser Initial Conditions
COE(2).ecc      = 0.00001; %True Anomaly (rad)
COE(2).SMA      = 7500;               %Semi-major Axis (km)
COE(2).inc      = 30;                  %Eccentricity
COE(2).RAAN     = 0;                 %Inclination (rad)
COE(2).AoP      = 0;                  %Right-ascension of the Ascending Node (rad)
COE(2).TA       = 0;                  %Argument of Perigee (rad)

%Converting COEs to State Vectors
[rt vt] = COEtoRV(COE(1), mu);
[rc vc] = COEtoRV(COE(2), mu);

%Calculating Relative Initial Conditions
r_rel = rc-rt;

omega_target =  cross(rt,vt)/norm(rt)^2;
v_rel = vc - vt - cross(omega_target,r_rel);

Lx      = rt/norm(rt);
Lz      = cross(rt, vt)/norm(cross(rt, vt));
Ly      = cross(Lz, Lx);

%Align unit vectors into rows
Lx      = Lx';
Ly      = Ly';
Lz      = Lz';

%Arranging vectrix for LVLH Frame
F_LVLH  = [ Lx      
            Ly
            Lz ];

%ECI Frame is defined as such
F_ECI   = [ 1 0 0 
            0 1 0
            0 0 1 ];
     
%Compute roration matrix - ECI to LVLH
C_LI    = F_LVLH*F_ECI';

%Apply rotation to obtain relative LVLH parameters
r_rel_LVLH = C_LI*r_rel; 
v_rel_LVLH = C_LI*v_rel;

fprintf('Target Spacecraft \n')
fprintf('    The initial position vector : R1 = [ %7.2f \t %7.2f \t %7.2f] km \n', rt)
fprintf('    The initial velocity vector : V1 = [ %7.2f \t %7.2f \t %7.2f] km/s \n', vt)
fprintf('    The initial orbital elements: \n')
fprintf('               Semi-major Axis:     a = %7.2f km  \n', COE(1).SMA)
fprintf('                  Eccentricity:     e = %7.2f     \n', COE(1).ecc)
fprintf('                  Inclincation:     i = %7.2f deg \n', COE(1).inc)
fprintf('                          RAAN:   Ohm = %7.2f deg \n', COE(1).RAAN)
fprintf('           Argument of Perigee: omega = %7.2f deg \n', COE(1).AoP)
fprintf('                  True Anomaly: theta = %7.2f deg \n', COE(1).TA)
fprintf('Chaser Spacecraft \n')
fprintf('    The final position vector : R2 = [ %7.2f \t %7.2f \t %7.2f] km \n', rc)
fprintf('    The final velocity vector : V2 = [ %7.2f \t %7.2f \t %7.2f] km/s \n', vc)
fprintf('    The final orbital elements: \n')
fprintf('               Semi-major Axis:     a = %7.2f km  \n', COE(2).SMA)
fprintf('                  Eccentricity:     e = %7.2f     \n', COE(2).ecc)
fprintf('                  Inclincation:     i = %7.2f deg \n', COE(2).inc)
fprintf('                          RAAN:   Ohm = %7.2f deg \n', COE(2).RAAN)
fprintf('           Argument of Perigee: omega = %7.2f deg \n', COE(2).AoP)
fprintf('                  True Anomaly: theta = %7.2f deg \n', COE(2).TA)
fprintf('Relative Motion - ECI Frame\n')
fprintf('    The final position vector : rho     = [ %12.4f \t %12.4f \t %12.4f] m \n', 1000*r_rel)
fprintf('    The final velocity vector : rho_dot = [ %12.4f \t %12.4f \t %12.4f] m/s \n', 1000*v_rel)
fprintf('Relative Motion - LVLH Frame\n')
fprintf('    The final position vector : rho     = [ %12.4f \t %12.4f \t %12.4f] m \n', 1000*r_rel_LVLH)
fprintf('    The final velocity vector : rho_dot = [ %12.4f \t %12.4f \t %12.4f] m/s \n', 1000*v_rel_LVLH)

fprintf('============================================================================== \n')


%Propogating the orbit
t0 = 0;
tf = 30000; 
time = t0:1:tf;
state_t = zeros(6, tf);
state_c = zeros(6, tf);

state_t0 = [ rt; vt];
state_c0 = [ rc; vc];

fprintf('Propagating the Two-Body Equations of Motion \n')
options = odeset('RelTol',1e-8,'AbsTol',1e-10);
[t_t,state_t]   = ode45(@TwoBody, time, state_t0, options);
[t_c,state_c]   = ode45(@TwoBody, time, state_c0, options);
fprintf('Simulation Completed \n\n')

rt_x  = state_t(:,1);
rt_y  = state_t(:,2);
rt_z  = state_t(:,3);
vt_x  = state_t(:,4);
vt_y  = state_t(:,5);
vt_z  = state_t(:,6);

rc_x  = state_c(:,1);
rc_y  = state_c(:,2);
rc_z  = state_c(:,3);
vc_x  = state_c(:,4);
vc_y  = state_c(:,5);
vc_z  = state_c(:,6);

del_x = rc_x-rt_x;
del_y = rc_y-rt_y;
del_z = rc_z-rt_z;
del_vx = vc_x-vt_x;
del_vy = vc_y-vt_y;
del_vz = vc_z-vt_z;


rt = [ rt_x rt_y rt_z ];
vt = [ vt_x vt_y vt_z ];
r_rel = [del_x del_y del_z];
v_rel = [del_vx del_vy del_vz];

for i = 1:length(v_rel)
    omega_target =  cross(rt(i,:),vt(i,:))/norm(rt(i,:))^2;
    v_rel(i,:) = v_rel(i,:) - cross(omega_target,r_rel(i,:));
end
figure(1)
hold on
plot(time,vt_x)
plot(time,vt_y)
plot(time,vt_z)
%plot3(rc_x, rc_y, rc_z)
grid on
xlabel(' X Position [km]')
ylabel(' Y Position [km]')
zlabel(' Z Position [km]')
legend('X', 'Y', 'Z')


figure(2)
hold on
plot3(del_x*1000, del_y*1000, del_z*1000)
xlabel(' Relative X Position [m]')
ylabel(' Relative Y Position [m]')
zlabel(' Relative Z Position [m]')
grid on

for i = 1:length(del_x)
    Lx = rt(i,:)/norm(rt(i,:));
    Lz = cross(rt(i,:), vt(i,:))/norm(cross(rt(i,:), vt(i,:)));
    Ly = cross(Lz, Lx);

    %Align unit vectors into rows
    Lx      = Lx';
    Ly      = Ly';
    Lz      = Lz';
    
    %Arranging vectrix for LVLH Frame
    F_LVLH  = [ Lx Ly Lz ]';
    
    %ECI Frame
    F_ECI   = [ 1 0 0 
                0 1 0
                0 0 1 ];
     
    %Compute roration matrix - ECI to LVLH
    C_LI    = F_LVLH*F_ECI';

    %Apply rotation to obtain relative LVLH parameters
    r_rel_LVLH(:,i) = C_LI*r_rel(i,:)'; 
    v_rel_LVLH(:,i) = C_LI*v_rel(i,:)';

end

r_rel_LVLH = r_rel_LVLH';
v_rel_LVLH = v_rel_LVLH';

figure(3)
hold on
plot(time,1000*r_rel_LVLH(:,1))
plot(time,1000*r_rel_LVLH(:,2))
plot(time,1000*r_rel_LVLH(:,3))
grid on
xlabel(' Time [s]')
ylabel(' Relative Position [m]')
legend('X', 'Y', 'Z')



figure(4)
hold on
plot(time, v_rel_LVLH(:,1)*1000)
plot(time, v_rel_LVLH(:,2)*1000)
plot(time, v_rel_LVLH(:,3)*1000)
xlabel(' Time [s]')
ylabel(' Relative Velocity [m/s]')
legend('X', 'Y', 'Z')
grid on
