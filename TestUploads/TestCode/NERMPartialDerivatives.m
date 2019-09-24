syms x y z theta rt x_dot y_dot z_dot theta_dot rt_dot mu

x_ddot = @(x, y, z, theta, rt, x_dot, y_dot, z_dot, theta_dot, rt_dot) theta_dot^2*x + 2*theta_dot*y_dot - 2*theta_dot*rt_dot*y/rt + mu/rt^2 - mu*(rt+x)*((rt+x)^2+y^2+z^2)^(-3/2);
y_ddot = @(x, y, z, theta, rt, x_dot, y_dot, z_dot, theta_dot, rt_dot) theta_dot^2*y - 2*theta_dot*x_dot + 2*theta_dot*rt_dot*x/rt - mu*y*((rt+x)^2+y^2+z^2)^(-3/2);
z_ddot = @(x, y, z, theta, rt, x_dot, y_dot, z_dot, theta_dot, rt_dot) -mu*z*((rt+x)^2+y^2+z^2)^(-3/2);
theta_ddot = @(x, y, z, theta, rt, x_dot, y_dot, z_dot, theta_dot, rt_dot) -2*rt_dot*theta_dot/rt;
rt_ddot = @(x, y, z, theta, rt, x_dot, y_dot, z_dot, theta_dot, rt_dot) theta_dot^2*rt - mu/rt^2;

part_ddot_x = jacobian(x_ddot, [x, y, z, theta, rt, x_dot, y_dot, z_dot, theta_dot, rt_dot]).';
part_ddot_y = jacobian(y_ddot, [x, y, z, theta, rt, x_dot, y_dot, z_dot, theta_dot, rt_dot]).';
part_ddot_z = jacobian(z_ddot, [x, y, z, theta, rt, x_dot, y_dot, z_dot, theta_dot, rt_dot]).';
part_ddot_theta = jacobian(theta_ddot, [x, y, z, theta, rt, x_dot, y_dot, z_dot, theta_dot, rt_dot]).'
part_ddot_rt = jacobian(rt_ddot, [x, y, z, theta, rt, x_dot, y_dot, z_dot, theta_dot, rt_dot]).'

%{
%Using direct equations from the symbolic editor
    %Derivatives of x-acceleration equation
    dxddot_dx           = theta_dot_priori^2 - mu/((rt_priori + x_priori)^2 + y_priori^2 + z_priori^2)^(3/2) + (3*mu*(2*rt_priori + 2*x_priori)*(rt_priori + x_priori))/(2*((rt_priori + x_priori)^2 + y_priori^2 + z_priori^2)^(5/2));
    dxddot_dy           = (3*mu*y_priori*(rt_priori + x_priori))/((rt_priori + x_priori)^2 + y_priori^2 + z_priori^2)^(5/2) - (2*rt_dot_priori*theta_dot_priori)/rt_priori;
    dxddot_dz           = (3*mu*z_priori*(rt_priori + x_priori))/((rt_priori + x_priori)^2 + y_priori^2 + z_priori^2)^(5/2);
    dxddot_dtheta       = 0;
    dxddot_drt          = (3*mu*(2*rt_priori + 2*x_priori)*(rt_priori + x_priori))/(2*((rt_priori + x_priori)^2 + y_priori^2 + z_priori^2)^(5/2)) - mu/((rt_priori + x_priori)^2 + y_priori^2 + z_priori^2)^(3/2) - (2*mu)/rt_priori^3 + (2*rt_dot_priori*theta_dot_priori*y_priori)/rt_priori^2;
    dxddot_dx_dot       = 0; 
    dxddot_dy_dot       = 2*theta_dot_priori;
    dxddot_dz_dot       = 0;
    dxddot_dtheta_dot   = 2*y_dot_priori + 2*theta_dot_priori*x_priori - (2*rt_dot_priori*y_priori)/rt_priori;
    dxddot_drt_dot      = -2*theta_dot_priori*y_priori/rt_priori;

    %Derivatives of y-acceleration equation
    dyddot_dx           = (2*rt_dot_priori*theta_dot_priori)/rt_priori + (3*mu*y_priori*(2*rt_priori + 2*x_priori))/(2*((rt_priori + x_priori)^2 + y_priori^2 + z_priori^2)^(5/2));
    dyddot_dy           = theta_dot_priori^2 - mu/((rt_priori + x_priori)^2 + y_priori^2 + z_priori^2)^(3/2) + (3*mu*y_priori^2)/((rt_priori + x_priori)^2 + y_priori^2 + z_priori^2)^(5/2);
    dyddot_dz           = (3*mu*y_priori*z_priori)/((rt_priori + x_priori)^2 + y_priori^2 + z_priori^2)^(5/2);
    dyddot_dtheta       = 0;
    dyddot_drt          = (3*mu*y_priori*(2*rt_priori + 2*x_priori))/(2*((rt_priori + x_priori)^2 + y_priori^2 + z_priori^2)^(5/2)) - (2*rt_dot_priori*theta_dot_priori*x_priori)/rt_priori^2;
    dyddot_dx_dot       = -2*theta_dot_priori; 
    dyddot_dy_dot       = 0;
    dyddot_dz_dot       = 0;
    dyddot_dtheta_dot   = 2*theta_dot_priori*y_priori - 2*x_dot_priori + (2*rt_dot_priori*x_priori)/rt_priori;
    dyddot_drt_dot      = (2*theta_dot_priori*x_priori)/rt_priori;

    %Derivatives of z-acceleration equation
    dzddot_dx           = (3*mu*z_priori*(2*rt_priori + 2*x_priori))/(2*((rt_priori + x_priori)^2 + y_priori^2 + z_priori^2)^(5/2));
    dzddot_dy           = (3*mu*y_priori*z_priori)/((rt_priori + x_priori)^2 + y_priori^2 + z_priori^2)^(5/2);
    dzddot_dz           = (3*mu*z_priori^2)/((rt_priori + x_priori)^2 + y_priori^2 + z_priori^2)^(5/2) - mu/((rt_priori + x_priori)^2 + y_priori^2 + z_priori^2)^(3/2);
    dzddot_dtheta       = 0;
    dzddot_drt          = (3*mu*z_priori*(2*rt_priori + 2*x_priori))/(2*((rt_priori + x_priori)^2 + y_priori^2 + z_priori^2)^(5/2));
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
    dtheta_ddot_drt         = (2*rt_dot_priori*theta_dot_priori)/rt_priori^2;
    dtheta_ddot_dx_dot      = 0; 
    dtheta_ddot_dy_dot      = 0;
    dtheta_ddot_dz_dot      = 0;
    dtheta_ddot_dtheta_dot  = -(2*rt_dot_priori)/rt_priori;
    dtheta_ddot_drt_dot     = -(2*theta_dot_priori)/rt_priori;
    
    %Derivatives of rt-acceleration equation
    drt_ddot_dx             = 0;
    drt_ddot_dy             = 0;
    drt_ddot_dz             = 0;
    drt_ddot_dtheta         = 0;
    drt_ddot_drt            = (2*mu)/rt_priori^3 + theta_dot_priori^2;
    drt_ddot_dx_dot         = 0; 
    drt_ddot_dy_dot         = 0;
    drt_ddot_dz_dot         = 0;
    drt_ddot_dtheta_dot     = 2*rt_priori*theta_dot_priori;
    drt_ddot_drt_dot        = 0;
    %}
