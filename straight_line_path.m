% Create a straight line path
function [S_true, Omega_true] = straight_line_path()
% Define path from start / end points
x_start = 0;  y_start = 0;
x_end = 60;  y_end = 60;

% Define (x,y) trajectory
t_sample = 0.1 ;
t_total = 6 ;
n_t = (t_total / t_sample) + 1 ;
x_traj = linspace(x_start, x_end, n_t-1);
y_traj = linspace(y_start, y_end, n_t-1);
x_traj = [x_start, x_traj]; % stay at initial point 1 time step, then go
y_traj = [y_start, y_traj];

% Obtain yaw angles
diff_x = diff(x_traj);
diff_y = diff(y_traj);
psi_traj = atan2(diff_y, diff_x);
psi_traj = [psi_traj, psi_traj(end)];

% Obtain velocities
x_vel_traj = (1/t_sample)*diff_x;
y_vel_traj = (1/t_sample)*diff_y;
x_vel_traj = [x_vel_traj, x_vel_traj(end)];
y_vel_traj = [y_vel_traj, y_vel_traj(end)];

% Obtain accelerations
diff_x_vel = diff(x_vel_traj);
diff_y_vel = diff(y_vel_traj);
x_accel_traj = nan*ones(1, n_t-1); % (only need n_t - 1)
y_accel_traj = nan*ones(1, n_t-1);
for i = 1:(n_t-1)  
    psi_i = psi_traj(i);
    rot_mat_I2B = [ cos(psi_i), -1*sin(psi_i); ...
                    sin(psi_i), cos(psi_i) ];
    diff_v = [diff_x_vel(i); diff_y_vel(i)];
    accel_i = (1/t_sample) * rot_mat_I2B * diff_v;
    x_accel_traj(i) = accel_i(1);
    y_accel_traj(i) = accel_i(2);
end

% True acceleration
S_true = [ x_accel_traj ; y_accel_traj ];

% Obtain true angular velocities
Omega_true = (1/t_sample) * diff(psi_traj); % (only need n_t - 1)
end

