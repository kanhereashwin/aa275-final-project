% Create a straight line path
function [S_true, Omega_true] = one_turn_path()
% Define path from start / end points
x_start = 0;  y_start = 0;
rad = 5; 
x_end = 60;  y_end = 60;
x_mid = x_start + rad; y_mid = y_end - rad;

% Define (x,y) trajectory
t_sample = 0.1 ;
t_total = 6 ;
n_t = (t_total / t_sample) + 1 ;
n_turn = 9;
n_line = (n_t - n_turn) / 2; 

% first leg
x_traj1 = x_start * ones(1, n_line - 1); % minus 1 to account for 1st pt
y_traj1 = linspace( y_start, y_mid, n_line - 1 ); 

% turn
theta_per_step = (pi / 2) / (n_turn + 1);
theta_vect = theta_per_step:theta_per_step:( (pi/2) - theta_per_step );
xp_turn = rad*cos(theta_vect);
yp_turn = rad*sin(theta_vect);
x_turn = x_start + rad - xp_turn;
y_turn = y_mid + yp_turn;

% second leg
x_traj2 = linspace( x_mid, x_end, n_line );
y_traj2 = y_end * ones(1, n_line);

% concatenate
x_traj = [x_start, x_traj1, x_turn, x_traj2]; % stay at initial point 1 time step, then go
y_traj = [y_start, y_traj1, y_turn, y_traj2];


%% Fix rest of path
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

