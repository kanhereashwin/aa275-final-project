%% parameters

% state dim
x_dim = 10; %(x, y, xdot, ydot, psi, b_ax, b_ay, b_psi_dot, cdt, cdt_dot)
m_dim = 6; % 6 satellites
u_dim = 3;

% timing
t_sample = 0.1 ;
t_total = 6 ; % seconds

% get W1 -- noise in IMU measurements
accel_meas_noise_var = 0.05;
ang_vel_meas_noise_var = 0.01;
W1 = diag([accel_meas_noise_var*ones(2,1); ang_vel_meas_noise_var]); 

% get W2 -- noise in IMU biases
accel_bias_noise_var = 0.01;
ang_vel_bias_noise_var = 0.002;
W2 = diag([accel_bias_noise_var*ones(2,1); ang_vel_bias_noise_var]); 

% get combined IMU noise matrix
W = [W1, zeros(u_dim); zeros(u_dim), W2];

% get Q for timing
Q_true_time = [0.0114, 0.0019; 0.0019, 0.0039];
Q_filt_time = [0.07, 0.04; 0.04, 0.08];

% get Q and R matrix
%Q = (1/10000).*make_random_covariance_matrix(x_dim) ;  % represents process noise in state from IMU noise
R_scale = 25; % corresponds to variance from each GPS pseudorange
%R = 5*eye(m_dim) ;  % represents noise from GPS measurements

%% Define satellite locations
prns = 1:6;
sat1 = [  100,   100, 100];
sat2 = [ -100,  -100, 100];
sat3 = [ -100,   100, 100];
sat4 = [  100,  -100, 100];
sat5 = [   50,     0, 100];
sat6 = [  -50,     0, 100];
sats_pos = [sat1; sat2; sat3; sat4; sat5; sat6];


%% run estimator
% set up time
t_vec = 0:t_sample:t_total ;
n_t = length(t_vec) ;

% initial state and covariance estimate ("new" corresponds to time step k,
% and "old" to time step k - 1) ;
init_pos_var = 1;  init_vel_var = 0.1;  init_psi_var = 0.01;
init_clk_var = 1;  init_clk_drift_var = 0.1;
x_init_est = [init_pos_var*randn(2,1) ; init_vel_var*randn(2,1); ...
    init_psi_var*randn(1,1); zeros(5,1)] ; %start 0 biases / clock drift
P_new = 2*diag([init_pos_var*ones(2,1); init_vel_var*ones(2,1); ...
    init_psi_var; ones(3,1); init_clk_var; init_clk_drift_var])  ;

% set up arrays to store true state
init_bias = [ 0.01 * rand(2,1); 0.01*rand(1,1) ]; % true bias (initial)
init_clk_states = mvnrnd(zeros(2,1), Q_true_time)';
x_init_true = [ zeros(5,1); init_bias; init_clk_states ];
x_new_true = x_init_true;
X_true = [x_init_true, nan(x_dim,n_t-1)] ;

% set up to store measurements and state estimate
exp_meas = expected_meas( x_init_true, prns, sats_pos );
init_meas = mvnrnd( exp_meas, R_scale*eye(length(prns)))' ;
X_meas = [init_meas, nan(m_dim,n_t-1)] ;
x_new_est = x_init_est;
X_est = [x_init_est, nan(x_dim, n_t-1)]; 

% Create "true" vehicle motion
% S_true = zeros(2,n_t-1);
% Omega_true = zeros(1, n_t-1);
[S_true, Omega_true] = straight_line_path();

%% Start simulation
for idx = 2:n_t
    % take a step with the true (noisy) dynamics
    x_old_true = X_true(:,idx-1) ;
    s_curr = S_true(:,idx-1) ;
    omega_curr = Omega_true(:,idx-1);
    u_k = [s_curr; omega_curr]; % "true" accel/ang. vel
    bias_k = x_old_true(6:8);
    time_states_k = x_old_true(9:10);
    u_k_bias = u_k - bias_k; % biased motion meas ("true" IMU meas)
    [A_true,B_true,~] = chimera_dyn(x_old_true, t_sample);
    x_new_true = A_true*x_old_true + B_true*u_k_bias ;
    
    % Update IMU bias
    imu_bias_noise = mvnrnd(zeros(u_dim,1), W2)';
    bias_new = bias_k + imu_bias_noise;
    x_new_true(6:8) = bias_new;
    
    % Update clock bias / drift
    time_states_noise = mvnrnd(zeros(2,1), Q_true_time)';
    time_states_new = time_states_k + time_states_noise;
    x_new_true(9:10) = time_states_new;
    
    % Model IMU measurement
    imu_meas_noise = mvnrnd(zeros(u_dim,1), W1)';
    u_k_meas = u_k_bias + imu_meas_noise; % received IMU meas
        
    % get measurement from noisy state
    R = R_scale*eye( length(prns) );
    meas_noise = mvnrnd(zeros(length(prns),1), R)' ;
    true_meas = expected_meas(x_new_true, prns, sats_pos);
    z_new = true_meas + meas_noise ; % new GPS measurement
    
    % PREDICT
    x_old_est = X_est(:,idx-1) ;
    [A_est,B_est,G_est] = chimera_dyn(x_old_est, t_sample);
    Q_est = G_est * W * G_est'; % estimated Q matrix
    Q_est(9:10, 9:10) = Q_filt_time;
    x_new_est = A_est*x_old_est + B_est*u_k_meas; % predicted state est
    P_new = A_est * P_new * A_est' + Q_est;

    % UPDATE
    y_new = z_new - expected_meas(x_new_est, prns, sats_pos) ; % innovation
    H = chimera_meas(x_new_est, prns, sats_pos);
    S_new = H * P_new * H' + R ; % innovation covariance
    K_new = P_new * H' * inv(S_new) ; % near-optimal Kalman gain
    x_new_est = x_new_est + K_new * y_new ; % updated state estimate
    P_new = (eye(x_dim) - K_new*H)*P_new ; % updated covariance estimate
    
    % save data
    X_true(:,idx) = x_new_true ;
    X_meas(:,idx) = z_new ;
    X_est(:,idx) = x_new_est ;
end

%% compute errors
%E_meas = vecnorm(X_true(1:2, :) - X_meas) ; % measurement error
E_est = vecnorm(X_true(1:2, :) - X_est(1:2, :)) ; % position estimate error

%% plotting
figure(1) ; clf

subplot(2,1,1) ; axis equal ; hold on ; grid on ;

plot_path(X_true(1:2,:),'g-')
%plot_path(X_meas(1:2,:),'r-')
plot_path(X_est(1:2,:),'b-')
legend('true','est')
xlabel('x_1')
ylabel('x_2')
title('trajectories in 2D space')


subplot(2,1,2) ; hold on ;
%plot(t_vec,E_meas,'r')
plot(t_vec,E_est,'b')
legend('est')
xlabel('time')
ylabel('error (2-norm)')


%% helper functions
function [A,B,G] = chimera_dyn(x,ts)
    A = eye(10);
    A(1:2, 3:4) = ts*eye(2);
    psi = x(5);
    rot_mat = [ cos(psi), sin(psi); -1*sin(psi), cos(psi) ];
    A(3:4, 6:7) = ts*rot_mat;
    A(5,8) = ts;
    A(9,10) = ts;
    
    B = zeros(10,3);
    B(3:4, 1:2) = ts*rot_mat;
    B(5,3) = ts;
    
    G = get_G(x, ts);

end

function G = get_G(x, ts)
    psi = x(5);
    rot_mat = [ cos(psi), sin(psi); -1*sin(psi), cos(psi) ];
    G = zeros(10,6);
    G(3:4, 1:2) = ts*rot_mat;
    G(5,3) = ts;
    G(6:8, 4:6) = eye(3);
end

function H = chimera_meas(x, prns, sats)
    m = length(prns);
    H = zeros(m, 10);
    sats_pos_prns = sats(prns, :);
    
    diff_pos = [x(1), x(2), 0] - sats_pos_prns; % states 1-2 are (x,y)
    cdt = x(9); 
    rho_vals = sqrt( sum(diff_pos .* diff_pos, 2) ) + cdt ;
    
    H(:, 1:2) = diff_pos(:, 1:2) ./ rho_vals;
    H(:, 9) = ones(m, 1);
end

function rho_vals = expected_meas(x, prns, sats)
    sats_pos_prns = sats(prns, :);
    diff_pos = [x(1), x(2), 0] - sats_pos_prns; % states 1-2 are (x,y)
    cdt = x(9); 
    rho_vals = sqrt( sum(diff_pos .* diff_pos, 2) ) + cdt ;
end
