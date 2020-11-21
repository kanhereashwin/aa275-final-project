%% parameters
% dynamics parameters
mu = 1.5 ;
x_max = 3.0 ;

% timing
t_sample = 0.01 ;
t_total = 15 ; % seconds

% noise (fixed for all time for now)
Q = (1/10000).*make_random_covariance_matrix(2) ; 
R = (1/100).*make_random_covariance_matrix(2) ;

%% automated from here
% set up symbolic state
x_sym = sym('x',[2 1]) ;

% create dynamics and Jacobian
f_sym = x_sym + t_sample.*van_der_pol_dyn(x_sym,mu) ;
F_sym = jacobian(f_sym) ;
f = matlabFunction(f_sym,'Vars',{x_sym}) ;
F = matlabFunction(F_sym,'Vars',{x_sym}) ;

% create measurement function and Jacobian
h_sym = van_der_pol_measurement(x_sym) ;
H_sym = jacobian(h_sym) ;
h = matlabFunction(h_sym,'Vars',{x_sym}) ;
H = matlabFunction(H_sym,'Vars',{x_sym}) ;

%% run estimator
% set up time
t_vec = 0:t_sample:t_total ;
n_t = length(t_vec) ;

% initial state and covariance estimate ("new" corresponds to time step k,
% and "old" to time step k - 1) ;
x_new = x_max.*rand(2,1) - x_max ./ 2 ;
P_new = eye(2) ;

% set up arrays to store true state
X_true = [x_new, nan(2,n_t-1)] ;

% set up to store measurements and state estimate
n_meas = mvnrnd(zeros(2,1),R)' ;
X_meas = [van_der_pol_measurement(x_new) + n_meas, nan(2,n_t-1)] ;
X_est = X_meas ;

for idx = 2:n_t
    % take a step with the true (noisy) dynamics
    x_old = X_true(:,idx-1) ;
    proc_noise = mvnrnd(zeros(2,1),Q)' ;
    x_new = f(x_old) + proc_noise ;
        
    % get measurement from noisy state
    meas_noise = mvnrnd(zeros(2,1),R)' ;
    z_new = h(x_new) + meas_noise ;
    
    % PREDICT
    x_hat_old = X_est(:,idx-1) ;
    x_hat_new = f(x_hat_old) ; % predicted state estimate
    F_new = F(x_hat_new) ;
    P_new = F_new * P_new * F_new' + Q ;
    
    % UPDATE
    y_new = z_new - h(x_hat_new) ; % innovation
    H_new = H(x_hat_new) ;
    S_new = H_new * P_new * H_new' + R ; % innovation covariance
    K_new = P_new * H_new' * inv(S_new) ; % near-optimal Kalman gain
    x_hat_new = x_hat_new + K_new * y_new ; % updated state estimate
    P_new = (eye(2) - K_new*H_new)*P_new ; % updated covariance estimate
    
    % save data
    X_true(:,idx) = x_new ;
    X_meas(:,idx) = z_new ;
    X_est(:,idx) = x_hat_new ;
end

%% compute errors
E_meas = vecnorm(X_true - X_meas) ; % measurement error
E_est = vecnorm(X_true - X_est) ; % state estimate error

%% plotting
figure(1) ; clf

subplot(2,1,1) ; axis equal ; hold on ; grid on ;

plot_path(X_true,'g-')
plot_path(X_meas,'r-')
plot_path(X_est,'b-')
legend('true','meas','est')
xlabel('x_1')
ylabel('x_2')
title('trajectories in phase space')


subplot(2,1,2) ; hold on ;
plot(t_vec,E_meas,'r')
plot(t_vec,E_est,'b')
legend('meas','est')
xlabel('time')
ylabel('error (2-norm)')

%% helper functions
function xd = van_der_pol_dyn(x,mu)
%     xd = [mu*(x(1,:) - (1/3).*x(1).^3 - x(2)) ;
%         (1/mu).*x(1)] ;
    xd = [x(2) ;
          mu*(1 - x(1)^2)*x(2) - x(1)] ;
end

function z = van_der_pol_measurement(x)
    z = x ;
end

function S = make_random_covariance_matrix(d)
    Q = rand(d);
    D = diag(rand(d,1)) ;
    S = Q*D*Q' ;
end