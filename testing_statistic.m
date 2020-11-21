tot_failures = 0;
tot_attempts = 0;

fail_vect_tot = [];
n_traj = 1000;
for i = 1:n_traj
    dual_est_chimera_imu_ekf_tight_withFused;
    tot_failures = tot_failures + length(fail_idcs);
    tot_attempts = tot_attempts + (n_t-1);
    if i == 1
        fail_vect_tot = fail_vect;
    else
        fail_vect_tot = fail_vect_tot + fail_vect;
    end
end

disp(['Detection rate: ', num2str(tot_failures / tot_attempts)]);

%%
figure();

%% Plot how fail vector changes with time
plot(t_vec(2:end),fail_vect_tot / n_traj); hold on;
plot(t_vec(2:end),fail_vect_tot / n_traj, '*');
grid on;
xlabel('time (s)', 'interpreter', 'latex');
% ylabel('detection rate', 'interpreter', 'latex');
% title(['Detection rate with time (across ', num2str(n_traj), ' sampled trajectories)'], 'interpreter', 'latex');

ylabel('false alarm rate', 'interpreter', 'latex');
title(['False alarm rate with time (across ', num2str(n_traj), ' sampled trajectories)'], 'interpreter', 'latex');
%%
ylim([-0.05, 1.05]);
