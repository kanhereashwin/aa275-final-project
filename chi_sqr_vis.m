x = 0:0.2:15;
y = chi2pdf(x,2);
thresh = 5.99; % 95-percentile statistic
figure;
plot(x,y);
xline(thresh, 'm--');
legend({'$\chi^2$ pdf','$95$-percentile threshold'}, 'interpreter', 'latex');
grid on;
xlabel('detection statistic (under null hypothesis $H_0$)', 'interpreter', 'latex');
ylabel('probability density ($2$ DOF)', 'interpreter', 'latex');
title('$\chi^2$ detection statistic pdf and threshold', 'interpreter', 'latex');