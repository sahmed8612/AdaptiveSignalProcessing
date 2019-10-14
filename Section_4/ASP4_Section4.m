%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% Section 4
%-----------------------------------------------------------------------
%
%-----------------------------------------------------------------------
% Question 1
%-----------------------------------------------------------------------
data = load('time-series.mat');
data = data.y;

% figure;
% yyaxis left; plot(data);
% ylabel('y[n]','fontsize',20);
% yyaxis right; plot([1, length(data)], [mean(data), mean(data)], 'linewidth', 2);
% ylabel('E\{y[n]\}','fontsize',20)
% xlabel('Sample, n','fontsize',20)
% set(gca,'fontsize',20);
% title('Timeseries','FontSize',20,'FontWeight','bold', 'HorizontalAlignment', 'center')
% set(gca,'fontsize',20)
% grid on
% %set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
% %print('../Images/Section_4/Section_4_1_1_M', '-depsc');
% 
% %zero-mean
% data_zero_mean = data-mean(data);
% model_order = 4;
% alpha = 1E-5;
% X = timeseries_to_matrix(data_zero_mean, model_order);
% y = data_zero_mean(model_order + 1 : end);
% [y_hat, w, e] = GLMS(X,y,alpha);
% 
% figure;
% hold on;
% plot(y);
% plot(y_hat);
% xlabel('Sample, n','fontsize',20)
% set(gca,'fontsize',20);
% title('Standard LMS (Zero Mean Signal)','FontSize',20,'FontWeight','bold', 'HorizontalAlignment', 'center')
% set(gca,'fontsize',20)
% grid on
% set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
% legend({'$y[n]$','$\hat{y}[n]$'},'Fontsize', 20,'interpreter','latex');
% %print('../Images/Section_4/Section_4_1_2_M', '-depsc');
% 
% figure;
% plot(e.^2)
% xlabel('Sample, n','fontsize',20)
% ylabel('Squared Error','fontsize',20)
% set(gca,'fontsize',20);
% title('Squared Error (LMS-Zero Mean Signal)','FontSize',20,'FontWeight','bold', 'HorizontalAlignment', 'center')
% set(gca,'fontsize',20)
% grid on
% set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
% %print('../Images/Section_4/Section_4_1_3_M', '-depsc');
% 
% % MSE
% mean_1_full_signal = mean(e.^2);
% mean_1_convergence_section = mean(e(180:end).^2);
% 
% % Prediction Gain
% rp_1_full_signal = 10*log10(var(y_hat)/var(e));
% rp_1_convergence_section = 10*log10(var(y_hat(200:end))/var(e(200:end)));

%-----------------------------------------------------------------------
% Question 2
%-----------------------------------------------------------------------
% %zero-mean plus non-linearity
% data_zero_mean = data-mean(data);
% model_order = 4;
% alpha = 1E-5;
% phi = @(x) tanh(x);
% X = timeseries_to_matrix(data_zero_mean, model_order);
% y = data_zero_mean(model_order + 1 : end);
% [y_hat, w, e] = GLMS(X,y,alpha,phi);
% 
% figure;
% hold on;
% plot(y);
% plot(y_hat);
% 
% xlabel('Sample, n','fontsize',20)
% set(gca,'fontsize',20);
% title('Nonlinear LMS (Zero Mean Signal)','FontSize',20,'FontWeight','bold', 'HorizontalAlignment', 'center')
% set(gca,'fontsize',20)
% grid on
% % set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
% legend({'$y[n]$','$\hat{y}[n]$'},'Fontsize', 20,'interpreter','latex');
% %print('../Images/Section_4/Section_4_2_1_M', '-depsc');
% 
% figure;
% plot(e.^2)
% xlabel('Sample, n','fontsize',20)
% ylabel('Squared Error','fontsize',20)
% set(gca,'fontsize',20);
% title('Squared Error (Non Linear LMS-Zero Mean Signal)','FontSize',20,'FontWeight','bold', 'HorizontalAlignment', 'center')
% set(gca,'fontsize',20)
% grid on
% set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
% %print('../Images/Section_4/Section_4_2_2_M', '-depsc');
% 
% % MSE
% mean_1_full_signal = mean(e.^2);
% mean_1_convergence_section = mean(e(20:end).^2);
% 
% % Prediction Gain
% rp_1_full_signal = 10*log10(var(y_hat)/var(e));
% rp_1_convergence_section = 10*log10(var(y_hat(20:end))/var(e(20:end)));


% %-----------------------------------------------------------------------
% % Question 3
% %-----------------------------------------------------------------------
% %zero-mean plus non-linearity
% data_zero_mean = data-mean(data);
% model_order = 4;
% alpha = 1E-5;
% phi = @(x) 45*tanh(x);
% X = timeseries_to_matrix(data_zero_mean, model_order);
% y = data_zero_mean(model_order + 1 : end);
% [y_hat, w, e] = GLMS(X,y,alpha,phi);
% 
% figure;
% hold on;
% plot(y);
% plot(y_hat);
% 
% xlabel('Sample, n','fontsize',20)
% set(gca,'fontsize',20);
% title('Amplified Nonlinear LMS (Zero Mean Signal)','FontSize',20,'FontWeight','bold', 'HorizontalAlignment', 'center')
% set(gca,'fontsize',20)
% grid on
% % set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
% legend({'$y[n]$','$\hat{y}[n]$'},'Fontsize', 20,'interpreter','latex');
% print('../Images/Section_4/Section_4_3_1_M', '-depsc');
% 
% figure;
% plot(e.^2)
% xlabel('Sample, n','fontsize',20)
% ylabel('Squared Error','fontsize',20)
% set(gca,'fontsize',20);
% title('Squared Error (Non Linear LMS-Zero Mean, Amplified)','FontSize',20,'FontWeight','bold', 'HorizontalAlignment', 'center')
% set(gca,'fontsize',20)
% grid on
% set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
% print('../Images/Section_4/Section_4_3_2_M', '-depsc');
% 
% % MSE
% mean_1_full_signal = mean(e.^2);
% mean_1_convergence_section = mean(e(20:end).^2);
% 
% % Prediction Gain
% rp_1_full_signal = 10*log10(var(y_hat)/var(e));
% rp_1_convergence_section = 10*log10(var(y_hat(20:end))/var(e(20:end)));

% %-----------------------------------------------------------------------
% % Question 4
% %-----------------------------------------------------------------------
%non-zero mean plus non-linearity
% model_order = 4;
% alpha = 1E-5;
% phi = @(x) 45*tanh(x);
% X = timeseries_to_matrix(data, model_order);
% y = data(model_order + 1 : end);
% [y_hat, w, e] = GLMS(X,y,alpha,phi);
% 
% figure;
% hold on;
% plot(y);
% plot(y_hat);
% 
% xlabel('Sample, n','fontsize',20)
% set(gca,'fontsize',20);
% title('Nonlinear LMS (Non-zero Mean Signal, Slow)','FontSize',20,'FontWeight','bold', 'HorizontalAlignment', 'center')
% set(gca,'fontsize',20)
% grid on
% % set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
% legend({'$y[n]$','$\hat{y}[n]$'},'Fontsize', 20,'interpreter','latex');
% print('../Images/Section_4/Section_4_4_1_M', '-depsc');
% 
% 
% figure;
% plot(e.^2)
% xlabel('Sample, n','fontsize',20)
% ylabel('Squared Error','fontsize',20)
% set(gca,'fontsize',20);
% title('Squared Error (Non Linear LMS-Zero Mean, Slow)','FontSize',20,'FontWeight','bold', 'HorizontalAlignment', 'center')
% set(gca,'fontsize',20)
% grid on
% set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
% %print('../Images/Section_4/Section_4_4_2_M', '-depsc');
% 
% % MSE
% mean_1_full_signal = mean(e.^2);
% mean_1_convergence_section = mean(e(180:end).^2);
% 
% % Prediction Gain
% rp_1_full_signal = 10*log10(var(y_hat)/var(e));
% rp_1_convergence_section = 10*log10(var(y_hat(180:end))/var(e(180:end)));
% 

% %-----------------------------------------------------------------------
% % Question 5
% %-----------------------------------------------------------------------
% %non-zero mean plus non-linearity plus pre-training
epochs = 100;
n_pretrain = 20;
model_order = 4;
alpha = 1E-5;
phi = @(x) 45*tanh(x);
X = timeseries_to_matrix(data, model_order);
y = data(model_order + 1 : end);

w_init = zeros(model_order, 1);
for i=1:epochs
    [~, w_init, ~] = GLMS(X(1:n_pretrain,:),y(1:n_pretrain),alpha,phi,w_init);
end
    
[y_hat, w, e] = GLMS(X,y,alpha,phi,w_init);

figure;
hold on;
plot(y);
plot(y_hat);

xlabel('Sample, n','fontsize',20)
set(gca,'fontsize',20);
title('Nonlinear LMS + Pretraining','FontSize',20,'FontWeight','bold', 'HorizontalAlignment', 'center')
set(gca,'fontsize',20)
grid on
% set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
legend({'$y[n]$','$\hat{y}[n]$'},'Fontsize', 20,'interpreter','latex');
%print('../Images/Section_4/Section_4_5_1_M', '-depsc');


figure;
plot(e.^2)
xlabel('Sample, n','fontsize',20)
ylabel('Squared Error','fontsize',20)
set(gca,'fontsize',20);
title('Squared Error (Pretraining)','FontSize',20,'FontWeight','bold', 'HorizontalAlignment', 'center')
set(gca,'fontsize',20)
grid on
set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
%print('../Images/Section_4/Section_4_5_2_M', '-depsc');

% MSE
mean_1_full_signal = mean(e.^2);
mean_1_convergence_section = mean(e(180:end).^2);

% Prediction Gain
rp_1_full_signal = 10*log10(var(y_hat)/var(e));
rp_1_convergence_section = 10*log10(var(y_hat(180:end))/var(e(180:end)));


