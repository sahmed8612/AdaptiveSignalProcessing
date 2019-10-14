%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% Section 2.1
%-----------------------------------------------------------------------
%
%-----------------------------------------------------------------------
% Question A  
%-----------------------------------------------------------------------
% rng(0)
% a=[0.1 0.8];
% variance=0.25;
% constant_value=0;
% samp_len=1000;
% no_realisations=100;
% model=arima('Constant',constant_value,'AR',a,'Variance',variance);
% x=simulate(model,samp_len,'samp_lenumPaths',no_realisations);
% error_mu_5=zeros(samp_len,no_realisations);
% error_mu_1=zeros(samp_len,no_realisations);
% order=2;
% 
% for i=1:no_realisations
%     mu=0.05;
%     [~, error_mu_5(:,i), ~] = ar_lms(x(:,i), mu, order);
%     mu=0.01;
%     [~, error_mu_1(:,i), ~] = ar_lms(x(:,i), mu, order);
% end
% 
% x_axis=1:samp_len;
% % 1 Realisation---------------------------------------
% figure(1)
% plot(x_axis, pow2db(error_mu_1(:,1).^2),'LineWidth',2,'color','red');
% hold on;
% plot(x_axis, pow2db(error_mu_5(:,1).^2),'LineWidth',2,'color','blue');
% hold off;
% plot_legend_func('1 Trial: Prediction Error^2', 'n', 'Error^2 (dB)', {'\mu=0.01','\mu=0.05'}, [0 1000 -70 30], 1);
% set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
% save_figure('../Images/Section_2/Section_2_1_A1_Wiener_Filter');
% 
% % 100 Realisation---------------------------------------
% figure(2)
% plot(x_axis, pow2db(mean(error_mu_1.^2,2)),'LineWidth',2,'color','red');
% hold on;
% plot(x_axis, pow2db(mean(error_mu_5.^2,2)),'LineWidth',2,'color','blue');
% hold off;
% plot_legend_func('100 Trials: Prediction Error^2','n', 'Error^2 (dB)', {'\mu=0.01 (\Sigmano_realisations mean)','\mu=0.05 (\Sigmano_realisations mean)'}, [0 1000 -10 0], 1);
% set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
% save_figure('../Images/Section_2/Section_2_1_A2_Wiener_Filter');
%
%-----------------------------------------------------------------------
% Question B
%-----------------------------------------------------------------------

% N = 1000;
% R = 100;
% a = [0.1 0.8];
% sigma_sec = 0.25;
% b = 0;
% p = length(a);
% ar = arima("Constant", b, "AR", a, "Variance", sigma_sec);
% x = simulate(ar, N, "NumPaths", R);
% 
% % step-size and error
% mu = [0.01 0.05];
% error = zeros(N, R, length(mu));
% 
% fig_realisation_1 = figure("Name", "1 Realisation");
% fig_ensmble = figure("Name", sprintf("%d Realisations", R));
% 
% for i=1:length(mu)
%     for j=1:R
%         % data preprocessing
%         [X, y] = Delta_XY(x(:, j), 1, p);
%         % LMS_2
%         [~, error(:, j, i), ~] = LMS_2(X, y, mu(i), 0);
%     end
%     figure(fig_realisation_1);
%     plot(pow2db(error(:, 1, i).^2), "DisplayName", sprintf("mu=%.2f", mu(i)),'LineWidth',2);
%     hold on;
%     figure(fig_ensmble);
%     plot(pow2db(mean(error(:, :, i).^2, 2)), "DisplayName", sprintf("mu=%.2f", mu(i)),'LineWidth',2);
%     hold on;
% end
% 
% figure(fig_realisation_1);
% title("Squared Prediction Error: 1 Realisation",'FontSize',20);
% set(gca,'fontsize',20);
% xlabel("n",'FontSize',20);
% ylabel("Squared Prediction Error (dB)",'FontSize',20);
% axis([0 1000 -70 30]);
% legend("show");
% set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
% save_figure('../Images/Section_2/Section_2_1_B1_time_avg_SS_F');
% 
% figure(fig_ensmble);
% title(sprintf("Squared Prediction Error: %d Realisations", R),'FontSize',20);
% set(gca,'fontsize',20);
% xlabel("n",'FontSize',20);
% ylabel("Squared Prediction Error (dB)",'FontSize',20);
% axis([0 1000 -10 0]);
% legend("show");
% set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
% save_figure('../Images/Section_2/Section_2_1_B2_time_avg_SS_F');

%end--------------------------------------------------------------------

% rng(0);
% a=[0.1 0.8];
% var=0.25;
% const=0;
% samp_len=1000;
% realisations=100;
% model=arima('Constant',const,'AR',a,'Variance',var);
% x=simulate(model,samp_len,'NumPaths',realisations);
% error_mu_5=zeros(samp_len,realisations);
% error_mu_1=zeros(samp_len,realisations);
% model_order=2;
% 
% for i=1:realisations
%     mu=0.05;
%     [~, error_mu_5(:,i), ~] = ar_lms(x(:,i), mu, model_order);
%     mu=0.01;
%     [~, error_mu_1(:,i), ~] = ar_lms(x(:,i), mu, model_order);
% end
% 
% % Calulcating means
% error_mu_1_mean=mean(error_mu_1.^2,2);
% error_mu_5_mean=mean(error_mu_5.^2,2);
% 
% start_index=700;
% x_axis=start_index:samp_len;
% error_mu_1_mean=error_mu_1_mean(start_index:end);
% error_mu_1_time_avg=ones(1,length(x_axis)).*mean(error_mu_1_mean);
% error_mu_5_mean=error_mu_5_mean(start_index:end);
% error_mu_5_time_avg=ones(1,length(x_axis)).*mean(error_mu_5_mean);
% 
% % Plotting Time Average for mu=0.01
% figure(1)
% plot(x_axis, pow2db(error_mu_1_mean),'LineWidth',2,'Color','red');
% hold on;
% plot(x_axis, pow2db(error_mu_1_time_avg),'LineWidth',2,'Color','blue');
% hold off;
% plot_legend_func('mu=0.01: Time Avg', 'n', 'Squared Prediction Error (dB)', {'\mu=0.05 (Total mean - 100 realisation)','Time Avg'}, [300 1000 -7 -5], 1);
% set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
% xlim([700, 1000])
% set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
% save_figure('../Images/Section_2/Section_2_1_B1_time_avg_mu');
% 
% % Plotting Time Average for mu=0.05
% figure(2)
% plot(x_axis, pow2db(error_mu_5_mean),'LineWidth',2,'Color','red');
% hold on;
% plot(x_axis, pow2db(error_mu_5_time_avg),'LineWidth',2,'Color','blue');
% hold off;
% plot_legend_func('mu=0.05: Time Avg', 'n', 'Squared Prediction Error (dB)', {'\mu=0.05 (Total mean - 100 realisation)','Time Avg'}, [300 1000 -7 -5], 1);
% set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
% xlim([700, 1000])
% set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
% save_figure('../Images/Section_2/Section_2_1_B2_time_avg_mu');
%-----------------------------------------------------------------------
% Question D
%-----------------------------------------------------------------------

% rng(0);
% a=[0.1 0.8];
% var=0.25;
% const=0;
% samp_len=1000;
% realisations=100;
% order=2;
% model=arima('Constant',const,'AR',a,'Variance',var);
% x=simulate(model,samp_len,'samp_lenumPaths',realisations);
% 
% % coefficient initialise
% coefficents_ar_mu_5=zeros(realisations*2,samp_len-order+1);
% coefficents_ar_mu_1=zeros(realisations*2,samp_len-order+1);
% order=2;
% 
% for i=1:realisations
%     mu=0.05;
%     [~,~,coefficents_ar_mu_5(1+(i-1)*2:i*2,:)] = ar_lms(x(:,i), mu, order);
%     mu=0.01;
%     [~,~,coefficents_ar_mu_1(1+(i-1)*2:i*2,:)] = ar_lms(x(:,i), mu, order);
% end
% 
% start_index=900;
% x_axis=1:samp_len-order+1;
% 
% coefficents_ar_mu_5_mean(1,:)=mean(coefficents_ar_mu_5(1:2:end,:));
% coefficents_ar_mu_5_mean(2,:)=mean(coefficents_ar_mu_5(2:2:end,:));
% coefficents_ar_mu_5_steady_state(1,:)=mean(coefficents_ar_mu_5_mean(1,start_index:end)).*ones(1,length(x_axis));
% coefficents_ar_mu_5_steady_state(2,:)=mean(coefficents_ar_mu_5_mean(2,start_index:end)).*ones(1,length(x_axis));
% coefficents_ar_mu_1_mean(1,:)=mean(coefficents_ar_mu_1(1:2:end,:));
% coefficents_ar_mu_1_mean(2,:)=mean(coefficents_ar_mu_1(2:2:end,:));
% coefficents_ar_mu_1_steady_state(1,:)=mean(coefficents_ar_mu_1_mean(1,start_index:end)).*ones(1,length(x_axis));
% coefficents_ar_mu_1_steady_state(2,:)=mean(coefficents_ar_mu_1_mean(2,start_index:end)).*ones(1,length(x_axis));
% 
% % mu=0.01, 1 Trial--------------
% figure(1)
% plot(x_axis,coefficents_ar_mu_1(1,:),'LineWidth',2,'color','red');
% hold on;
% plot(x_axis,coefficents_ar_mu_1(2,:),'LineWidth',2,'color','green');
% hold off;
% plot_legend_func('\mu=0.01: 1 Trial', 'n', 'Coefficents Magnitude', {'a_1','a_2'}, [0 1000 -0.5 1.5], 1);
% set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
% save_figure('../Images/Section_2/Section_2_1_D1_time_avg');
% 
% % mu=0.05, 1 Trial--------------
% figure(2)
% plot(x_axis,coefficents_ar_mu_5(1,:),'LineWidth',2,'color','red');
% hold on;
% plot(x_axis,coefficents_ar_mu_5(2,:),'LineWidth',2,'color','green');
% hold off;
% plot_legend_func('\mu=0.05: 1 Trial', 'n', 'Coefficents Magnitude', {'a_1','a_2'}, [0 1000 -0.5 1.5], 1);
% set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
% save_figure('../Images/Section_2/Section_2_1_D2_time_avg');
% 
% % mu=0.01, 100 Trials--------------
% figure(3)
% plot(x_axis,coefficents_ar_mu_1_mean(1,:),'LineWidth',2,'color','red');
% hold on;
% plot(x_axis,coefficents_ar_mu_1_mean(2,:),'LineWidth',2,'color','green');
% hold off;
% plot_legend_func('\mu=0.01: 100 Avg Trials', 'n', 'Coefficents Magnitude', {'a_1','a_2'}, [0 1000 -0.5 1.5], 1);
% set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
% save_figure('../Images/Section_2/Section_2_1_D3_time_avg');
% 
% % mu=0.05, 100 Trials--------------
% figure(4)
% plot(x_axis,coefficents_ar_mu_5_mean(1,:),'LineWidth',2,'color','red');
% hold on;
% plot(x_axis,coefficents_ar_mu_5_mean(2,:),'LineWidth',2,'color','green');
% hold off;
% plot_legend_func('\mu=0.05: 100 Avg Trials', 'n', 'Coefficents Magnitude', {'a_1','a_2'}, [0 1000 -0.5 1.5], 1);
% set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
% save_figure('../Images/Section_2/Section_2_1_D4_time_avg');
% 
% % a_1 for mu=0.01 and mu=0.05 SS values
% figure(5)
% plot(x_axis,coefficents_ar_mu_1_mean(1,:),'LineWidth',2,'color','green');
% hold on;
% plot(x_axis,coefficents_ar_mu_1_steady_state(1,:),'LineWidth',2*2,'color','red');
% plot(x_axis,coefficents_ar_mu_5_mean(1,:),'LineWidth',2,'color','black');
% plot(x_axis,coefficents_ar_mu_5_steady_state(1,:),'LineWidth',2*2,'color','blue');
% hold off;
% plot_legend_func('a_1 Steady State, Varying \mu', 'n', 'Coefficents Magnitude', {'\mu=0.01','\mu=0.01 (Time Avg)','\mu=0.05','\mu=0.05 (Time Avg)'}, [900 1000 0.05 .12], 1);
% set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
% save_figure('../Images/Section_2/Section_2_1_D5_time_avg');
% 
% % a_2 for mu=0.01 and mu=0.05 SS values
% figure(6)
% plot(x_axis,coefficents_ar_mu_1_mean(2,:),'LineWidth',2,'color','green');
% hold on;
% plot(x_axis,coefficents_ar_mu_1_steady_state(2,:),'LineWidth',2*2,'color','red');
% plot(x_axis,coefficents_ar_mu_5_mean(2,:),'LineWidth',2,'color','black');
% plot(x_axis,coefficents_ar_mu_5_steady_state(2,:),'LineWidth',2*2,'color','blue');
% hold off;
% plot_legend_func('a_2 Steady State, Varying \mu', 'n', 'Coefficents Magnitude', {'\mu=0.01','\mu=0.01, (Time Avg)','\mu=0.05','\mu=0.05 (Time Avg)'}, [900 1000 0.7 0.8], 1);
% set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
% save_figure('../Images/Section_2/Section_2_1_D6_time_avg');

%-----------------------------------------------------------------------
% Question F
%----------------------------------------------------------------------
% rng(0);
% load Data/colors.mat;
% 
% gamma=[0.1 0.5 0.9];
% a=[0.1 0.8];
% var=0.25;
% const=0;
% samp_len=1001;
% realisations=100;
% order=2;
% 
% model=arima('Constant',const,'AR',a,'Variance',var);
% x=simulate(model,samp_len,'NumPaths',realisations);
% 
% % coefficient initialise
% coefficents_ar_mu_5_gamma_1=zeros(realisations*2,samp_len-order+1);
% coefficents_ar_mu_5_gamma_5=zeros(realisations*2,samp_len-order+1);
% coefficents_ar_mu_5_gamma_9=zeros(realisations*2,samp_len-order+1);
% coefficents_ar_mu_1_gamma_1=zeros(realisations*2,samp_len-order+1);
% coefficents_ar_mu_1_gamma_5=zeros(realisations*2,samp_len-order+1);
% coefficents_ar_mu_1_gamma_9=zeros(realisations*2,samp_len-order+1);
% 
% 
% order=2;
% 
% for i=1:realisations
%     mu=0.05;
%     [~,~,coefficents_ar_mu_5_gamma_1(1+(i-1)*2:i*2,:)] = leaky_lms(x(:,i), mu, gamma(1), order);
%     [~,~,coefficents_ar_mu_5_gamma_5(1+(i-1)*2:i*2,:)] = leaky_lms(x(:,i), mu, gamma(2), order);
%     [~,~,coefficents_ar_mu_5_gamma_9(1+(i-1)*2:i*2,:)] = leaky_lms(x(:,i), mu, gamma(3), order);
%     mu=0.01;
%     [~,~,coefficents_ar_mu_1_gamma_1(1+(i-1)*2:i*2,:)] = leaky_lms(x(:,i), mu, gamma(1), order);
%     [~,~,coefficents_ar_mu_1_gamma_5(1+(i-1)*2:i*2,:)] = leaky_lms(x(:,i), mu, gamma(2), order);
%     [~,~,coefficents_ar_mu_1_gamma_9(1+(i-1)*2:i*2,:)] = leaky_lms(x(:,i), mu, gamma(3), order);
% end
% 
% % calc means
% coefficents_ar_mu_5_gamma_1_mean(1,:)=mean(coefficents_ar_mu_5_gamma_1(1:2:end,:));
% coefficents_ar_mu_5_gamma_1_mean(2,:)=mean(coefficents_ar_mu_5_gamma_1(2:2:end,:));
% coefficents_ar_mu_5_gamma_5_mean(1,:)=mean(coefficents_ar_mu_5_gamma_5(1:2:end,:));
% coefficents_ar_mu_5_gamma_5_mean(2,:)=mean(coefficents_ar_mu_5_gamma_5(2:2:end,:));
% coefficents_ar_mu_5_gamma_9_mean(1,:)=mean(coefficents_ar_mu_5_gamma_9(1:2:end,:));
% coefficents_ar_mu_5_gamma_9_mean(2,:)=mean(coefficents_ar_mu_5_gamma_9(2:2:end,:));
% 
% coefficents_ar_mu_1_gamma_1_mean(1,:)=mean(coefficents_ar_mu_1_gamma_1(1:2:end,:));
% coefficents_ar_mu_1_gamma_1_mean(2,:)=mean(coefficents_ar_mu_1_gamma_1(2:2:end,:));
% coefficents_ar_mu_1_gamma_5_mean(1,:)=mean(coefficents_ar_mu_1_gamma_5(1:2:end,:));
% coefficents_ar_mu_1_gamma_5_mean(2,:)=mean(coefficents_ar_mu_1_gamma_5(2:2:end,:));
% coefficents_ar_mu_1_gamma_9_mean(1,:)=mean(coefficents_ar_mu_1_gamma_9(1:2:end,:));
% coefficents_ar_mu_1_gamma_9_mean(2,:)=mean(coefficents_ar_mu_1_gamma_9(2:2:end,:));
% 
% x_axis=1:samp_len-order+1;
% a_1_ref=ones(1,length(x_axis))*a(1);
% a_2_ref=ones(1,length(x_axis))*a(2);
% 
% % mu=0.01, gamma=0.1
% figure(1)
% plot(x_axis,coefficents_ar_mu_1_gamma_1_mean(1,:),'LineWidth',2,'color','red');
% hold on;
% plot(x_axis,coefficents_ar_mu_1_gamma_1_mean(2,:),'LineWidth',2,'color','blue');
% plot(x_axis,a_1_ref,'--','LineWidth',2,'Color','green');
% plot(x_axis,a_2_ref,'--','LineWidth',2,'Color','black');
% hold off;
% plot_legend_func('\mu=0.01, \gamma=0.1', 'n', 'Coefficents Magnitude', {'a_1','a_2','a_1 Ref','a_2 Ref'}, [0 1000 0 1], 1);
% set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
% save_figure('../Images/Section_2/Section_2_1_F1_leaky_LMS_2');
% 
% % mu=0.01, gamma=0.5
% figure(2)
% plot(x_axis,coefficents_ar_mu_1_gamma_5_mean(1,:),'LineWidth',2,'color','red');
% hold on;
% plot(x_axis,coefficents_ar_mu_1_gamma_5_mean(2,:),'LineWidth',2,'color','blue');
% plot(x_axis,a_1_ref,'--','LineWidth',2,'Color','green');
% plot(x_axis,a_2_ref,'--','LineWidth',2,'Color','black');
% hold off;
% plot_legend_func('\mu=0.01, \gamma=0.5', 'n', 'Coefficents Magnitude', {'a_1','a_2','a_1 Ref','a_2 Ref'}, [0 1000 0 1], 1);
% set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
% save_figure('../Images/Section_2/Section_2_1_F2_leaky_LMS_2');
% 
% % mu=0.01, gamma=0.9
% figure(3)
% plot(x_axis,coefficents_ar_mu_1_gamma_9_mean(1,:),'LineWidth',2,'color','red');
% hold on;
% plot(x_axis,coefficents_ar_mu_1_gamma_9_mean(2,:),'LineWidth',2,'color','blue');
% plot(x_axis,a_1_ref,'--','LineWidth',2,'Color','green');
% plot(x_axis,a_2_ref,'--','LineWidth',2,'Color','black');
% hold off;
% plot_legend_func('\mu=0.01, \gamma=0.9', 'n', 'Coefficents Magnitude', {'a_1','a_2','a_1 Ref','a_2 Ref'}, [0 1000 0 1], 1);
% set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
% save_figure('../Images/Section_2/Section_2_1_F3_leaky_LMS_2');
% 
% % mu=0.05, gamma=0.1
% figure(4)
% plot(x_axis,coefficents_ar_mu_5_gamma_1_mean(1,:),'LineWidth',2,'color','red');
% hold on;
% plot(x_axis,coefficents_ar_mu_5_gamma_1_mean(2,:),'LineWidth',2,'color','blue');
% plot(x_axis,a_1_ref,'--','LineWidth',2,'Color','green');
% plot(x_axis,a_2_ref,'--','LineWidth',2,'Color','black');
% hold off;
% plot_legend_func('\mu=0.05, \gamma=0.1', 'n', 'Coefficents Magnitude', {'a_1','a_2','a_1 Ref','a_2 Ref'}, [0 1000 0 1], 1);
% set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
% save_figure('../Images/Section_2/Section_2_1_F4_leaky_LMS_2');
% 
% % mu=0.05, gamma=0.5
% figure(5)
% plot(x_axis,coefficents_ar_mu_5_gamma_5_mean(1,:),'LineWidth',2,'color','red');
% hold on;
% plot(x_axis,coefficents_ar_mu_5_gamma_5_mean(2,:),'LineWidth',2,'color','blue');
% plot(x_axis,a_1_ref,'--','LineWidth',2,'Color','green');
% plot(x_axis,a_2_ref,'--','LineWidth',2,'Color','black');
% hold off;
% plot_legend_func('\mu=0.05, \gamma=0.5', 'n', 'Coefficents Magnitude', {'a_1','a_2','a_1 Ref','a_2 Ref'}, [0 1000 0 1], 1);
% set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
% save_figure('../Images/Section_2/Section_2_1_F5_leaky_LMS_2');
% 
% % mu=0.05, gamma=0.9
% figure(6)
% plot(x_axis,coefficents_ar_mu_5_gamma_9_mean(1,:),'LineWidth',2,'color','red');
% hold on;
% plot(x_axis,coefficents_ar_mu_5_gamma_9_mean(2,:),'LineWidth',2,'color','blue');
% plot(x_axis,a_1_ref,'--','LineWidth',2,'Color','green');
% plot(x_axis,a_2_ref,'--','LineWidth',2,'Color','black');
% hold off;
% plot_legend_func('\mu=0.05, \gamma=0.9', 'n', 'Coefficents Magnitude', {'a_1','a_2','a_1 Ref','a_2 Ref'}, [0 1000 0 1], 1);
% set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
% save_figure('../Images/Section_2/Section_2_1_F6_leaky_LMS_2');

