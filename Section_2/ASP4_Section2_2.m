%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% Section 2.2
%-----------------------------------------------------------------------
%
%-----------------------------------------------------------------------
% Question A  
%-----------------------------------------------------------------------
% rng(1);
% rho = 0.005;
% mu = 0;
% b = [1 0.9];
% a = 1;
% noise_power=0.5;
% samp_len = 1000;
% realisations=100;
% 
% coeffecients_ben=zeros(samp_len, realisations);
% coeffecients_ang=zeros(samp_len, realisations);
% coeffecients_mat=zeros(samp_len, realisations);
% coeffecients_mu_1=zeros(samp_len, realisations);
% coeffecients_mu_01=zeros(samp_len, realisations);
% 
% % model order parameter
% order = 1;
% 
% % run realisations
% for i=1:realisations
%     % generate filtered noise
%     w=noise_generator(samp_len,noise_power);
%     y=filter(b,a,w);
%     % obtain coefficients using different algorithms
%     [~, ~,coeffecients_ben(:,i)] = gass_lms(w,y, mu, rho, order,'ben');
%     [~, ~,coeffecients_ang(:,i)] = gass_lms(w,y, mu, rho, order,'ang');
%     [~, ~,coeffecients_mat(:,i)] = gass_lms(w,y, mu, rho, order,'mat');
%     mu=0.1;
%     [~, ~,coeffecients_mu_1(:,i)] = lms(w,y, mu, order);
%     mu=0.01;
%     [~, ~,coeffecients_mu_01(:,i)] = lms(w,y, mu, order);
% end
% 
% % mean cal
% coeffecients_ben_mean=mean(coeffecients_ben,2);
% coeffecients_ang_mean=mean(coeffecients_ang,2);
% coeffecients_mat_mean=mean(coeffecients_mat,2);
% coeffecients_mu_1_mean=mean(coeffecients_mu_1,2);
% coeffecients_mu_01_mean=mean(coeffecients_mu_01,2);
% 
% figure(1)
% plot(b(2)-coeffecients_mu_01(:,1),'LineWidth',2,'color','red');
% hold on;
% plot(b(2)-coeffecients_mu_1(:,1),'LineWidth',2,'color','blue');
% plot(b(2)-coeffecients_ben(:,1),'LineWidth',2,'color','green');
% plot(b(2)-coeffecients_ang(:,1),'LineWidth',2,'color','magenta');
% plot(b(2)-coeffecients_mat(:,1),'LineWidth',2,'color','cyan');
% hold off;
% plot_legend_func('Coefficients Convergence', 'n', 'Weight Error', {'\mu=0.001 (LMS)', '\mu=0.01 (LMS)', 'Benveniste', 'Ang and Farhang', 'Matthew and Xie'}, [0 1000 -1 1.5], 1);
% set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
% save_figure('../Images/Section_2/Section_2_2_A1_Weights_Converge');
% 
% figure(2)
% plot(pow2db((b(2)-coeffecients_mu_01(:,1)).^2),'LineWidth',2,'color','red');
% hold on;
% plot(pow2db((b(2)-coeffecients_mu_1(:,1)).^2),'LineWidth',2,'color','blue');
% plot(pow2db((b(2)-coeffecients_ben(:,1)).^2),'LineWidth',2,'color','green');
% plot(pow2db((b(2)-coeffecients_ang(:,1)).^2),'LineWidth',2,'color','magenta');
% plot(pow2db((b(2)-coeffecients_mat(:,1)).^2),'LineWidth',2,'color','cyan');
% hold off;
% plot_legend_func('Error^2 (dB)', 'n', 'Weight Error', {'\mu=0.001 (LMS)', '\mu=0.01 (LMS)', 'Benveniste', 'Ang and Farhang', 'Matthew and Xie'}, [0 1000 -300 50], 1);
% set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
% save_figure('../Images/Section_2/Section_2_2_A2_Weights_Converge');
% 
% figure(3)
% plot(b(2)-coeffecients_mu_01_mean,'LineWidth',2,'color','red');
% hold on;
% plot(b(2)-coeffecients_mu_1_mean,'LineWidth',2,'color','blue');
% plot(b(2)-coeffecients_ben_mean,'LineWidth',2,'color','green');
% plot(b(2)-coeffecients_ang_mean,'LineWidth',2,'color','magenta');
% plot(b(2)-coeffecients_mat_mean,'LineWidth',2,'color','cyan');
% hold off;
% plot_legend_func('100 Trials: Avg Weights', 'n', 'Weight Error', {'\mu=0.001 (LMS)', '\mu=0.01 (LMS)', 'Benveniste', 'Ang and Farhang', 'Matthew and Xie'}, [0 1000 -0.2 1], 1);
% set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
% save_figure('../Images/Section_2/Section_2_2_A3_Weights_Converge');
% 
% figure(4)
% plot(pow2db((b(2)-coeffecients_mu_01_mean).^2),'LineWidth',2,'color','red');
% hold on;
% plot(pow2db((b(2)-coeffecients_mu_1_mean).^2),'LineWidth',2,'color','blue');
% plot(pow2db((b(2)-coeffecients_ben_mean).^2),'LineWidth',2,'color','green');
% plot(pow2db((b(2)-coeffecients_ang_mean).^2),'LineWidth',2,'color','magenta');
% plot(pow2db((b(2)-coeffecients_mat_mean).^2),'LineWidth',2,'color','cyan');
% hold off;
% plot_legend_func('100 Trials: Error^2 (dB)', 'n', 'Weight Error', {'\mu=0.001 (LMS)', '\mu=0.01 (LMS)', 'Benveniste', 'Ang and Farhang', 'Matthew and Xie'}, [0 1000 -300 50], 1);
% set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
% save_figure('../Images/Section_2/Section_2_2_A4_Weights_Converge');

%-----------------------------------------------------------------------
%Question B
%-----------------------------------------------------------------------
% rng(1);
% order=1;
% mu=0.01;
% rho=0.005;
% algo='ben';
% N=1000;
% noise_power=0.5;
% b=[1 0.9];
% a=1;
% w=noise_generator(N,noise_power);
% y=filter(b,a,w);
% 
% [~,~,coeffecients_ben] = gass_lms(w,y,1,rho,order,'ben');
% [~,~,coeffecients_GNGD] = GNGD_lms(w,y,10,rho,order);
% 
% figure(1);
% plot(coeffecients_ben,'LineWidth',2,'color','red');
% hold on;
% plot(coeffecients_GNGD,'LineWidth',2,'color','green');
% hold off;
% plot_legend_func('Comparison between GNGD & GASS', 'n', 'Coefficients Magnitude', {'\mu=1 (Benveniste), ', '\mu=10 (GNGD)'}, [0 30 -0.2 1.2], 1);
% set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
% save_figure('../Images/Section_2/Section_2_2_B1_Compare_GNGD');