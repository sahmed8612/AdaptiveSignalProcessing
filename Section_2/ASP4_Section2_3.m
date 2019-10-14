%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% Section 2.3
%-----------------------------------------------------------------------
%
%-----------------------------------------------------------------------
% Question A  
%-----------------------------------------------------------------------
% rng(0);
% load Data/colors.mat
% 
% mu=0.01;
% samp_len=1000;
% t=1:samp_len;
% freq_0=0.005;
% sine_signalwave_sig = sin(2*pi*freq_0*t);
% noise_power = 1;
% b=[1 0 0.5];
% a=1;
% realisations=100;
% y=zeros(samp_len,realisations);
% 
% x_hat_delay_1=zeros(realisations,samp_len);
% x_hat_delay_2=zeros(realisations,samp_len);
% x_hat_delay_3=zeros(realisations,samp_len);
% x_hat_delay_4=zeros(realisations,samp_len);
% 
% error_delay_1=zeros(realisations,samp_len-1);
% error_delay_2=zeros(realisations,samp_len-1);
% error_delay_3=zeros(realisations,samp_len-1);
% error_delay_4=zeros(realisations,samp_len-1);
% 
% mpse=zeros(realisations,4);
% 
% % model order
% order=1;
% 
% for i=1:realisations
%     w=noise_generator(samp_len,noise_power);
%     filtered_noise=filter(b,a,w);
%     y(:,i)=sine_signalwave_sig'+filtered_noise;
%     
%     delay=1;
%     [x_hat_delay_1(i,:), error_delay_1(i,:)] = ale_lms(y(:,i),mu,delay,order);
%     mpse(i,1)=sum((sine_signalwave_sig-x_hat_delay_1(i,:)).^2)/samp_len;
%     delay=2;
%     [x_hat_delay_2(i,:), error_delay_2(i,:)] = ale_lms(y(:,i),mu,delay,order);
%     mpse(i,2)=sum((sine_signalwave_sig-x_hat_delay_2(i,:)).^2)/samp_len;
%     delay=3;
%     [x_hat_delay_3(i,:), error_delay_3(i,:)] = ale_lms(y(:,i),mu,delay,order);
%     mpse(i,3)=sum((sine_signalwave_sig-x_hat_delay_3(i,:)).^2)/samp_len;
%     delay=4;
%     [x_hat_delay_4(i,:), error_delay_4(i,:)] = ale_lms(y(:,i),mu,delay,order);
%     mpse(i,4)=sum((sine_signalwave_sig-x_hat_delay_4(i,:)).^2)/samp_len;
% end
% 
% % means calc
% delay_1_mean_x_hat=mean(x_hat_delay_1);
% delay_2_mean_x_hat=mean(x_hat_delay_2);
% delay_3_mean_x_hat=mean(x_hat_delay_3);
% delay_4_mean_x_hat=mean(x_hat_delay_4);
% 
% delay_1_mean_error=mean(error_delay_1);
% delay_2_mean_error=mean(error_delay_2);
% delay_3_mean_error=mean(error_delay_3);
% delay_4_mean_error=mean(error_delay_4);
% 
% mpse_mean=mean(mpse);
% 
% % x-hat for Delay = 1--------------------------------
% figure(1)
% plot(t,delay_1_mean_x_hat,'LineWidth',2,'color','red');
% plot_func('100 Realisations Avg: \Delta = 1','n','Amplitude', [0 1000 -1 1], 1);
% set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
% save_figure('../Images/Section_2/Section_2_3_A1_Ale_Delay');
% 
% % x-hat for Delay = 2--------------------------------
% figure(2)
% plot(t,delay_2_mean_x_hat,'LineWidth',2,'color','red');
% plot_func('100 Realisations Avg: \Delta = 2','n','Amplitude', [0 1000 -1 1], 1);
% set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
% save_figure('../Images/Section_2/Section_2_3_A2_Ale_Delay');
% 
% % x-hat for Delay = 3--------------------------------
% figure(3)
% plot(t,delay_3_mean_x_hat,'LineWidth',2,'color','red');
% plot_func('100 Realisations Avg: \Delta = 3','n','Amplitude', [0 1000 -1 1], 1);
% set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
% save_figure('../Images/Section_2/Section_2_3_A3_Ale_Delay');
% 
% % x-hat for Delay = 4--------------------------------
% figure(4)
% plot(t,delay_4_mean_x_hat,'LineWidth',2,'color','red');
% plot_func('100 Realisations Avg: \Delta = 4','n','Amplitude', [0 1000 -1 1], 1);
% set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
% save_figure('../Images/Section_2/Section_2_3_A4_Ale_Delay');
% 
% % MPSE for Delay = 1--------------------------------
% figure(5)
% for i=1:realisations
%     sine_signal_noisey=plot(y(:,i), 'LineWidth', 2, 'Color', 'blue');
%     if i==1
%         hold on;
%     end
% end
% for i=1:realisations
%     sine_signal_denoised=plot(x_hat_delay_1(i,:), 'LineWidth', 2, 'Color', 'red');
% end
% sine_signal_clean=plot(sine_signalwave_sig,'LineWidth',2*2, 'Color', 'green');
% hold off;
% str=sprintf('Average MPSE = %.3f', mpse_mean(1));
% plot_func(str,'n','Amplitude', [0 1000 -6 6], 1);
% legend([sine_signal_noisey sine_signal_denoised sine_signal_clean],{'Sine-Noisy','Sine-Denoised','Sine-Clean'},'FontSize',15);
% set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
% save_figure('../Images/Section_2/Section_2_3_A5_Ale_Delay');
% 
% % MPSE for Delay = 2--------------------------------
% figure(6)
% for i=1:realisations
%     sine_signal_noisey=plot(y(:,i), 'LineWidth', 2, 'Color', 'blue');
%     if i==1
%         hold on;
%     end
% end
% for i=1:realisations
%     sine_signal_denoised=plot(x_hat_delay_2(i,:), 'LineWidth', 2, 'Color', 'red');
% end
% sine_signal_clean=plot(sine_signalwave_sig,'LineWidth',2*2, 'Color', 'green');
% hold off;
% str=sprintf('Average MPSE = %.3f', mpse_mean(2));
% plot_func(str,'n','Amplitude', [0 1000 -6 6], 1);
% legend([sine_signal_noisey sine_signal_denoised sine_signal_clean],{'Sine-Noisy','Sine-Denoised','Sine-Clean'},'FontSize',15);
% set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
% save_figure('../Images/Section_2/Section_2_3_A6_Ale_Delay');
% 
% % MPSE for Delay = 3--------------------------------
% figure(7)
% for i=1:realisations
%     sine_signal_noisey=plot(y(:,i), 'LineWidth', 2, 'Color', 'blue');
%     if i==1
%         hold on;
%     end
% end
% for i=1:realisations
%     sine_signal_denoised=plot(x_hat_delay_3(i,:), 'LineWidth', 2, 'Color', 'red');
% end
% sine_signal_clean=plot(sine_signalwave_sig,'LineWidth',2*2, 'Color', 'green');
% hold off;
% str=sprintf('Average MPSE = %.3f', mpse_mean(3));
% plot_func(str,'n','Amplitude', [0 1000 -6 6], 1);
% legend([sine_signal_noisey sine_signal_denoised sine_signal_clean],{'Sine-Noisy','Sine-Denoised','Sine-Clean'},'FontSize',15);
% set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
% save_figure('../Images/Section_2/Section_2_3_A7_Ale_Delay');
% 
% % MPSE for Delay = 4--------------------------------
% figure(8)
% for i=1:realisations
%     sine_signal_noisey=plot(y(:,i), 'LineWidth', 2, 'Color', 'blue');
%     if i==1
%         hold on;
%     end
% end
% for i=1:realisations
%     sine_signal_denoised=plot(x_hat_delay_4(i,:), 'LineWidth', 2, 'Color', 'red');
% end
% sine_signal_clean=plot(sine_signalwave_sig,'LineWidth',2*2, 'Color', 'green');
% hold off;
% str=sprintf('Average MPSE = %.3f', mpse_mean(4));
% plot_func(str,'n','Amplitude', [0 1000 -6 6], 1);
% legend([sine_signal_noisey sine_signal_denoised sine_signal_clean],{'Sine-Noisy','Sine-Denoised','Sine-Clean'},'FontSize',15);
% set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
% save_figure('../Images/Section_2/Section_2_3_A8_Ale_Delay');
%-----------------------------------------------------------------------
% Question B  
%-----------------------------------------------------------------------
% rng(0);
% order=[5,10,15,20];
% delay_len=1:25;
% mu=0.01;
% N=1000;
% t=1:N;
% freq_0=0.005;
% sine_signalwave_sig = sin(2*pi*freq_0*t);
% noise_power = 1;
% b=[1 0 0.5];
% a=1;
% realisations=100;
% 
% mpse=zeros(realisations,length(delay_len),length(order));
% 
% for i=1:realisations
%     w=noise_generator(N,noise_power);
%     filtered_noise=filter(b,a,w);
%     y=sine_signalwave_sig'+filtered_noise;
%     for j=1:length(order)
%         for k=1:length(delay_len)
%            x_hat=ale_lms(y,mu,delay_len(k),order(j));
%            mpse(i,k,j)=sum((sine_signalwave_sig-x_hat).^2)./N;
%         end
%     end
% end
% mpse_mean=mean(mpse);
% 
% % MPSE and Delays
% figure(1);
% plot(delay_len,mpse_mean(:,:,1),'-','LineWidth',2,'color','red');
% hold on;
% plot(delay_len,mpse_mean(:,:,2),'-','LineWidth',2,'color','blue');
% plot(delay_len,mpse_mean(:,:,3),'-','LineWidth',2,'color','green');
% plot(delay_len,mpse_mean(:,:,4),'-','LineWidth',2,'color','cyan');
% hold off;
% plot_legend_func('MPSE When Varying \Delta','\Delta','MPSE', {'M=5','M=10','M=15','M=20'}, [1 25 0.2 0.8], 1);
% set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
% save_figure('../Images/Section_2/Section_2_3_B1_Delay_vs_MPSE');

% %For large lag --------------------------------
% rng(0);
% order=5;
% mu=0.01;
% N=1000;
% t=1:N;
% freq_0=0.005;
% sine_signal = sin(2*pi*freq_0*t);
% 
% noise_power = 1;
% b=[1 0 0.5];
% a=1;
% realisations=100;
% y=zeros(N,realisations);
% delay_5_x_hat=zeros(realisations,N);
% delay_25_x_hat=zeros(realisations,N);
% 
% mpse=zeros(realisations,2);
% 
% for i=1:realisations
%     w=noise_generator(N,noise_power);
%     filtered_noise=filter(b,a,w);
%     y(:,i)=sine_signal'+filtered_noise;
%     
%     delay=5;
%     delay_5_x_hat(i,:) = ale_lms(y(:,i),mu,delay,order);
%     mpse(i,1)=sum((sine_signal-delay_5_x_hat(i,:)).^2)/N;
%     delay=25;
%     delay_25_x_hat(i,:) = ale_lms(y(:,i),mu,delay,order);
%     mpse(i,2)=sum((sine_signal-delay_25_x_hat(i,:)).^2)/N;
% 
% end
% mpse_mean=mean(mpse);
% 
% % MPSE for Delay = 5--------------------------------
% figure(1)
% for i=1:10
%     sine_signal_noisey=plot(y(:,i), 'LineWidth', 2, 'Color', 'blue');
%     if i==1
%         hold on;
%     end
% end
% for i=1:10
%     sine_signal_denoised=plot(delay_5_x_hat(i,:), 'LineWidth', 2, 'Color', 'red');
% end
% sine_signal_clean=plot(sine_signal,'LineWidth',2*2, 'Color', 'green');
% hold off;
% str=sprintf('Increasing Delay Effects:  MPSE = %.3f', mpse_mean(1));
% str=strcat('M=5, \Delta = 5,',str);
% plot_func(str,'n','Amplitude', [0 1000 -6 6], 1);
% legend([sine_signal_noisey sine_signal_denoised sine_signal_clean],{'Sine-Noisy','Sine-Denoised','Sine-Clean'},'FontSize',15);
% set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
% save_figure('../Images/Section_2/Section_2_3_B2_Delay_vs_MPSE')
% 
% % MPSE for Delay = 25--------------------------------
% figure(2)
% for i=1:10
%     sine_signal_noisey=plot(y(:,i), 'LineWidth', 2, 'Color', 'blue');
%     if i==1
%         hold on;
%     end
% end
% for i=1:10
%     sine_signal_denoised=plot(delay_25_x_hat(i,:), 'LineWidth', 2, 'Color', 'red');
% end
% sine_signal_clean=plot(sine_signal,'LineWidth',2*2, 'Color', 'green');
% hold off;
% str=sprintf(' Increasing Delay Effects: MPSE = %.3f', mpse_mean(2));
% str=strcat('M=5, \Delta = 25,',str);
% plot_func(str,'n','Amplitude', [0 1000 -6 6], 1);
% legend([sine_signal_noisey sine_signal_denoised sine_signal_clean],{'Sine-Noisy','Sine-Denoised','Sine-Clean'},'FontSize',15);
% set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
% save_figure('../Images/Section_2/Section_2_3_B3_Delay_vs_MPSE')

% %For Appropriate lag --------------------------------

% rng(0);
% order=1:15;
% delays=3;
% mu=[0.005,0.01,0.015,0.020,0.025];
% 
% N=1000;
% t=1:N;
% freq_0=0.005;
% sine_signal = sin(2*pi*freq_0*t);
% noise_power = 1;
% b=[1 0 0.5];
% a=1;
% realisations=100;
% 
% mpse=zeros(realisations,length(order),length(mu));
% 
% for i=1:realisations
%     w=noise_generator(N,noise_power);
%     filtered_noise=filter(b,a,w);
%     y=sine_signal'+filtered_noise;
%     for j=1:length(order)
%         for k=1:length(mu)
%            x_hat=ale_lms(y,mu(k),delays,order(j));
%            mpse(i,j,k)=sum((sine_signal-x_hat).^2)./N;
%         end
%     end
% end
% mpse_mean=mean(mpse);
% 
% % MPSE & Delays Realtionship
% figure(1);
% plot(order,mpse_mean(:,:,1),'-','LineWidth',2,'color','red');
% hold on;
% plot(order,mpse_mean(:,:,2),'-','LineWidth',2,'color','blue');
% plot(order,mpse_mean(:,:,3),'-','LineWidth',2,'color','green');
% plot(order,mpse_mean(:,:,4),'-','LineWidth',2,'color','yellow');
% plot(order,mpse_mean(:,:,5),'-','LineWidth',2,'color','cyan');
% hold off;
% plot_legend_func('MPSE & M (Delays) Relationship','M','MPSE', {'/mu=0.005','/mu=0.010','/mu=0.015','/mu=0.020','/mu=0.025'}, [1 15 0 1.5], 1);
% set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
% save_figure('../Images/Section_2/Section_2_3_B4_Delay_vs_MPSE')

%-----------------------------------------------------------------------
% Question C
%-----------------------------------------------------------------------

% rng(0);
% load Data/colors.mat
% order=6;
% delay=3;
% mu=0.01;
% N=1000;
% t=1:N;
% freq_0=0.005;
% sine_signal = sin(2*pi*freq_0*t);
% 
% noise_power = 1;
% b=[1 0 0.5];
% a=1;
% 
% 
% realisations=100;
% y=zeros(N,realisations);
% ale_x_hat=zeros(realisations,N);
% anc_x_hat=zeros(realisations,N);
% ale_errors=zeros(realisations,N-1);
% anc_errors=zeros(realisations,N-1);
% 
% mpse=zeros(realisations,2);
% 
% for i=1:realisations
%     w=noise_generator(N,noise_power);
%     filtered_noise=filter(b,a,w);
%     y(:,i)=sine_signal'+filtered_noise;
%     
%    [ale_x_hat(i,:),ale_errors(i,:)]=ale_lms(y(:,i),mu,delay,order);
%    mpse(i,1)=sum((sine_signal-ale_x_hat(i,:)).^2)./N;
%    
%    [n_hat,anc_errors(i,:)]=anc_lms(y(:,i),filtered_noise,mu,order);
%    anc_x_hat(i,:)=y(:,i)'-n_hat;
%    mpse(i,2)=sum((sine_signal-anc_x_hat(i,:)).^2)./N;
%    
% end
% 
% % get means
% ale_x_hat_mean=mean(ale_x_hat);
% anc_x_hat_mean=mean(anc_x_hat);
% mpse_mean=mean(mpse);
% 
% % generate x-axis for plotting
% x_axis=1:N;
% 
% % Plot The Average Signals Obtained using ALE and ANC
% figure(1);
% plot(x_axis,ale_x_hat_mean,'LineWidth',2, 'Color', 'blue');
% hold on;
% plot(x_axis,anc_x_hat_mean,'LineWidth',2, 'Color', 'red');
% hold off;
% plot_legend_func('ALE & ANC Comparison', 'n', 'Amplitude', {'ALE Average Output','ANC Average Output'}, [0 1000 -2 2], 1);
% set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
% save_figure('../Images/Section_2/Section_2_3_C1_ALE_ANC_Comparison')
% 
% % Plotting of MPSE for ALE, Averaged Over 100 Realisations
% figure(2)
% for i=1:10
%     sine_signal_noisey=plot(y(:,i), 'LineWidth', 2, 'Color', 'blue');
%     if i==1
%         hold on;
%     end
% end
% for i=1:10
%     sine_signal_denoised=plot(ale_x_hat(i,:), 'LineWidth', 2, 'Color', 'red');
% end
% sine_signal_clean=plot(sine_signal,'LineWidth',2*2, 'Color', 'green');
% hold off;
% str=sprintf(' When MPSE = %.3f', mpse_mean(1));
% str=strcat('ALE', str);
% plot_func(str,'n','Amplitude', [0 1000 -6 6], 1);
% legend([sine_signal_noisey sine_signal_denoised sine_signal_clean],{'Sine-Noisy','Sine-Denoised','Sine-Clean'},'FontSize',15);
% set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
% save_figure('../Images/Section_2/Section_2_3_C2_ALE_ANC_Comparison')
% 
% % ANC for MPSE
% figure(3)
% for i=1:10
%     sine_signal_noisey=plot(y(:,i), 'LineWidth', 2, 'Color', 'blue');
%     if i==1
%         hold on;
%     end
% end
% for i=1:10
%     sine_signal_denoised=plot(anc_x_hat(i,:), 'LineWidth', 2, 'Color', 'red');
% end
% sine_signal_clean=plot(sine_signal,'LineWidth',2*2, 'Color', 'green');
% hold off;
% str=sprintf(' When MPSE = %.3f', mpse_mean(2));
% str=strcat('ANC', str);
% plot_func(str,'n','Amplitude', [0 1000 -6 6], 1);
% legend([sine_signal_noisey sine_signal_denoised sine_signal_clean],{'Sine-Noisy','Sine-Denoised','Sine-Clean'},'FontSize',15);
% set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
% save_figure('../Images/Section_2/Section_2_3_C3_ALE_ANC_Comparison')

%-----------------------------------------------------------------------
% Question D
%-----------------------------------------------------------------------

rng(0);
load Data/EEG_Data/EEG_Data/EEG_Data_Assignment2.mat;
load Data/colors.mat

freq_samp = 1200;
mu=[0.025 0.01 0.005 0.001];
order=[10 15 20 25];
K=16382;
N=length(POz);
L=4096;
overlap_const=0.8; 
t=1:N;
freq_0=50;
noise_power=1;
% to make sampling period is 1, need normalize

f=freq_0/freq_samp;
sine_signal=sin(2*pi*f*t);
w=noise_generator(N,noise_power);
y=sine_signal'+w;
POz=POz-mean(POz);

P0z_noiseless=zeros(N,length(mu),length(order));

for i=1:length(mu)
    for j=1:length(order)
        [noise_predictor]=anc_lms(POz,y,mu(i),order(j));
        P0z_noiseless(:,i,j)=POz-noise_predictor';
    end
end

% orginal spectrogram
figure(1);
spectrogram(POz, rectwin(L),round(overlap_const*L),K,freq_samp,'yaxis');
ylim([0 100]);
set(gca,'fontsize',20);
title('POz Original', 'FontSize', 20);
grid off;
save_figure('../Images/Section_2/Section_2_3_D1_EEG_Denoising')

% mu =  0.01 and M=10
figure(2);
spectrogram(P0z_noiseless(:,2,1), rectwin(L),round(overlap_const*L),K,freq_samp,'yaxis');
ylim([0 100]);
set(gca,'fontsize',20);
title('\mu = 0.01, M = 10', 'FontSize', 20*32/24);
grid off;
save_figure('../Images/Section_2/Section_2_3_D2_EEG_Denoising')

% mu =  0.01 and M=15
figure(3);
spectrogram(P0z_noiseless(:,2,2), rectwin(L),round(overlap_const*L),K,freq_samp,'yaxis');
ylim([0 100]);
set(gca,'fontsize',20);
title('\mu = 0.01, M = 15', 'FontSize', 20*32/24);
grid off;
save_figure('../Images/Section_2/Section_2_3_D3_EEG_Denoising')

% plot spectrogram for mu =  0.01 and M=20
figure(4);
spectrogram(P0z_noiseless(:,2,3), rectwin(L),round(overlap_const*L),K,freq_samp,'yaxis');
ylim([0 100]);
set(gca,'fontsize',20);
title('\mu = 0.01, M = 20', 'FontSize', 20*32/24);
grid off;
save_figure('../Images/Section_2/Section_2_3_D4_EEG_Denoising')

% mu =  0.01 and M=25
figure(5);
spectrogram(P0z_noiseless(:,2,4), rectwin(L),round(overlap_const*L),K,freq_samp,'yaxis');
ylim([0 100]);
set(gca,'fontsize',20);
title('\mu = 0.01, M = 25', 'FontSize', 20*32/24);
grid off;
save_figure('../Images/Section_2/Section_2_3_D5_EEG_Denoising')

% mu =  0.025 and M=25
figure(6);
spectrogram(P0z_noiseless(:,1,4), rectwin(L),round(overlap_const*L),K,freq_samp,'yaxis');
ylim([0 100]);
set(gca,'fontsize',20);
title('\mu = 0.025, M = 25', 'FontSize', 20*32/24);
grid off;
save_figure('../Images/Section_2/Section_2_3_D6_EEG_Denoising')

% mu =  0.01 and M=25
figure(7);
spectrogram(P0z_noiseless(:,2,4), rectwin(L),round(overlap_const*L),K,freq_samp,'yaxis');
ylim([0 100]);
set(gca,'fontsize',20);
title('\mu = 0.01, M = 25', 'FontSize', 20*32/24);
grid off;
save_figure('../Images/Section_2/Section_2_3_D7_EEG_Denoising')

% mu =  0.005 and M=25
figure(8);
spectrogram(P0z_noiseless(:,3,4), rectwin(L),round(overlap_const*L),K,freq_samp,'yaxis');
ylim([0 100]);
set(gca,'fontsize',20);
title('\mu = 0.005, M = 25', 'FontSize', 20*32/24);
grid off;
save_figure('../Images/Section_2/Section_2_3_D8_EEG_Denoising')

% mu =  0.001 and M=25
figure(9);
spectrogram(P0z_noiseless(:,4,4), rectwin(L),round(overlap_const*L),K,freq_samp,'yaxis');
ylim([0 100]);
set(gca,'fontsize',20);
title('\mu = 0.001, M = 25', 'FontSize', 20*32/24);
grid off;
save_figure('../Images/Section_2/Section_2_3_D9_EEG_Denoising')

% Bartlett Periodograms, averaged over 10seconds, for both Cleaned and Noisy POz Components
figure(10)
num_seconds_average_over=2;
L = freq_samp*num_seconds_average_over;
noverlap_const = 0;
[pxx,w]= pwelch(POz,rectwin(L),noverlap_const,K,freq_samp,'onesided');
[pxx_cleaned,w_cleaned]= pwelch(P0z_noiseless(:,4,4),rectwin(L),noverlap_const,K,freq_samp,'onesided');
[pxx_semi_cleaned,w_semi_cleaned]= pwelch(P0z_noiseless(:,1,4),rectwin(L),noverlap_const,K,freq_samp,'onesided');
plot(w,pow2db(pxx),'LineWidth',2,'Color','red');
hold on;
plot(w_cleaned,pow2db(pxx_cleaned),'LineWidth',2,'Color','blue');
plot(w_semi_cleaned,pow2db(pxx_semi_cleaned),'LineWidth',2,'Color','green')
hold off;
plot_legend_func('Bartlett EEG Periodogram Spectrum', 'Frequency (Hz)', 'Power/Freq (dB/Hz)', {'POz Original', 'mu=0.001 (M=25)', 'mu=0.025 (M=25)'}, [0 60 -140 -90], 1);
set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
save_figure('../Images/Section_2/Section_2_3_D10_EEG_Denoising')

% Bartlett Periodograms, mu=0.001, M=25
figure(11)
error_mu_001_M_25=(pxx-pxx_cleaned).^2;
plot(w,pow2db(error_mu_001_M_25),'LineWidth',2, 'Color', 'red');
plot_func('Error^2 When \mu=0.001 & M=25', 'Frequency (Hz)', 'Error Power (dB)', [0 60 -310 -190], 1);
set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
save_figure('../Images/Section_2/Section_2_3_D11_EEG_Denoising')

% Bartlett Periodograms when mu=0.025, M=25
figure(12)
error_mu_01_M_25=(pxx-pxx_semi_cleaned).^2;
plot(w,pow2db(error_mu_01_M_25),'LineWidth',2, 'Color', 'red');
plot_func('Error^2 When \mu=0.025 & M=25', 'Frequency (Hz)', 'Error Power (dB)', [0 60 -310 -190], 1);
set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
save_figure('../Images/Section_2/Section_2_3_D12_EEG_Denoising')