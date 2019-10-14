%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% Section 3.3
%-----------------------------------------------------------------------
%
%-----------------------------------------------------------------------
% Question C 
%-----------------------------------------------------------------------
% rng(0);
% mu=1;
% K=1024;
% samp_len=1500;
% samp_freq=1600;
% noise_power=0.05;
% 
% % signal frequencies
% f=zeros(samp_len,1);
% f(1:500)=100;
% f(501:1000)=100+(1:500)./2;
% f(1001:1500)=100+((1:500)./25).^2;
% phase=cumsum(f);
% 
% x=exp(1j*(phase*2*pi/samp_freq));
% w=wgn(samp_len,1,pow2db(0.05),'complex');
% y=x+w;
% [~,~,dft_clms_non_leaky]=dft_clms(y, mu, K, 0);
% [~,~,dft_clms_leaky]=dft_clms(y, mu, K, 0.05);
% 
% % dft_clms plot
% figure(1)
% surf(1:samp_len, [0:K-1].*samp_freq/K, abs(dft_clms_non_leaky), 'LineStyle', 'none');
% view(2);
% plot_func('Time-Frequency Plot: DFT-CLMS', 'n', 'Frequency (Hz)', [0 1500 0 600], 1);
% set(gca, 'Color', 'none');
% set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
% save_figure('../Images/Section_3/Section_3_3_C1_DFT_timefreq');
% 
% % dft_clms leaky plot
% figure(2)
% surf(1:samp_len, [0:K-1].*samp_freq/K, abs(dft_clms_leaky), 'LineStyle', 'none');
% view(2);
% plot_func('Time-Frequency Plot: Leaky DFT-CLMS', 'n', 'Frequency (Hz)', [0 1500 0 600], 1);
% set(gca, 'Color', 'none');
% set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
% save_figure('../Images/Section_3/Section_3_3_C2_DFT_timefreq');
%-----------------------------------------------------------------------
% Question D
%-----------------------------------------------------------------------
% rng(0);
% load Data/colors.mat
% load Data/EEG_Data/EEG_Data/EEG_Data_Assignment2.mat
% 
% % ref signal
% initial_val=1000;
% samp_len=1200;
% y=POz(initial_val:initial_val+samp_len-1);
% mu=1;
% K=8192;
% 
% [~,~,dft_clms_non_leaky]=dft_clms(y, mu, K, 0);
% [~,~,dft_clms_leaky]=dft_clms(y, mu, K, 0.001);
% 
% % dft_clms for EEG
% figure(1)
% surf(1:samp_len, [0:K-1].*fs/K, abs(dft_clms_non_leaky), 'LineStyle', 'none');
% view(2);
% plot_func('Time-Frequency Plot: DFT-CLMS (EEG Data)', 'n', 'Frequency (Hz)', [0 1200 0 100], 1);
% caxis([0 0.0000009])
% set(gca, 'Color', 'none');
% set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
% save_figure('../Images/Section_3/Section_3_3_D1_EEG_DFT_timefreq');
% 
% % Leaky dft_clms for EEG
% figure(2)
% surf(1:samp_len,[0:K-1].*fs/K, abs(dft_clms_leaky), 'LineStyle', 'none');
% view(2);
% plot_func('Time-Frequency Plot: Leaky DFT-CLMS (EEG Data)', 'n', 'Frequency (Hz)', [0 1200 0 100], 1);
% caxis([0 0.0000009])
% set(gca, 'Color', 'none');
% set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
% save_figure('../Images/Section_3/Section_3_3_D2_EEG_DFT_timefreq');

