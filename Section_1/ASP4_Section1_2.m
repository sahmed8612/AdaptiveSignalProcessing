%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% Section 1.2
%-----------------------------------------------------------------------
%
%-----------------------------------------------------------------------
% Question A  
%-----------------------------------------------------------------------

% % load sunspot.dat;
% % year = sunspot(:,1);
% % relNums = sunspot(:,2);
% % 
% % % length of fft
% % K=2048;
% % % no pre-processing
% % xf=abs(fftshift(fft([relNums' zeros(1, K-length(relNums))])));
% % Pxx1=pow2db(xf.^2/(length(relNums)*2*pi))';
% % % removing the mean
% % relNums_without_mean=relNums-mean(relNums);
% % xf=abs(fftshift(fft([relNums_without_mean' zeros(1, K-length(relNums_without_mean))])));
% % Pxx2=pow2db(xf.^2/(length(relNums_without_mean)*2*pi))';
% % % removing the trend
% % relNums_without_trend=detrend(relNums);
% % xf=abs(fftshift(fft([relNums_without_trend' zeros(1, K-length(relNums_without_trend))])));
% % Pxx3=pow2db(xf.^2/(length(relNums_without_trend)*2*pi))';
% % % log then subtract mean
% % relNums=relNums+0.001;
% % relNums_with_log=log(relNums);
% % relNums_with_log_without_mean=relNums_with_log-mean(relNums_with_log);
% % xf=abs(fftshift(fft([relNums_with_log_without_mean' zeros(1, K-length(relNums_with_log_without_mean))])));
% % Pxx4=pow2db(xf.^2/(length(relNums_with_log_without_mean)*2*pi))';
% % fs=-1:2/K:1-1/K;
% % 
% % % plot periodogram
% % figure(1)
% % plot(fs,Pxx1,fs,Pxx2,fs,Pxx3,fs,Pxx4,'LineWidth',2);
% % set(gca,'fontsize',20);
% % axis([0 1 -40 70])
% % title('Peridogram of Sunspot Series','FontSize',20);
% % xlabel('Normalised Frequency (x \pi rad/sample)', 'FontSize', 20);
% % ylabel('Power/Frequency (dB/rad/sample)', 'FontSize', 20);
% % legend('\fontsize{25}No Preprocessing','\fontsize{25}Removing Mean','\fontsize{25}Removing Trend','\fontsize{25}Log, Removing Mean')
% % % graph_saving('../report/images/part1/sunspot_periodogram');

%%---------------------------------------

% load sunspot.dat
% % (zero-padded) signal length 
% K = 2.^12;
% x = sunspot(:, 2);
% t = sunspot(:, 1);
% N = length(x);
% z = detrend(x-mean(x));
% l_ = log(x + eps);
% l = l_ - mean(l_);
% 
% % plot time series
% plot(t, x,'LineWidth',2,'Color','red');
% hold on;
% g = plot(t, z,'LineWidth',2,'Color','blue');
% hold on;
% plot(t, l,'LineWidth',2,'Color','green');
% title("Sunpots Time Series",'FontSize', 20);
% xlabel("Time (years)",'FontSize', 20);
% ylabel("Number of Sunspots",'FontSize', 20);
% legend(["Sunspots (Raw)", "Mean with detrend", "Mean with log"]);
% %ylim([-20 60]);
% set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
% save_figure('../Images/Section_1/Section_1_2_A1_Sunspot');
% 
% % Windowing used
% windows = [ones(N, 1), chebwin(N)];
% 
% for i = 1:size(windows, 2)
%     
%     h = windows(:, i);
%     label = labels(i, 1);
%     Y = abs(fftshift(fft([(x .* h)' zeros(1, K-N)])));
%     Power_yy = pow2db(Y.^2 / (N * 2 * pi))';
%     Power_yy_one = Power_yy((length(Power_yy)/2 + 1):end, 1);
%     Z = abs(fftshift(fft([(z .* h)' zeros(1, K-N)])));
%     Power_zz = pow2db(Z.^2 / (N * 2 * pi))';
%     % Periodogram one-side
%     Power_zz_one = Power_zz((length(Power_zz)/2 + 1):end, 1);
%     L = abs(fftshift(fft([(l .* h)' zeros(1, K-N)])));
%     PSD = pow2db(L.^2 / (N * 2 * pi))';
%     PSD_one = PSD((length(PSD)/2 + 1):end, 1);
%     f = 0:2/K:1-1/K;
%  
%     fig = figure("Name", sprintf("Sunspots Periodogram %s", label));
%     plot(f, Power_yy_one,'LineWidth',2,'Color','red');
%     hold on;
%     g = plot(f, Power_zz_one,'LineWidth',2,'Color','blue');
%     hold on;
%     plot(f, PSD_one,'LineWidth',2,'Color','green');
%     title('Sunpots Periodogram (Chebyshev Window)','FontSize', 20);
%     xlabel("Normalised Frequency",'FontSize', 20);
%     ylabel("Power Density (dB)",'FontSize', 20);
%     legend(["Sunspots (Raw)", "Mean with detrend", "Mean with log"]);
%     ylim([-20 60]);
%     set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
% end

% %-----------------------------------------------------------------------
% % Question B
% %-----------------------------------------------------------------------
% 
% eeg = load('Data/EEG_Data/EEG_Data/EEG_Data_Assignment2.mat');
% Samp_freq = eeg.fs;
% N = length(eeg.POz);
% K = Samp_freq * 10;
% time_diff = [10 5 1];
% 
% % mean removal
% POz = eeg.POz - mean(eeg.POz);
% % Standard Periodogram
% [Pxx, wx] = periodogram(POz, rectwin(N), K, Samp_freq, 'onesided');
% 
% fig1 = figure("Name", "EEG Standard Periodogram");
% plot(wx, pow2db(Pxx),'LineWidth',2,'Color','red');
% title("EEG Periodogram (Rectangular Window) ",'FontSize', 20);
% set(gca,'fontsize',20);
% xlabel("Frequency (Hz)",'FontSize', 20);
% ylabel("Power Density (dB)",'FontSize', 20);
% xticks(0:10:60);
% set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);axis([0 60 -150 -100]);
% 
% % Bartlett Averaged Periodogram
% fig2 = figure("Name", "EEG Standard Periodogram");
% for i = 1:length(time_diff)
% 
%     % segment samples length
%     L = Samp_freq * time_diff(i);
%     % periodogram
%     [Pzz, wz]= pwelch(POz, rectwin(L), 0, K, Samp_freq, 'onesided');
%     % figure
%     plot(wx, pow2db(Pzz), "DisplayName", sprintf("Time Diff = %d s", time_diff(i)),'LineWidth',2);
%     hold on;
%     
% end
% title("EEG Periodogram (Bartlett)",'FontSize', 20);
% set(gca,'fontsize',20);
% xlabel("Frequency (Hz)",'FontSize', 20);
% ylabel("Power Density (dB)",'FontSize', 20);
% xticks(0:10:60);
% set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);axis([0 60 -135 -90]);
% legend("show", "Location", "north");
% 
% for i = 1:length(time_diff)
%     fig3 = figure("Name", sprintf("EEG Standard vs Bartlett Periodogram time_diff=%d", time_diff(i)));
%     L = Samp_freq * time_diff(i);
%  
%     [Pzz, wz]= pwelch(POz, rectwin(L), 0, K, Samp_freq, 'onesided');
%     % figure Generation
%     plot(wx, pow2db(Pxx),'LineWidth',2,'Color','red', "DisplayName", "Standard");
%     hold on;
%     plot(wx, pow2db(Pzz),'LineWidth',2,'Color','blue', "DisplayName", sprintf("Time Diff = %d s", time_diff(i)));
%     title(sprintf("EEG Periodograms"),'FontSize', 20);
%     set(gca,'fontsize',20);
%     xlabel("Frequency (Hz)",'FontSize', 20);
%     ylabel("Power Density (dB)",'FontSize', 20);
%     xticks(0:10:60);
%     set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
%     axis([0 60 -135 -90]);
%     legend("show", "Location", "North");
%     
% end
%-----------------------------------------------------------------------

% 
% load('data/EEG_Data/EEG_Data/EEG_Data_Assignment2.mat');
% 
% % Length of DFT, based on Hint of 5 DFT Samples per Hz
% K=fs*10;
% 
% % remove the meqan
% POz=POz-mean(POz);
% 
% % standard periodogram
% [pxx,w]=periodogram(POz,rectwin(length(POz)),K,fs,'twosided');
% figure(1);
% plot(w,pow2db(pxx),'LineWidth',2);
% axis([0 60 -140 -100]);
% set(gca,'fontsize',20);
% title('EEG Periodogram','FontSize',20);
% xlabel('Frequency (Hz)', 'FontSize', 20);
% ylabel('Power/Frequency (dB/rad/sample)', 'FontSize', 20);
% set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
% save_figure('../Images/Section_1/Section_1_2_B1_EEG_Periodogram');
% 
% % Bartlett ---------------------------------------
% 
% % Bartlett Averaged Periodogram - Rectwindow, No Overlap
% load('data/EEG_Data/EEG_Data/EEG_Data_Assignment2.mat');
% 
% % Length of DFT, based on Hint of 5 DFT Samples per Hz
% K=fs*10;
% 
% % remove the mean
% POz=POz-mean(POz);
% 
% % averaged over 10s
% segmentLength = 12000;
% noverlap = 0;
% [pxx,w]= pwelch(POz,rectwin(segmentLength),noverlap,K,fs,'twosided');
% figure(2);
% plot(w,pow2db(pxx),'LineWidth',2,'Color','red');
% hold on;
% 
% % averaged over 5s
% segmentLength = 6000;
% noverlap = 0;
% [pxx,w]= pwelch(POz,rectwin(segmentLength),noverlap,K,fs,'twosided');
% plot(w,pow2db(pxx),'LineWidth',2,'Color','blue');
% 
% % averaged over 1s
% segmentLength = 1200;
% noverlap = 0;
% [pxx,w]= pwelch(POz,rectwin(segmentLength),noverlap,K,fs,'twosided');
% plot(w,pow2db(pxx),'LineWidth',2,'Color','green');
% hold off;
% 
% axis([0 60 -140 -100]);
% set(gca,'fontsize',20);
% title('Bartlett Periodogram of EEG','FontSize',20);
% xlabel('Frequency (Hz)', 'FontSize', 20);
% ylabel('Power/Frequency (dB/rad/sample)', 'FontSize', 20);
% legend('\fontsize{25}Window Length: 10s','\fontsize{25}Window Length: 5s','\fontsize{25}Window Length: 1s');
% % graph_saving('../report/images/part1/POz_welch_periodogram');
% 
% 
% % Bartlett Periodogram only 1 second
% load('data/EEG_Data/EEG_Data/EEG_Data_Assignment2.mat');
% 
% % Length of DFT, based on Hint of 5 DFT Samples per Hz
% K=fs*10;
% 
% % center th mean
% POz=POz-mean(POz);
% 
% % averaged over 1s
% segmentLength = 1200;
% noverlap = 0;
% [pxx,w]= pwelch(POz,rectwin(segmentLength),noverlap,K,fs,'twosided');
% plot(w,pow2db(pxx),'LineWidth',2,'Color',[0.4940 0.1840 0.5560]);
% 
% axis([0 60 -140 -100]);
% set(gca,'fontsize',20);
% title('Bartlett Periodogram of EEG','FontSize',20);
% xlabel('Frequency (Hz)', 'FontSize', 20);
% ylabel('Power/Frequency (dB/rad/sample)', 'FontSize', 20);
% legend('\fontsize{25}Window Length: 1s');
% save_figure('../Images/Section_1/Section_1_2_B2_EEG_Periodogram');
