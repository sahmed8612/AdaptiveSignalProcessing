%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% Section 1.5
%-----------------------------------------------------------------------
%
%-----------------------------------------------------------------------
% Question A  
%-----------------------------------------------------------------------
% RRI Trial 1 data---------------------------
% load('Data/trial1.mat');
% load('Data/RRI_fs.mat');
% x=detrend(xRRI1);
% x=x.*hann(length(x))';
% N=length(x);
% K=2048;
% data_f=abs(fftshift(fft([x zeros(1, K-N)])));
% Power_xx =pow2db(data_f.^2/(N*2*pi));
% fs=-1:2/K:1-1/K;
% 
% figure(1)
% plot(fs*RRI_fs./2,Power_xx,'LineWidth',2,'color','red');
% axis([0 2 -100 0])
% set(gca,'fontsize',20);
% title('Trial 1 - Periodogram','FontSize',20);
% xlabel('Frequency (Hz)', 'FontSize', 20);
% ylabel('Power/Freq (dB/Hz)', 'FontSize', 20);
% set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
% save_figure('../Images/Section_1/Section_1_5_A1_Periodgram_RRI');
% 
% % RRI Trial 2 data---------------------------
% load('Data/trial2.mat');
% load('Data/RRI_fs.mat');
% x=detrend(xRRI2);
% x=x.*hann(length(x))';
% K=2048;
% N=length(x);
% data_f=abs(fftshift(fft([x zeros(1, K-N)])));
% Power_xx =pow2db(data_f.^2/(N*2*pi));
% fs=-1:2/K:1-1/K;
% 
% figure(2)
% plot(fs*RRI_fs./2,Power_xx,'LineWidth',2,'color','red');
% axis([0 2 -100 0])
% set(gca,'fontsize',20);
% title('Trial 2 - Periodogram','FontSize',20);
% xlabel('Frequency (Hz)', 'FontSize', 20);
% ylabel('Power/Freq (dB/Hz)', 'FontSize', 20);
% set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
% save_figure('../Images/Section_1/Section_1_5_A2_Periodgram_RRI');
% 
% % RRI Trial 3 data---------------------------
% load('Data/trial3.mat');
% load('Data/RRI_fs.mat');
% x=detrend(xRRI3);
% x=x.*hann(length(x))';
% N=length(x);
% K=2048;
% data_f=abs(fftshift(fft([x zeros(1, K-N)])));
% Power_xx =pow2db(data_f.^2/(N*2*pi));
% fs=-1:2/K:1-1/K;
% 
% figure(3)
% plot(fs*RRI_fs./2,Power_xx,'LineWidth',2,'color','red');
% axis([0 2 -100 0])
% set(gca,'fontsize',20);
% title('Trial 3 - Periodogram','FontSize',20);
% xlabel('Frequency (Hz)', 'FontSize', 20);
% ylabel('Power/Freq (dB/Hz)', 'FontSize', 20);
% set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
% save_figure('../Images/Section_1/Section_1_5_A3_Periodgram_RRI');

%-----------------------------------------------------------------------
% Question B
%-----------------------------------------------------------------------
%Trial 1----------------------
load('Data/trial1.mat');
load('Data/RRI_fs.mat');
K=2048;
x=detrend(xRRI1);
% 200 Samples Averaged----------------
Length_seg = 200;
length_overlap = 0;
[Power_xx,w]= pwelch(x,hann(Length_seg),length_overlap,K,fsRRI,'onesided');
figure(1);
plot(w,pow2db(Power_xx),'LineWidth',2,'Color','red');
hold on;
% 100 Samples Averaged----------------
Length_seg = 100;
length_overlap = 0;
[Power_xx,w]= pwelch(x,rectwin(Length_seg),length_overlap,K,fsRRI,'onesided');
plot(w,pow2db(Power_xx),'LineWidth',2,'Color','blue');
% 50 Samples Averaged----------------
Length_seg = 50;
length_overlap = 0;
[Power_xx,w]= pwelch(x,rectwin(Length_seg),length_overlap,K,fsRRI,'onesided');
plot(w,pow2db(Power_xx),'LineWidth',2,'Color','green');
hold off;

axis([0 2 -100 0]);
set(gca,'fontsize',20);
title('Trial 1 - Averaged Periodogram','FontSize',20);
xlabel('Frequecy (Hz)', 'FontSize', 20);
ylabel('Power/Freq (dB/rad/sample)', 'FontSize', 20);
legend({'200 Samples','100 Samples','50 Samples','color'},'FontSize', 15);
set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
save_figure('../Images/Section_1/Section_1_5_B1_Avg_Periodogram_RRI_SS');

% Trial 2----------------------
load('Data/trial2.mat');
load('Data/RRI_fs.mat');
x=detrend(xRRI2);
% 200 Samples Averaged----------------
Length_seg = 200;
length_overlap = 0;
[Power_xx,w]= pwelch(x,hann(Length_seg),length_overlap,K,fsRRI,'onesided');
figure(2);
plot(w,pow2db(Power_xx),'LineWidth',2,'Color','red');
hold on;
% 100 Samples Averaged----------------
Length_seg = 100;
length_overlap = 0;
[Power_xx,w]= pwelch(x,rectwin(Length_seg),length_overlap,K,fsRRI,'onesided');
plot(w,pow2db(Power_xx),'LineWidth',2,'Color','blue');
% 50 Samples Averaged----------------
Length_seg = 50;
length_overlap = 0;
[Power_xx,w]= pwelch(x,rectwin(Length_seg),length_overlap,K,fsRRI,'onesided');
plot(w,pow2db(Power_xx),'LineWidth',2,'Color','green');
hold off;

axis([0 2 -100 0]);
set(gca,'fontsize',20);
title('Trial 2 - Averaged Periodogram','FontSize',20);
xlabel('Frequecy (Hz)', 'FontSize', 20);
ylabel('Power/Freq (dB/rad/sample)', 'FontSize', 20);
legend({'200 Samples','100 Samples','50 Samples','color'},'FontSize', 15);
set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
save_figure('../Images/Section_1/Section_1_5_B2_Avg_Periodogram_RRI_SS');
% Trial 3----------------------

load('Data/trial3.mat');
load('Data/RRI_fs.mat');
x=detrend(xRRI3);
% 200 Samples Averaged----------------
Length_seg = 200;
length_overlap = 0;
[Power_xx,w]= pwelch(x,hann(Length_seg),length_overlap,K,fsRRI,'onesided');
figure(3);
plot(w,pow2db(Power_xx),'LineWidth',2,'Color','red');
hold on;
% 100 Samples Averaged----------------
Length_seg = 100;
length_overlap = 0;
[Power_xx,w]= pwelch(x,rectwin(Length_seg),length_overlap,K,fsRRI,'onesided');
plot(w,pow2db(Power_xx),'LineWidth',2,'Color','blue');

% 50 Samples Averaged----------------
Length_seg = 50;
length_overlap = 0;
[Power_xx,w]= pwelch(x,rectwin(Length_seg),length_overlap,K,fsRRI,'onesided');
plot(w,pow2db(Power_xx),'LineWidth',2,'Color','green');
hold off;

axis([0 2 -100 0]);
set(gca,'fontsize',20);
title('Trial 3 - Averaged Periodogram','FontSize',20);
xlabel('Frequecy (Hz)', 'FontSize', 20);
ylabel('Power/Freq (dB/rad/sample)', 'FontSize', 20);
legend({'200 Samples','100 Samples','50 Samples','color'},'FontSize', 15);
set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
save_figure('../Images/Section_1/Section_1_5_B3_Avg_Periodogram_RRI_SS');

%-----------------------------------------------------------------------
% Question C
%-----------------------------------------------------------------------
% % Trial 1-----------------------------------
% clear magnitudes;
% load('Data/trial1.mat');
% load('Data/RRI_fs.mat');
% 
% models=[2,4,10];
% x=detrend(xRRI1);
% x=x.*hann(length(x))';
% N=length(x);
% K=2048;
% signal_f=abs(fftshift(fft([x zeros(1, K-N)])));
% Power_xx=pow2db(signal_f.^2/(N*2*pi));
% fs=-1:2/K:1-1/K;
% 
% for i = 1:length(models)
%     [a_predicted, noise_var] = aryule(x, models(i));
%     [magnitudes(:, i), w] = freqz(noise_var^(1/2), a_predicted, length(x),fsRRI);
%     legend_str{1, i+1} = char(sprintf('Model - Based Method: Order %d', models(i))); 
%     
% end
% 
% % When Model order = 2
% figure(1)
% plot(fs*fsRRI./2,Power_xx,'LineWidth',2,'Color','red');
% hold on;
% plot(w, pow2db(abs(magnitudes(:, 1)).^2), 'LineWidth', 2,'Color','blue');
% legend({'Standard', 'Model','color'}, 'FontSize', 15);
% hold off;
% axis([0 2 -100 0]);
% set(gca,'fontsize',20);
% title('AR Spectrum Trial 1: Order=2','FontSize',20);
% xlabel('Frequency (Hz)', 'FontSize', 20);
% ylabel('Power/Freq (dB/Hz)', 'FontSize', 20);
% set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
% save_figure('../Images/Section_1/Section_1_5_C1_AR_Specttrum_RRI');
% 
% % When Model order = 4
% figure(2)
% plot(fs*fsRRI./2,Power_xx,'LineWidth',2,'Color','red');
% hold on;
% plot(w, pow2db(abs(magnitudes(:, 2)).^2), 'LineWidth', 2,'Color','blue');
% legend({'Standard', 'Model','color'}, 'FontSize', 15);
% hold off;
% axis([0 2 -100 0]);
% set(gca,'fontsize',20);
% title('AR Spectrum Trial 1: Order=4','FontSize',20);
% xlabel('Frequency (Hz)', 'FontSize', 20);
% ylabel('Power/Freq (dB/Hz)', 'FontSize', 20);
% set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
% save_figure('../Images/Section_1/Section_1_5_C2_AR_Specttrum_RRI');
% 
% % When Model order = 10
% figure(3)
% plot(fs*fsRRI./2,Power_xx,'LineWidth',2,'Color','red');
% hold on;
% plot(w, pow2db(abs(magnitudes(:, 3)).^2), 'LineWidth', 2, 'Color','blue');
% legend({'Standard', 'Model','color'}, 'FontSize', 15);
% hold off;
% axis([0 2 -100 0]);
% set(gca,'fontsize',20);
% title('AR Spectrum Trial 1: Order=10','FontSize',20);
% xlabel('Frequency (Hz)', 'FontSize', 20);
% ylabel('Power/Freq (dB/Hz)', 'FontSize', 20);
% set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
% save_figure('../Images/Section_1/Section_1_5_C3_AR_Specttrum_RRI');
% 
% 
% % Trial 2-----------------------------------
% clear magnitudes;
% load('Data/trial2.mat');
% load('Data/RRI_fs.mat');
% 
% models=[8,12,20];
% x=detrend(xRRI2);
% 
% x=x.*hann(length(x))';
% N=length(x);
% K=2048;
% signal_f=abs(fftshift(fft([x zeros(1, K-N)])));
% Power_xx=pow2db(signal_f.^2/(N*2*pi));
% fs=-1:2/K:1-1/K;
% 
% for i = 1:length(models)
%     % use aryule function to determine coefficients as well as estimated error
%     [a_predicted, noise_var] = aryule(x, models(i));
% 
%     % generate model based estimate of PSD
%     [magnitudes(:, i), w] = freqz(noise_var^(1/2), a_predicted, length(x),fsRRI);
%     
%     % create string for legend
%     legend_str{1, i+1} = char(sprintf('Model - Based Method: Order %d', models(i))); 
%     
% end
% 
% % When Model order = 8
% 
% figure(1)
% plot(fs*fsRRI./2,Power_xx,'LineWidth',2,'Color', 'red');
% hold on;
% plot(w, pow2db(abs(magnitudes(:, 1)).^2), 'LineWidth', 2,'Color', 'blue');
% legend({'Standard', 'Model','color'}, 'FontSize', 15);
% hold off;
% axis([0 2 -100 0]);
% set(gca,'fontsize',20);
% title('AR Spectrum Trial 2: Order=8','FontSize',20);
% xlabel('Frequency (Hz)', 'FontSize', 20);
% ylabel('Power/Freq (dB/Hz)', 'FontSize', 20);
% set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
% save_figure('../Images/Section_1/Section_1_5_C4_AR_Specttrum_RRI');
% % When Model orde = 12
% 
% figure(2)
% plot(fs*fsRRI./2,Power_xx,'LineWidth',2,'Color', 'red');
% hold on;
% plot(w, pow2db(abs(magnitudes(:, 2)).^2), 'LineWidth', 2, 'Color', 'blue');
% legend({'Standard', 'Model','color'}, 'FontSize', 15);
% hold off;
% axis([0 2 -100 0]);
% set(gca,'fontsize',20);
% title('AR Spectrum Trial 2: Order=12','FontSize',20);
% xlabel('Frequency (Hz)', 'FontSize', 20);
% ylabel('Power/Freq (dB/Hz)', 'FontSize', 20);
% set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
% save_figure('../Images/Section_1/Section_1_5_C5_AR_Specttrum_RRI');
% 
% %  When Model order = 20
% 
% figure(3)
% plot(fs*fsRRI./2,Power_xx,'LineWidth',2,'Color', 'red');
% hold on;
% plot(w, pow2db(abs(magnitudes(:, 3)).^2), 'LineWidth', 2, 'Color','blue');
% legend({'Standard', 'Model','color'}, 'FontSize', 15);
% hold off;
% axis([0 2 -100 0]);
% set(gca,'fontsize',20);
% title('AR Spectrum Trial 2: Order=20','FontSize',20);
% xlabel('Frequency (Hz)', 'FontSize', 20);
% ylabel('Power/Freq (dB/Hz)', 'FontSize', 20);
% set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
% save_figure('../Images/Section_1/Section_1_5_C6_AR_Specttrum_RRI');
% 
% % Trial 3-----------------------------------
% clear magnitudes;
% load('Data/trial3.mat');
% load('Data/RRI_fs.mat');
% 
% models=[2,13,35];
% x=detrend(xRRI3);
% 
% % periodogram plot 
% x=x.*hann(length(x))';
% N=length(x);
% K=2048;
% signal_f=abs(fftshift(fft([x zeros(1, K-N)])));
% Power_xx=pow2db(signal_f.^2/(N*2*pi));
% fs=-1:2/K:1-1/K;
% 
% for i = 1:length(models)
%     % use aryule function to determine coefficients as well as estimated error
%     [a_predicted, noise_var] = aryule(x, models(i));
% 
%     % generate model based estimate of PSD
%     [magnitudes(:, i), w] = freqz(noise_var^(1/2), a_predicted, length(x),fsRRI);
%     
%     % create string for legend
%     legend_str{1, i+1} = char(sprintf('Model - Based Method: Order %d', models(i))); 
%     
% end
% 
% % When Model orde = 2
% 
% figure(1)
% plot(fs*fsRRI./2,Power_xx,'LineWidth',2,'Color', 'red');
% hold on;
% plot(w, pow2db(abs(magnitudes(:, 1)).^2), 'LineWidth', 2,'Color', 'blue');
% legend({'Standard', 'Model','color'}, 'FontSize', 15);
% hold off;
% axis([0 2 -100 0]);
% set(gca,'fontsize',20);
% title('AR Spectrum Trial 3: Order=2','FontSize',20);
% xlabel('Frequency (Hz)', 'FontSize', 20);
% ylabel('Power/Freq (dB/Hz)', 'FontSize', 20);
% set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
% save_figure('../Images/Section_1/Section_1_5_C7_AR_Specttrum_RRI');
% 
% % When Model orde = 13
% 
% figure(2)
% plot(fs*fsRRI./2,Power_xx,'LineWidth',2,'Color', 'red');
% hold on;
% plot(w, pow2db(abs(magnitudes(:, 2)).^2), 'LineWidth', 2, 'Color', 'blue');
% legend({'Standard', 'Model','color'}, 'FontSize', 15);
% hold off;
% axis([0 2 -100 0]);
% set(gca,'fontsize',20);
% title('AR Spectrum Trial 3: Order=13','FontSize',20);
% xlabel('Frequency (Hz)', 'FontSize', 20);
% ylabel('Power/Freq (dB/Hz)', 'FontSize', 20);
% set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
% save_figure('../Images/Section_1/Section_1_5_C8_AR_Specttrum_RRI');
% 
% % When Model orde = 35
% 
% figure(3)
% plot(fs*fsRRI./2,Power_xx,'LineWidth',2,'Color', 'red');
% hold on;
% plot(w, pow2db(abs(magnitudes(:, 3)).^2), 'LineWidth', 2, 'Color', 'blue');
% legend({'Standard', 'Model','color'}, 'FontSize', 15);
% hold off;
% axis([0 2 -100 0]);
% set(gca,'fontsize',20);
% title('AR Spectrum Trial 3: Order=35','FontSize',20);
% xlabel('Frequency (Hz)', 'FontSize', 20);
% ylabel('Power/Freq (dB/Hz)', 'FontSize', 20);
% set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
% save_figure('../Images/Section_1/Section_1_5_C9_AR_Specttrum_RRI');
