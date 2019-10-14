%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% Section 1.3
%-----------------------------------------------------------------------
%
%-----------------------------------------------------------------------
% Question A  
%-----------------------------------------------------------------------
% Biased & Unbiased Correlogram using WGN -------------------------------
% rng(0);
% samples=1024;
% power_xx=1;
% x_wgn=wgn(samples,1,power_xx);
% [acf_biased,~,power_xx_biased,~]=correlation_est(x_wgn,'biased');
% [acf_unbiased,time_diff,power_xx_unbiased,samp_f]=correlation_est(x_wgn,'unbiased');
% 
% figure(1)
% plot(time_diff,acf_unbiased,'LineWidth',2,'Color','red')
% hold
% plot(time_diff,acf_biased,'LineWidth',2,'Color','blue')
% set(gca,'fontsize',20);
% title('ACF Using WGN','FontSize',45,'FontWeight','bold', 'HorizontalAlignment', 'center')
% set(gca,'fontsize',20)
% xlabel('k','fontsize',20)
% ylabel('r[k]','fontsize',25)
% axis([time_diff(1)-1 time_diff(end)+1 -1 2]);
% set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
% legend({'Unbiased','Biased','color'},'Fontsize', 15);
% save_figure('../Images/Section_1/Section_1_3_A1_WGN_estimates');
% 
% figure(2)
% plot(samp_f,real(power_xx_unbiased),'LineWidth',2,'Color','red');
% hold
% plot(samp_f,real(power_xx_biased),'LineWidth',2,'Color','blue');
% set(gca,'fontsize',20);
% title('Correlogram Using WGN','FontSize',45,'FontWeight','bold', 'HorizontalAlignment', 'center')
% set(gca,'fontsize',20)
% xlabel('Normalised Freq','fontsize',20)
% ylabel('Power','fontsize',25)
% axis([-1 1 -3 5]);
% set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
% legend({'Unbiased','Biased','color'},'Fontsize', 15);
% save_figure('../Images/Section_1/Section_1_3_A2_WGN_correlogram');
% 
% 
% %Biased & Unbiased Correlogram Using a sinewave--------------------------------
% rng(0);
% samples=1024;
% x_sine=sin(2*pi*1:samples);
% [acf_biased,~,power_xx_biased,~]=correlation_est(x_sine,'biased');
% [acf_unbiased,time_diff,power_xx_unbiased,samp_f]=correlation_est(x_sine,'unbiased');
% 
% figure(3)
% plot(time_diff,acf_unbiased,'LineWidth',2,'Color','red')
% hold
% plot(time_diff,acf_biased,'LineWidth',2,'Color','blue')
% set(gca,'fontsize',20);
% title('ACF Using Sinewave','FontSize',45,'FontWeight','bold', 'HorizontalAlignment', 'center')
% set(gca,'fontsize',20)
% xlabel('k','fontsize',20)
% ylabel('r[k]','fontsize',25)
% axis([time_diff(1)-1 time_diff(end)+1 -1 1]);
% set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
% legend({'Unbiased','Biased','color'},'Fontsize', 15);
% save_figure('../Images/Section_1/Section_1_3_A3_Sinewave_estimates');
% 
% figure(4)
% plot(samp_f,real(power_xx_unbiased),'LineWidth',2,'Color','red');
% hold
% plot(samp_f,real(power_xx_biased),'LineWidth',2,'Color','blue');
% set(gca,'fontsize',20);
% title('Correlogram Using WGN','FontSize',45,'FontWeight','bold', 'HorizontalAlignment', 'center')
% set(gca,'fontsize',20)
% xlabel('Normalised Freq','fontsize',20)
% ylabel('Power','fontsize',25)
% axis([-1 1 -2 4]);
% set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
% legend({'Unbiased','Biased','color'},'Fontsize', 15);
% save_figure('../Images/Section_1/Section_1_3_A4_Sinewave_correlogram');
% 
% 
% %Biased & Unbiased Correlogram Using a Filtered Sinewave ----------------------------------
% rng(0);
% samples=1024;
% power_xx=1;
% 
% a = 1;
% b = [1/4 1/4 1/4 1/4];
% x_wgn=wgn(samples,1,power_xx);
% x_filt=filter(b,a,x_wgn);
% [acf_biased,~,power_xx_biased,~]=correlation_est(x_filt,'biased');
% [acf_unbiased,time_diff,power_xx_unbiased,samp_f]=correlation_est(x_filt,'unbiased');
% 
% figure(5)
% plot(time_diff,acf_unbiased,'LineWidth',2,'Color','red')
% hold
% plot(time_diff,acf_biased,'LineWidth',2,'Color','blue')
% set(gca,'fontsize',20);
% title('ACF Using Filtered WGN','FontSize',45,'FontWeight','bold', 'HorizontalAlignment', 'center')
% set(gca,'fontsize',20)
% xlabel('k','fontsize',20)
% ylabel('r[k]','fontsize',25)
% axis([time_diff(1)-1 time_diff(end)+1 -0.3 0.5]);
% set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
% legend({'Unbiased','Biased','color'},'Fontsize', 15);
% save_figure('../Images/Section_1/Section_1_3_A5_Filtered_WGN__estimates');
% 
% 
% figure(6)
% plot(samp_f,real(power_xx_unbiased),'LineWidth',2,'Color','red');
% hold
% plot(samp_f,real(power_xx_biased),'LineWidth',2,'Color','blue');
% set(gca,'fontsize',20);
% title('Correlogram Using Filtered WGN','FontSize',45,'FontWeight','bold', 'HorizontalAlignment', 'center')
% set(gca,'fontsize',20)
% xlabel('Normalised Freq','fontsize',20)
% ylabel('Power','fontsize',25)
% axis([-1 1 -2 4]);
% set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
% legend({'Unbiased','Biased','color'},'Fontsize', 15);
% save_figure('../Images/Section_1/Section_1_3_A6_Filtered_WGN_correlogram');

%-----------------------------------------------------------------------
% Question B
%-----------------------------------------------------------------------
% Random process with 2 sinewaves

% F_Samp=5;
% n = 0:1/(2*F_Samp):10;
% freq_0=0.4;
% freq_1=0.8;
% amp=1;
% var=1;
% iter=100;
% power_xx_biased=zeros(iter,length(n)*2-1);
% 
% figure(7);
% for i=1:iter
%     rng(i)
%     w = var.*randn(length(n),1);
%     y = amp.*sin(2*pi*n*freq_0)+amp.*sin(2*pi*n*freq_1)+w';
%     [~,~,power_xx_biased(i,:),fs]=correlation_est(y,'biased');
%     s=plot(fs*F_Samp,real(power_xx_biased(i,:)),'Color','red','LineWidth',2);
%     hold on;
% end
% power_xx_mean=mean(real(power_xx_biased));
% m=plot(fs*F_Samp,power_xx_mean,'Color','green','LineWidth', 2);
% hold off;
% set(gca,'fontsize',20);
% axis([fs(1)*F_Samp fs(end)*F_Samp 0 10]);
% title('PSD Estimate: Freq=0.4,0.8','FontSize',20);
% xlabel('Frequency', 'FontSize', 20);
% ylabel('Power', 'FontSize', 20);
% set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
% legend([s m],{'Random Signal','Total Mean of Siganls','color'},'FontSize',15);
% %save_figure('../Images/Section_1/Section_1_3_B1_Sinewave_correlogram_Var');
% 
% 
% figure(8);
% power_xx_std=std(real(power_xx_biased));
% m=plot(fs*F_Samp,power_xx_std,'Color','green','LineWidth', 2);
% hold off;
% set(gca,'fontsize',20);
% axis([fs(1)*F_Samp fs(end)*F_Samp 0 3]);
% title('PSD Estimate Std: Freq=0.4,0.8','FontSize',20);
% xlabel('Frequency', 'FontSize', 20);
% ylabel('Power', 'FontSize', 20);
% set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
% %save_figure('../Images/Section_1/Section_1_3_B2_Sinewave_correlogram_std');
% 
% %Random process with 2 sinewaves, f=1.2,2.4
% 
% F_Samp=5;
% n = 0:1/(2*F_Samp):10;
% freq_0=1.2;
% freq_1=2.4;
% amp=1;
% var=1;
% iter=100;
% power_xx_biased=zeros(iter,length(n)*2-1);
% 
% figure(9);
% for i=1:iter
%     rng(i)
%     w = var.*randn(length(n),1);
%     y = amp.*sin(2*pi*n*freq_0)+amp.*sin(2*pi*n*freq_1) + w';
%     [~,~,power_xx_biased(i,:),fs]=correlation_est(y,'biased');
%     s=plot(fs*F_Samp,real(power_xx_biased(i,:)),'Color','red','LineWidth',2);
%     hold on;
% end
% power_xx_mean=mean(real(power_xx_biased));
% m=plot(fs*F_Samp,power_xx_mean,'Color','green','LineWidth', 2);
% hold off;
% set(gca,'fontsize',20);
% axis([fs(1)*F_Samp fs(end)*F_Samp 0 10]);
% title('PSD Estimate: Freq=1.2,2.4','FontSize',20);
% xlabel('Frequency', 'FontSize', 20);
% ylabel('Power', 'FontSize', 20);
% set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
% legend([s m],{'Random Signal','Total Mean of Siganls','color'},'FontSize',15);
% %save_figure('../Images/Section_1/Section_1_3_B3_Sinewave_correlogram_Var_f2');
% 
% figure(10);
% power_xx_std=std(real(power_xx_biased));
% m=plot(fs*F_Samp,power_xx_std,'Color','green','LineWidth', 2);
% hold off;
% set(gca,'fontsize',20);
% axis([fs(1)*F_Samp fs(end)*F_Samp 0 3]);
% title('PSD Estimate Std: Freq=0.4,0.8','FontSize',20);
% xlabel('Frequency', 'FontSize', 20);
% ylabel('Power', 'FontSize', 20);
% set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
% %save_figure('../Images/Section_1/Section_1_3_B4_Sinewave_correlogram_std_f2');
% 
%-----------------------------------------------------------------------
% Question C
%-----------------------------------------------------------------------
% Ploting figures in dB
% F_Samp=5;
% n = 0:1/(2*F_Samp):10;
% freq_0=1.2;
% freq_1=2.4;
% amp=1;
% var=1;
% iter=100;
% power_xx_biased=zeros(iter,length(n)*2-1);
% 
% figure(11);
% for i=1:iter
%     rng(i)
%     w = var.*randn(length(n),1);
%     y = amp.*sin(2*pi*n*freq_0)+amp.*sin(2*pi*n*freq_1) + w';
%     [~,~,power_xx_biased(i,:),fs]=correlation_est(y,'biased');
%     s=plot(fs*F_Samp,pow2db(real(power_xx_biased(i,:))),'Color','red','LineWidth',2);
%     hold on;
% end
% power_xx_mean=mean(real(power_xx_biased));
% m=plot(fs*F_Samp,pow2db(power_xx_mean),'Color','green','LineWidth', 2);
% hold off;
% set(gca,'fontsize',20);
% axis([fs(1)*F_Samp fs(end)*F_Samp -40 20]);
% title('PSD Estimate(dB): Freq=1.2,2.4','FontSize',20);
% xlabel('Frequency', 'FontSize', 20);
% ylabel('Power', 'FontSize', 20);
% set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
% legend([s m],{'Random Signal','Total Mean of Siganls','color'},'FontSize',15);
% %save_figure('../Images/Section_1/Section_1_3_C1_Sinewave_correlogram_Var_db');
% 
% 
% figure(12);
% power_xx_std=std(real(power_xx_biased));
% m=plot(fs*F_Samp,pow2db(power_xx_std),'Color','green','LineWidth', 2);
% hold off;
% set(gca,'fontsize',20);
% axis([fs(1)*F_Samp fs(end)*F_Samp -15 5]);
% title('Std of dB plots','FontSize',20);
% xlabel('Frequency', 'FontSize', 20);
% ylabel('Power', 'FontSize', 20);
% set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
% %save_figure('../Images/Section_1/Section_1_3_C2_Std_of_db');
%
%-----------------------------------------------------------------------
% Question D
%-----------------------------------------------------------------------
% Calculating when n=30;

% K=128;
% n = 0:30;
% noise_signal = 0.2/sqrt(2)*(randn(size(n))+1j*randn(size(n)));
% x = exp(1j*2*pi*0.3*n)+exp(1j*2*pi*0.32*n)+ noise_signal;
% power_xx=fft([x zeros(1,K-length(x))])./(1*length(n));
% fs=0:1/K:1-1/K;
% figure(1);
% plot(fs,pow2db(abs(power_xx)),'LineWidth',2,'color','red');
% set(gca,'fontsize',20);
% axis([0 1 -30 5]);
% title('Complex exp Periodogram: n=30','FontSize',20);
% xlabel('Frequency (Hz)', 'FontSize', 20);
% ylabel('Power/Freq (dB/Hz)', 'FontSize', 20);
% set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
% save_figure('../Images/Section_1/Section_1_3_D1_Complex_Periodogram');
% 
% 
% % Calculating when n=50;
% 
% K=128;
% n = 0:50;
% noise_signal = 0.2/sqrt(2)*(randn(size(n))+1j*randn(size(n)));
% x = exp(1j*2*pi*0.3*n)+exp(1j*2*pi*0.32*n)+ noise_signal;
% power_xx=fft([x zeros(1,K-length(x))])./(1*length(n));
% fs=0:1/K:1-1/K;
% figure(2);
% plot(fs,pow2db(abs(power_xx)),'LineWidth',2,'color','red');
% set(gca,'fontsize',20);
% axis([0 1 -30 5]);
% title('Complex exp Periodogram: n=50','FontSize',20);
% xlabel('Frequency (Hz)', 'FontSize', 20);
% ylabel('Power/Freq (dB/Hz)', 'FontSize', 20);
% set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
% save_figure('../Images/Section_1/Section_1_3_D2_Complex_Periodogram');
% 
% % Calculating when n=70;
% 
% K=128;
% n = 0:70;
% noise_signal = 0.2/sqrt(2)*(randn(size(n))+1j*randn(size(n)));
% x = exp(1j*2*pi*0.3*n)+exp(1j*2*pi*0.32*n)+ noise_signal;
% power_xx=fft([x zeros(1,K-length(x))])./(1*length(n));
% fs=0:1/K:1-1/K;
% figure(3);
% plot(fs,pow2db(abs(power_xx)),'LineWidth',2,'color','red');
% set(gca,'fontsize',20);
% axis([0 1 -30 5]);
% title('Complex exp Periodogram: n=70','FontSize',20);
% xlabel('Frequency (Hz)', 'FontSize', 20);
% ylabel('Power/Freq (dB/Hz)', 'FontSize', 20);
% set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
% save_figure('../Images/Section_1/Section_1_3_D3_Complex_Periodogram');

%-----------------------------------------------------------------------
% Question E
%-----------------------------------------------------------------------
% When frequencies are f0=0.3,0.32, as n=30-----------------------------
rng(1);

n = 0:30;
noise = 0.2/sqrt(2)*(randn(size(n))+1j*randn(size(n)));
x = exp(1j*2*pi*0.3*n)+exp(1j*2*pi*0.32*n)+ noise;
[~,R] = corrmtx(x,14,'mod');
[S,F] = pmusic(R,2,[ ],1,'corr');
figure(1);
plot(F,S,'linewidth',2,'color','red'); set(gca,'xlim',[0.25 0.40]);
set(gca,'fontsize',20);
title('MUSIC Algorithm: n=30 & f=0.30,0.32','FontSize',20);
xlabel('Freq (Hz)', 'FontSize', 20);
ylabel('Psuedospectrum', 'FontSize', 20);
set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
save_figure('../Images/Section_1/Section_1_3_E1_MUSIC');


% Resoulution when f0=0.3,0.32, n=30
rng(1);
iter=100;
S=zeros(iter,256);
figure(2);
for i=1:iter
n = 0:30;
noise = 0.2/sqrt(2)*(randn(size(n))+1j*randn(size(n)));
x = exp(1j*2*pi*0.3*n)+exp(1j*2*pi*0.32*n)+ noise;
[~,R] = corrmtx(x,14,'mod');
[S(i,:),F] = pmusic(R,2,[ ],1,'corr');
s=plot(F,S(i,:),'linewidth',2,'Color','blue');
if(i==2)
    hold on;
end
end
m=plot(F,mean(S),'linewidth',2,'Color','red');
hold off;
set(gca,'xlim',[0.25 0.40]);
set(gca,'fontsize',20);
title('Resolution MUSIC: n=30 & f=0.30,0.32','FontSize',20);
xlabel('Freq (Hz)', 'FontSize', 20);
ylabel('Psuedospectrum', 'FontSize', 20);
legend([s m],{'Random Signal','Total Mean of Signals','color'},'FontSize',15);
set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
save_figure('../Images/Section_1/Section_1_3_E2_MUSIC');



% f0=0.30,0.32,0.35 and when n=30-------------------------------------------
rng(1);

n = 0:30;
noise = 0.2/sqrt(2)*(randn(size(n))+1j*randn(size(n)));
x = exp(1j*2*pi*0.3*n)+exp(1j*2*pi*0.32*n)+exp(1j*2*pi*0.35*n)+noise;

[~,R] = corrmtx(x,14,'mod');
[S,F] = pmusic(R,3,[ ],1,'corr');
figure(3);
plot(F,S,'linewidth',2,'color','red'); set(gca,'xlim',[0.25 0.40]);
set(gca,'fontsize',20);
title('MUSIC Algorithm: n=30 & f=0.30,0.32,0.35','FontSize',20);
xlabel('Freq (Hz)', 'FontSize', 20);
ylabel('Psuedospectrum', 'FontSize', 20);
set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
save_figure('../Images/Section_1/Section_1_3_E3_MUSIC');



% Resolution when f0=0.30, f1=0.32 f2=0.35, n=30
rng(1);

iter=100;
S=zeros(iter,256);
figure(4);
for i=1:iter
n = 0:30;
noise = 0.2/sqrt(2)*(randn(size(n))+1j*randn(size(n)));
x = exp(1j*2*pi*0.3*n)+exp(1j*2*pi*0.32*n)+exp(1j*2*pi*0.35*n)+noise;
[~,R] = corrmtx(x,14,'mod');
[S(i,:),F] = pmusic(R,3,[ ],1,'corr');
s=plot(F,S(i,:),'linewidth',2,'Color','blue');
if(i==2)
    hold on;
end
end
m=plot(F,mean(S),'linewidth',2,'Color','red');
hold off;
set(gca,'xlim',[0.25 0.40]);
set(gca,'fontsize',20);
title('Resolution MUSIC: n=30 & f=0.30,0.32,0.35','FontSize',20);
xlabel('Freq (Hz)', 'FontSize', 20);
ylabel('Psuedospectrum', 'FontSize', 20);
legend([s m],{'Random Signal','Total Mean of Signals','color'},'FontSize',15);
set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
save_figure('../Images/Section_1/Section_1_3_E4_MUSIC');


% f0=0.30, f1=0.32 f2=0.35, n=70-------------------------------------------
rng(1);

n = 0:70;
noise = 0.2/sqrt(2)*(randn(size(n))+1j*randn(size(n)));
x = exp(1j*2*pi*0.3*n)+exp(1j*2*pi*0.32*n)+exp(1j*2*pi*0.35*n)+noise;

[~,R] = corrmtx(x,14,'mod');
[S,F] = pmusic(R,3,[ ],1,'corr');
figure(5)
plot(F,S,'linewidth',2,'color','red'); set(gca,'xlim',[0.25 0.40]);
set(gca,'fontsize',20);
title('MUSIC Algorithm: n=70 & f=0.30,0.32,0.35','FontSize',20);
xlabel('Freq (Hz)', 'FontSize', 20);
ylabel('Psuedospectrum', 'FontSize', 20);
set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
save_figure('../Images/Section_1/Section_1_3_E5_MUSIC');

% Resolution when f0=0.30, f1=0.32 f2=0.35, n=70

iter=100;
S=zeros(iter,256);
figure(6);
for i=1:iter
n = 0:70;
noise = 0.2/sqrt(2)*(randn(size(n))+1j*randn(size(n)));
x = exp(1j*2*pi*0.3*n)+exp(1j*2*pi*0.32*n)+exp(1j*2*pi*0.35*n)+noise;
[~,R] = corrmtx(x,14,'mod');
[S(i,:),F] = pmusic(R,3,[ ],1,'corr');
s=plot(F,S(i,:),'linewidth',2,'Color','blue');
if(i==2)
    hold on;
end
end
m=plot(F,mean(S),'linewidth',2,'Color','red');
hold off;
set(gca,'xlim',[0.25 0.40]);
set(gca,'fontsize',20);
title('Resolution MUSIC : n=70 & f=0.30,0.32,0.35','FontSize',20);
xlabel('Freq (Hz)', 'FontSize', 20);
ylabel('Psuedospectrum', 'FontSize', 20);
legend([s m],{'Random Signal','Total Mean of Signals','color'},'FontSize',15);
set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
save_figure('../Images/Section_1/Section_1_3_E6_MUSIC');
