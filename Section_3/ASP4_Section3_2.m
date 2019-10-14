%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% Section 3.2
%-----------------------------------------------------------------------
%
%-----------------------------------------------------------------------
% Question A  
%-----------------------------------------------------------------------

rng(0);
load Data/colors.mat
order=1;
K=1024;
N=1500;
freq_samp=1600;
noise_power=0.05;
f=zeros(N,1);
f(1:500)=100;
f(501:1000)=100+(1:500)./2;
f(1001:1500)=100+((1:500)./25).^2;
phase=cumsum(f);

% FM signal creation
x=exp(1j*(phase*2*pi/freq_samp));
w=wgn(N,1,pow2db(0.05),'complex');
y=x+w;
a=aryule(y,order);

% Block-Based Estimate from AR Coefficients
[h_aryule, w_aryule]=freqz(1, a, N, freq_samp);

% AR Coefficients estimated from Adaptive filter
mu=0.1;
[~,~,h_clms_mu_1] = ar_clms(y, mu, order);
mu=0.01;
[~,~,h_clms_mu_01] = ar_clms(y, mu, order);
mu=0.0025;
[~,~,h_clms_mu_0025] = ar_clms(y, mu, order);

H_mu_1=zeros(K,N,length(mu));
H_mu_01=zeros(K,N,length(mu));
H_mu_0025=zeros(K,N,length(mu));


% calculate the power spectrum for each time instance
for n = 1:N
    [h ,~]= freqz(1 , [1; -conj(h_clms_mu_1(n))], K, freq_samp);
    H_mu_1(:, n) = abs(h).^2; 
    [h ,~]= freqz(1 , [1; -conj(h_clms_mu_01(n))], K, freq_samp);
    H_mu_01(:, n) = abs(h).^2; 
    [h ,w]= freqz(1 , [1; -conj(h_clms_mu_0025(n))], K, freq_samp);
    H_mu_0025(:, n) = abs(h).^2; 
end
H_median=50*median(median(H_mu_1));
H_mu_1(H_mu_1>H_median)=H_median;
H_median=50*median(median(H_mu_01));
H_mu_01(H_mu_01>H_median)=H_median;
H_median=50*median(median(H_mu_0025));
H_mu_0025(H_mu_0025>H_median)=H_median;


% Original Frequency
figure(1)
plot(f, 'LineWidth', 2,'color','red');
plot_func('FM Signal Frequency', 'n', 'Frequency (Hz)', [0 1500 0 600], 1);
set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
save_figure('../Images/Section_3/Section_3_2_A1_time_freq_est');

% ArYule Spectrum
figure(2)
plot(w_aryule,mag2db(abs(h_aryule)), 'LineWidth', 2,'color','red');
plot_func('Aryule Magnitude Spectrum', 'Frequency (Hz)', 'Magnitude (dB)', [0 800 -20 30], 1);
set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
save_figure('../Images/Section_3/Section_3_2_A2_time_freq_est');


%-----------------------------------------------------------------------
% Question B  
%-----------------------------------------------------------------------

% Whne mu=0.1
figure(3)
surf(1:N, w, H_mu_1, 'LineStyle', 'none');
view(2);
plot_func('Frequency Change When \mu=0.1', 'n', 'Frequency (Hz)', [0 1500 0 800], 1);
set(gca, 'Color', 'none');
set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
save_figure('../Images/Section_3/Section_3_2_A3_time_freq_est');

% When mu=0.01
figure(4)
surf(1:N, w, H_mu_01, 'LineStyle', 'none');
view(2);
plot_func('Frequency Change When \mu=0.01', 'n', 'Frequency (Hz)', [0 1500 0 800], 1);
set(gca, 'Color', 'none');
set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
save_figure('../Images/Section_3/Section_3_2_A4_time_freq_est');

% When mu=0.0025
figure(5)
surf(1:N, w, H_mu_0025, 'LineStyle', 'none');
view(2);
plot_func('Frequency Change When \mu=0.0025', 'n', 'Frequency (Hz)', [0 1500 0 800], 1);
set(gca, 'Color', 'none');
set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
save_figure('../Images/Section_3/Section_3_2_A5_time_freq_est');

