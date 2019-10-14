%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% Section 1.4
%-----------------------------------------------------------------------
%
%-----------------------------------------------------------------------
% Question A  
%-----------------------------------------------------------------------
% 
% rng(0)
% a=[2.76, -3.81, 2.65, -0.92];
% model=arima('Constant',0,'AR',a,'Variance',1);
% x=simulate(model,1000);
% x=x(501:end);
% AR_model=[4,8,10,14];
% legend_str{1, 1} = char('Ideal Spectrum'); 
% 
% for i = 1:length(AR_model)
%     [a_predicted, noise_var] = aryule(x, AR_model(i));
%     [magnitudes(:, i), ~] = freqz(noise_var^(1/2), a_predicted, length(x));
%     legend_str{1, i+1} = char(sprintf('Model - Based Method: Order %d', AR_model(i))); 
%     
% end
% [ideal, w] = freqz(1^(1/2), [1 -a], length(x));
% 
% When the model order is = 4----------------
% figure(1);
% plot(w/pi, abs(ideal).^2, 'LineWidth', 2,'color','red');
% hold on
% plot(w/pi, abs(magnitudes(:, 1)).^2, 'LineWidth',2,'color','green');
% hold off
% set(gca,'fontsize',20);
% title('AR Spectrum for n=500: Order=4','FontSize',20*32/24);
% xlabel('Normalised Freq', 'FontSize', 20);
% ylabel('Power/Freq', 'FontSize', 20);
% legend({'Ideal','Model','color'},'FontSize',15);
% set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
% save_figure('../Images/Section_1/Section_1_4_A1_AR_Specttrum_Est');
% 
% 
% When the model order is = 8----------------
% figure(2)
% plot(w/pi, abs(ideal).^2, 'LineWidth', 2,'color','red');
% hold on
% plot(w/pi, abs(magnitudes(:, 2)).^2, 'LineWidth', 2,'color','green');
% hold off
% set(gca,'fontsize',20);
% title('AR Spectrum for n=500: Order=8','FontSize',20*32/24);
% xlabel('Normalised Freq', 'FontSize', 20);
% ylabel('Power/Freq', 'FontSize', 20);
% legend({'Ideal','Model','color'},'FontSize',15);
% set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
% save_figure('../Images/Section_1/Section_1_4_A2_AR_Specttrum_Est');
% 
% When the model order is = 10----------------
% figure(3)
% plot(w/pi, abs(ideal).^2, 'LineWidth', 2,'color','red');
% hold on
% plot(w/pi, abs(magnitudes(:, 3)).^2, 'LineWidth', 2,'color','green');
% hold off
% set(gca,'fontsize',20);
% title('AR Spectrum for n=500: Order=10','FontSize',20*32/24);
% xlabel('Normalised Freq', 'FontSize', 20);
% ylabel('Power/Freq', 'FontSize', 20);
% legend({'Ideal','Model','color'},'FontSize',15);
% set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
% save_figure('../Images/Section_1/Section_1_4_A3_AR_Specttrum_Est');
% 
% When the model order is = 14-----------------
% figure(4)
% plot(w/pi, abs(ideal).^2, 'LineWidth', 2,'color','red');
% hold on
% plot(w/pi, abs(magnitudes(:, 4)).^2, 'LineWidth', 2,'color','green');
% hold off
% set(gca,'fontsize',20);
% title('AR Spectrum for n=500: Order=14','FontSize',20*32/24);
% xlabel('Normalised Freq', 'FontSize', 20);
% ylabel('Power/Freq', 'FontSize', 20);
% legend({'Ideal','Model','color'},'FontSize',15);
% set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
% save_figure('../Images/Section_1/Section_1_4_A4_AR_Specttrum_Est');

%-----------------------------------------------------------------------
%Question B
%-----------------------------------------------------------------------

% rng(0)
% a=[2.76, -3.81, 2.65, -0.92];
% model=arima('Constant',0,'AR',a,'Variance',1);
% x=simulate(model,10000);
% x=x(501:end);
% AR_model=[2,4,12];
% legend_str{1, 1} = char('Ideal Spectrum'); 
% 
% for i = 1:length(AR_model)
%     [a_predicted, noise_var] = aryule(x, AR_model(i));
%     [magnitudes(:, i), ~] = freqz(noise_var^(1/2), a_predicted, length(x));
%     legend_str{1, i+1} = char(sprintf('Model - Based Method: Order %d', AR_model(i))); 
%     
% end
% [ideal, w] = freqz(1^(1/2), [1 -a], length(x));
% 
% When the model order is = 2----------------
% figure(1);
% plot(w/pi, abs(ideal).^2, 'LineWidth', 2,'color','red');
% hold on
% plot(w/pi, abs(magnitudes(:, 1)).^2, 'LineWidth', 2,'color','green');
% hold off
% set(gca,'fontsize',20);
% title('AR Spectrum for n=9500: Order=2','FontSize',20);
% xlabel('Normalised Freq', 'FontSize', 20);
% ylabel('Power/Freq', 'FontSize', 20);
% legend({'Ideal','Model','color'},'FontSize',15);
% set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
% save_figure('../Images/Section_1/Section_1_4_B1_AR_Specttrum_Est');
% 
% When the model order is = 4----------------
% figure(2)
% plot(w/pi, abs(ideal).^2, 'LineWidth', 2,'color','red');
% hold on
% plot(w/pi, abs(magnitudes(:, 2)).^2, 'LineWidth', 2,'color','green');
% hold off
% axis([0 1 0 40000])
% set(gca,'fontsize',20);
% title('AR Spectrum for n=9500: Order=4','FontSize',20);
% xlabel('Normalised Freq', 'FontSize', 20);
% ylabel('Power/Freq', 'FontSize', 20);
% legend({'Ideal','Model','color'},'FontSize',15);
% set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
% save_figure('../Images/Section_1/Section_1_4_B2_AR_Specttrum_Est');
% 
% When the model order is = 12----------------
% figure(3)
% plot(w/pi, abs(ideal).^2, 'LineWidth', 2,'color','red');
% hold on
% plot(w/pi, abs(magnitudes(:, 3)).^2, 'LineWidth', 2,'color','green');
% hold off
% axis([0 1 0 40000])
% set(gca,'fontsize',20);
% title('AR Spectrum for n=9500: Order=12','FontSize',20);
% xlabel('Normalised Freq', 'FontSize', 20);
% ylabel('Power/Freq', 'FontSize', 20);
% legend({'Ideal','Model','color'},'FontSize',15);
% set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
% save_figure('../Images/Section_1/Section_1_4_B3_AR_Specttrum_Est');

%-----------------------------------------------------------------------
%Question C
%-----------------------------------------------------------------------
load('Data/PCR/PCAPCR.mat');
[~,S,~]=svd(X);
S=diag(S(1:10,:));
[~,Snoise_signal,~]=svd(Xnoise);
Snoise_signal=diag(Snoise_signal(1:10,:));

%Plotting X---------------------
figure(1)
plot(S,'LineWidth', 2,'color','red');
axis([1 10 0 40]);
set(gca,'fontsize',20);
title('X - Singular Values','FontSize',20);
xlabel('Index Number', 'FontSize', 20);
ylabel('Magnitude of Single Value', 'FontSize', 20);
set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
%save_figure('../Images/Section_1/Section_1_4_C1_SVD_PCR');

% Ploting Xnoise---------------------

figure(2);
plot(Snoise_signal,'LineWidth', 2,'color','red');
axis([1 10 0 40]);
set(gca,'fontsize',20);
title('X with Noise - Singular Values','FontSize',20);
xlabel('Index Number', 'FontSize', 20);
ylabel('Magnitude of Single Value', 'FontSize', 20);
set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
%save_figure('../Images/Section_1/Section_1_4_C2_SVD_PCR');

% Plotting Error in Values---------------------
figure(3);
error=(S-Snoise_signal).^2;
plot(error,'LineWidth', 2,'color','red');
axis([1 10 0 300]);
set(gca,'fontsize',20);
title('Errors in Singular Values','FontSize',20);
xlabel('Index Number', 'FontSize', 20);
ylabel('Magnitude of Single Value', 'FontSize', 20);
set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
%save_figure('../Images/Section_1/Section_1_4_C3_SVD_PCR');
