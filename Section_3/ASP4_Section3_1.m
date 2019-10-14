%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% Section 3.1
%-----------------------------------------------------------------------
%
%-----------------------------------------------------------------------
% Question A  
%-----------------------------------------------------------------------
% rng(0);
% load Data/colors.mat
% 
% mu=0.1;
% order=1;
% realisations=100;
% N=1000;
% noise_power=1;
% b=[1.5+1i, 2.5-0.5i];
% a=1;
% clms_error=complex(zeros(realisations,N));
% aclms_error=complex(zeros(realisations,N));
% w=complex(zeros(realisations,N));
% y=complex(zeros(realisations,N));
% 
% % run realisations
% for i=1:realisations
% 
%     % get complex circular wgn
%     w(i,:)=wgn(1,N,pow2db(noise_power),'complex');
%     y(i,:)=complex_func(b,a,w(i,:));
%     
%     [~, e] = clms(y(i,:), w(i,:), mu, order);
%     clms_error(i,:)=abs(e).^2;
%     [~, e] = a_clms(y(i,:), w(i,:), mu, order);
%     aclms_error(i,:)=abs(e).^2;
%     
% end
% x_axis=1:N;
% start=700;
% 
% clms_error_mean=mean(clms_error);
% clms_error_steady_state=mean(clms_error_mean(start:end));
% aclms_error_mean=mean(aclms_error);
% aclms_error_steady_state=mean(aclms_error_mean(start:end));
% w_mean=mean(w);
% y_mean=mean(y);
% 
% % Learning Curves 1
% figure(1)
% plot(x_axis,pow2db(clms_error_mean),'LineWidth',2,'color','red');
% hold on;
% plot(x_axis,pow2db(aclms_error_mean),'LineWidth',2,'color','blue');
% hold off;
% legend1=sprintf('Complex, SState = %.1f',pow2db(clms_error_steady_state));
% legend2=sprintf('Augmented Complex, SState = %.1f',pow2db(aclms_error_steady_state));
% plot_legend_func('CLMS & ACLMS Learning Curves ', 'n', 'Error Power (dB)', {legend1, legend2}, [0 1000 -350 100], 1);
% set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
% save_figure('../Images/Section_3/Section_3_1_A1_ALCMS_CLMS_Comparison');
% 
% % Learning Curves 2
% figure(2)
% plot(real(w_mean),imag(w_mean),'*','LineWidth',2,'Color','red');
% hold on;
% plot(real(y_mean),imag(y_mean),'*','LineWidth',2,'Color','blue');
% hold off;
% plot_legend_func('Noise and WLMA(1) Distribution', 'Re', 'Im', {'Circular WGN', 'Widely Linear Moving Average Process (1st Order)'}, [-1.1 1.1 -1.1 1.1], 1);
% set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
% save_figure('../Images/Section_3/Section_3_1_A2_ALCMS_CLMS_Comparison');

%----------------------------
%Aryule Figures
N = 1500;
t = 1:N;
samp_freq = 1500;

sigma_squared = 0.05;
error_ta = wgn(N, 1, pow2db(0.05), 'complex');
f = arrayfun(@input_frequency, t)';
phi = cumtrapz(f);
y = exp(1j * (2 * pi * phi / samp_freq)) + error_ta;

for p = [1 5 10]
    a = aryule(y, p);
    [h, w] = freqz(1, a, N, samp_freq);
    psd = mag2db(abs(h));

    % figure - AR(p)
    fig = figure("Name", sprintf("AR(%d)", p));
    plot(w, psd,'LineWidth',2)
    title(sprintf("AR(%d) using Aryule (Complete FM)", p),'FontSize', 20);
    xlabel("Frequency (Hz)",'FontSize', 20);
    ylabel("Power Spectral Density (dB)",'FontSize', 20);
    set(gca,'fontsize',20);
    grid on;
    set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
end

for i = 1:3

    for p = [1]
     
        % starting index
        start = 1 + (i-1)*500;
        a = aryule(y(start:start+499), p);
        [h, w] = freqz(1, a, N, samp_freq);
        psd = mag2db(abs(h));

        fig = figure("Name", sprintf("AR(%d)", p));
        plot(w, psd, "Color", 'red','LineWidth',2)
        title(sprintf("AR(%d) using Aryule (FM Parts)", p),'FontSize', 20);
        xlabel("Frequency (Hz)",'FontSize', 20);
        ylabel("Power Spectral Density (dB)",'FontSize', 20);
        set(gca,'fontsize',20);
        grid on;
        set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);

    end

end

%-----------------------------------------------------------------------
% Question B
%-----------------------------------------------------------------------
% % complex signal for all wind types
% load Data/colors.mat
% load Data/wind-dataset/low-wind.mat
% wind_low=complex(v_east,v_north);
% load Data/wind-dataset/medium-wind.mat
% wind_med=complex(v_east,v_north);
% load Data/wind-dataset/high-wind.mat
% wind_high=complex(v_east,v_north);
% max_lim=10;
% 
% % circularity Values generation
% wind_low_mean_real=mean(real(wind_low));
% wind_low_vert_x=[wind_low_mean_real wind_low_mean_real];
% wind_low_vert_y=[-max_lim max_lim];
% wind_low_mean_imag=mean(imag(wind_low));
% wind_low_hor_x=[-max_lim max_lim];
% wind_low_hor_y=[wind_low_mean_imag wind_low_mean_imag];
% wind_low_mean = abs(mean((wind_low).^2)/mean(abs(wind_low).^2));
% 
% wind_med_mean =abs( mean((wind_med).^2)/mean(abs(wind_med).^2));
% wind_med_mean_real=mean(real(wind_med));
% wind_med_vert_x=[wind_med_mean_real wind_med_mean_real];
% wind_med_vert_y=[-max_lim max_lim];
% wind_med_mean_imag=mean(imag(wind_med));
% wind_med_hor_x=[-max_lim max_lim];
% wind_med_hor_y=[wind_med_mean_imag wind_med_mean_imag];
% 
% wind_high_mean =abs( mean((wind_high).^2)/mean(abs(wind_high).^2));
% wind_high_mean_real=mean(real(wind_high));
% wind_high_vert_x=[wind_high_mean_real wind_high_mean_real];
% wind_high_vert_y=[-max_lim max_lim];
% wind_high_mean_imag=mean(imag(wind_high));
% wind_high_hor_x=[-max_lim max_lim];
% wind_high_hor_y=[wind_high_mean_imag wind_high_mean_imag];
% 
% % length of data
% N=length(wind_low);
% x_vert=[0 0]; y_vert=[-10 10];
% x_hor=[-10 10]; y_hor=[0 0];
% mu=[0.001 0.0025 0.005 0.0075 0.0001];
% order=1:30;
% mpse_clms=zeros(3,length(order),length(mu));
% mpse_aclms=zeros(3,length(order),length(mu));
% 
% for i=1:length(mu)
%     for j=1:length(order)
%         [~,e]=ar_clms(wind_low,mu(i),order(j));
%         mpse_clms(1,j,i)=sum(abs(e).^2)./N;
%         [~,e]=ar_clms(wind_med,mu(i),order(j));
%         mpse_clms(2,j,i)=sum(abs(e).^2)./N;
%         [~,e]=ar_clms(wind_high,mu(i),order(j));
%         mpse_clms(3,j,i)=sum(abs(e).^2)./N;
%         
%         [~,e]=ar_aclms(wind_low,mu(i),order(j));
%         mpse_aclms(1,j,i)=sum(abs(e).^2)./N;
%         [~,e]=ar_aclms(wind_med,mu(i),order(j));
%         mpse_aclms(2,j,i)=sum(abs(e).^2)./N;
%         [~,e]=ar_aclms(wind_high,mu(i),order(j));
%         mpse_aclms(3,j,i)=sum(abs(e).^2)./N;
%     end
% end
% 
% % Low Wind
% figure(1);
% plot(real(wind_low),imag(wind_low),'o','LineWidth', 2, 'color','red');
% hold on;
% plot([-5 5],[0 0],'--','LineWidth', 2, 'Color', 'green');
% plot([0 0],[-5 5],'--','LineWidth', 2, 'Color', 'green');
% plot(wind_low_hor_x,wind_low_hor_y,'--','LineWidth', 2, 'Color', 'blue');
% plot(wind_low_vert_x,wind_low_vert_y,'--','LineWidth', 2, 'Color', 'blue');
% hold off;
% str=sprintf('Circularity Value =%.3f (Low Wind)', wind_low_mean);
% plot_func(str, 'v_{east}[n]', 'v_{north}[n]', [-4 4 -4 4], 1);
% set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
% save_figure('../Images/Section_3/Section_3_1_B1_Wind');
% 
% % Medium Wind
% figure(2);
% plot(real(wind_med),imag(wind_med),'o','LineWidth', 2,'color','red');
% hold on;
% plot([-5 5],[0 0],'--','LineWidth', 2*2, 'Color', 'green');
% plot([0 0],[-5 5],'--','LineWidth', 2*2, 'Color', 'green');
% plot(wind_med_hor_x,wind_med_hor_y,'--','LineWidth', 2, 'Color', 'blue');
% plot(wind_med_vert_x,wind_med_vert_y,'--','LineWidth', 2, 'Color', 'blue');
% hold off;
% str=sprintf('Circularity Value =%.3f (Medium Wind)', wind_med_mean);
% plot_func(str, 'v_{east}[n]', 'v_{north}[n]', [-4 4 -4 4], 1);
% set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
% save_figure('../Images/Section_3/Section_3_1_B2_Wind');
% 
% % High Wind
% figure(3);
% plot(real(wind_high),imag(wind_high),'o','LineWidth', 2,'color','red');
% hold on;
% plot([-5 5],[0 0],'--','LineWidth', 2*2, 'Color', 'green');
% plot([0 0],[-5 5],'--','LineWidth', 2*2, 'Color', 'green');
% plot(wind_high_hor_x,wind_high_hor_y,'--','LineWidth', 2, 'Color', 'blue');
% plot(wind_high_vert_x,wind_high_vert_y,'--','LineWidth', 2, 'Color', 'blue');
% hold off;
% str=sprintf('Circularity Value =%.3f (High Wind)', wind_high_mean);
% plot_func(str, 'v_{east}[n]', 'v_{north}[n]', [-4 4 -4 4], 1);
% set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
% save_figure('../Images/Section_3/Section_3_1_B3_Wind');
% 
% % Low Wind (Filter Orders)
% figure(4)
% plot(order,pow2db(mpse_clms(1,:,5)),'-*','LineWidth',2, 'Color', 'red');
% hold on;
% plot(order,pow2db(mpse_aclms(1,:,5)),'-*','LineWidth',2, 'Color', 'green');
% hold off;
% plot_legend_func('\mu=0.0001 (Low Wind)', 'Order of Filter', 'MPSE (dB)', {'Complex', 'Augmented Complex'},[0 30 -16 -14], 1);
% set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
% save_figure('../Images/Section_3/Section_3_1_B4_Wind');
% 
% % Medium Wind (Filter Orders)
% figure(5)
% plot(order,pow2db(mpse_clms(2,:,3)),'-*','LineWidth',2, 'Color', 'red');
% hold on;
% plot(order,pow2db(mpse_aclms(2,:,3)),'-*','LineWidth',2, 'Color', 'green');
% hold off;
% plot_legend_func('\mu=0.0001 (Medium Wind)', 'Order of Filter', 'MPSE (dB)', {'Complex', 'Augmented Complex'},[0 30 -12 -11.4], 1);
% set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
% save_figure('../Images/Section_3/Section_3_1_B5_Wind');
% 
% % High Wind (Filter Orders)
% figure(6)
% plot(order,pow2db(mpse_clms(3,:,1)),'-*','LineWidth',2, 'Color', 'red');
% hold on;
% plot(order,pow2db(mpse_aclms(3,:,1)),'-*','LineWidth',2, 'Color', 'green');
% hold off;
% plot_legend_func('\mu=0.0001 (High Wind)', 'Order of Filter', 'MPSE (dB)', {'Complex', 'Augmented Complex'},[0 30 -5.75 -4.75], 1);
% set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
% save_figure('../Images/Section_3/Section_3_1_B6_Wind');

%-----------------------------------------------------------------------
% Question C
%-----------------------------------------------------------------------
% rng(0);
% load Data/colors.mat
% 
% freq_0=50;
% freq_samp=10000;
% N=1000;
% 
% distort_phi_phase=0;
% distortion_b=0;
% distortion_c=0;
% 
% % phases amplitudes 
% V_a=ones(1,N);
% V_b=ones(1,N);
% V_c=ones(1,N);
% v_balanced=clarkson_v_transform(N,freq_0,freq_samp,distort_phi_phase,V_a,V_b,V_c,distortion_b,distortion_c);
% 
% % unbalance because of phase
% distortion_b_unbalanced = [pi/2 pi/3 pi/5 pi/7 pi/4];
% distortion_c_unbalanced = [pi/3 pi/4 pi/6 pi/8 pi/4];
% distort_v_unbalanced = complex(zeros(N,length(distortion_b)));
% 
% for i=1:length(distortion_c_unbalanced)
%     distort_v_unbalanced(:,i)=clarkson_v_transform(N,freq_0,freq_samp,distort_phi_phase,V_a,V_b,V_c,distortion_b_unbalanced(i),distortion_c_unbalanced(i));
% end
% 
% voltage_multipliers=[0.8 1 1.2;0.6 1 1.4;0.4 1 1.6;0.2 1 1.8];
% v_unbalanced_magnitude = complex(zeros(N,length(voltage_multipliers)));
% 
% for i=1:length(voltage_multipliers)
%     v_unbalanced_magnitude(:,i)=clarkson_v_transform(N,freq_0,freq_samp,distort_phi_phase,V_a*voltage_multipliers(i,1),V_b*voltage_multipliers(i,2),V_c*voltage_multipliers(i,3),distortion_b,distortion_c);
% end
% 
% % balanced signal
% figure(1)
% plot(real(v_balanced), imag(v_balanced), 'o','LineWidth',2,'color','red');
% plot_func('Complex Signal: Balanced', 'Re', 'Imag', [-2.5 2.5 -2.5 2.5], 1);
% set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
% save_figure('../Images/Section_3/Section_3_1_C1_Clarkson_V_Transform');
% 
% % Unbalanced Signals Distortion
% figure(2)
% plot(real(distort_v_unbalanced(:,1)), imag(distort_v_unbalanced(:,1)), 'o','LineWidth',2,'color','red');
% hold on;
% plot(real(distort_v_unbalanced(:,2)), imag(distort_v_unbalanced(:,2)), 'o','LineWidth',2,'color','blue');
% plot(real(distort_v_unbalanced(:,3)), imag(distort_v_unbalanced(:,3)), 'o','LineWidth',2,'color','green');
% plot(real(distort_v_unbalanced(:,4)), imag(distort_v_unbalanced(:,4)), 'o','LineWidth',2,'color','cyan');
% hold off;
% plot_func('Complex Signal: Unbalanced phase', 'Re', 'Imag', [-2.5 2.5 -2.5 2.5], 1);
% h = legend({'$\Delta_{b}=\frac{\pi}{2}, \ \Delta_{c}=\frac{\pi}{3}$', '$\Delta_{b}=\frac{\pi}{3}, \ \Delta_{c}=\frac{\pi}{4}$', '$\Delta_{b}=\frac{\pi}{5}, \ \Delta_{c}=\frac{\pi}{6}$', '$\Delta_{b}=\frac{\pi}{7}, \ \Delta_{c}=\frac{\pi}{8}$'}, 'FontSize', 10);
% set(h,'Interpreter','latex');
% set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
% save_figure('../Images/Section_3/Section_3_1_C2_Clarkson_V_Transform');
% 
% % Unbalanced Signals Magnitude
% figure(3)
% plot(real(v_unbalanced_magnitude(:,1)), imag(distort_v_unbalanced(:,1)), 'o','LineWidth',2,'color','red');
% hold on;
% plot(real(v_unbalanced_magnitude(:,2)), imag(distort_v_unbalanced(:,2)), 'o','LineWidth',2,'color','blue');
% plot(real(v_unbalanced_magnitude(:,3)), imag(distort_v_unbalanced(:,3)), 'o','LineWidth',2,'color','green');
% plot(real(v_unbalanced_magnitude(:,4)), imag(distort_v_unbalanced(:,4)), 'o','LineWidth',2,'color','cyan');
% hold off;
% plot_func('Complex Signal: Unbalanced Magnitude', 'Re', 'Imag', [-2.5 2.5 -2.5 2.5], 1);
% h = legend({'$V_a = 0.8, \ V_b = 1, \ V_c = 1.2$', '$V_a = 0.6, \ V_b = 1, \ V_c = 1.4$', '$V_a = 0.4, \ V_b = 1, \ V_c = 1.4$', '$V_a = 0.2, \ V_b = 1, \ V_c = 1.8$' }, 'FontSize', 10);
% set(h,'Interpreter','latex');
% set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
% save_figure('../Images/Section_3/Section_3_1_C3_Clarkson_V_Transform');
% 
% 
% figure(4)
% plot(real(distort_v_unbalanced(:,5)), imag(distort_v_unbalanced(:,5)), '*','LineWidth',2,'color','red');
% hold on;
% plot(real(v_unbalanced_magnitude(:,4)), imag(distort_v_unbalanced(:,4)), '*','LineWidth',2,'color','blue');
% hold off;
% plot_func('Complex Signal: Unbalanced Magnitudes', 'Re', 'Imag', [-2.5 2.5 -2.5 2.5], 1);
% h = legend({'$\Delta_{b}=\frac{\pi}{4}, \ \Delta_{c}=\frac{\pi}{4}$','$V_a = 0.2, \ V_b = 1, \ V_c = 1.8$' }, 'FontSize', 10);
% set(h,'Interpreter','latex');
% set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
% save_figure('../Images/Section_3/Section_3_1_C4_Clarkson_V_Transform');
%-----------------------------------------------------------------------
% Question D
%-----------------------------------------------------------------------
% rng(0);
% mu=0.01;
% order=1;
% 
% % parameters of sinewave
% freq_0=50;
% freq_samp=1000;
% N=10000;
% distort_phi_phase=0;
% distortion_b=0;
% distortion_c=0;
% 
% % phase amplitudes
% V_a=ones(1,N);
% V_b=ones(1,N);
% V_c=ones(1,N);
% 
% % balanced signal
% v_balanced=clarkson_v_transform(N,freq_0,freq_samp,distort_phi_phase,V_a,V_b,V_c,distortion_b,distortion_c);
% 
% % distoring signal
% distortion_b_unbalanced = pi/4;
% distortion_c_unbalanced = pi/4;
% distort_v_unbalanced=clarkson_v_transform(N,freq_0,freq_samp,distort_phi_phase,V_a,V_b,V_c,distortion_b_unbalanced,distortion_c_unbalanced);
% voltage_multipliers=[0.2 1 1.8];
% v_unbalanced_magnitude=clarkson_v_transform(N,freq_0,freq_samp,distort_phi_phase,V_a*voltage_multipliers(1),V_b*voltage_multipliers(2),V_c*voltage_multipliers(3),distortion_b,distortion_c);
% 
% 
% [~,~,h_balanced_clms]=ar_clms(v_balanced,mu,order);
% freq_0_balanced_clms=(freq_samp/(2*pi))*atan(imag(conj(h_balanced_clms))./real(conj(h_balanced_clms)));
% [~,~,h_balanced_aclms,g_balanced_aclms]=ar_aclms(v_balanced,mu,order);
% freq_0_balanced_aclms=(freq_samp/(2*pi))*atan((sqrt( (imag(h_balanced_aclms).^2) - (abs(g_balanced_aclms).^2) ))./ (real(h_balanced_aclms)));
% 
% [~,~,h_unbalanced_clms]=ar_clms(distort_v_unbalanced,mu,order);
% freq_0_unbalanced_clms=(freq_samp/(2*pi))*atan(imag(conj(h_unbalanced_clms))./real(conj(h_unbalanced_clms)));
% [~,~,h_unbalanced_aclms,g_unbalanced_aclms]=ar_aclms(distort_v_unbalanced,mu,order);
% freq_0_unbalanced_aclms=(freq_samp/(2*pi))*atan((sqrt( (imag(h_unbalanced_aclms).^2) - (abs(g_unbalanced_aclms).^2) ))./ (real(h_unbalanced_aclms)));
% 
% [~,~,h_unbalanced_magnitude_clms]=ar_clms(v_unbalanced_magnitude,mu,order);
% freq_0_unbalanced_magnitude_clms=(freq_samp/(2*pi))*atan(imag(conj(h_unbalanced_magnitude_clms))./real(conj(h_unbalanced_magnitude_clms)));
% [~,~,h_unbalanced_magnitude_aclms,g_unbalanced_magnitude_aclms]=ar_aclms(distort_v_unbalanced,mu,order);
% freq_0_unbalanced_magnitude_aclms=(freq_samp/(2*pi))*atan((sqrt( (imag(h_unbalanced_magnitude_aclms).^2) - (abs(g_unbalanced_magnitude_aclms).^2) ))./ (real(h_unbalanced_magnitude_aclms)));
% 
% % Balanced System freq
% figure(1)
% start=3;
% plot(freq_0_balanced_clms(start:end),'LineWidth',2,'color','red');
% hold on;
% plot(abs(freq_0_balanced_aclms(start:end)),'LineWidth',2,'color','green');
% hold off;
% plot_legend_func('Complex Signal: Balanced','n','Frequency, f_{o}',{'Complex', 'Augmented Complex'},[1 420 0 300], 1);
% set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
% save_figure('../Images/Section_3/Section_3_1_D1_Convergence_aclms_clms');
% 
% % Unbalanced System Freq, Phase
% figure(2)
% plot(freq_0_unbalanced_clms(start:end),'LineWidth',2,'color','red');
% hold on;
% plot(abs(freq_0_unbalanced_aclms(start:end)),'LineWidth',2,'color','green');
% hold off;
% plot_legend_func('Complex Signal: Phase Unbalanced (\Delta_{b}=\pi/4, \Delta_{c}=\pi/4)','n','Frequency, f_{o}',{'Complex', 'Augmented Complex'},[1 600 -10 250], 1);
% set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
% save_figure('../Images/Section_3/Section_3_1_D2_Convergence_aclms_clms');
% 
% % Unbalanced System Freq, Magnitude
% figure(3)
% plot(freq_0_unbalanced_magnitude_clms(start:end),'LineWidth',2,'color','red');
% hold on;
% plot(abs(freq_0_unbalanced_magnitude_aclms(start:end)),'LineWidth',2,'color','green');
% hold off;
% plot_legend_func('Complex Signal: Mag Unbalanced (V_a=0.2, V_b=1, V_c=1.8)','n','Frequency, f_{o}',{'Complex', 'Augmented Complex'},[1 600 -10 250], 1);
% set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
% save_figure('../Images/Section_3/Section_3_1_D3_Convergence_aclms_clms');

