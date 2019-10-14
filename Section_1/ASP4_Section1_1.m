%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% Section 1.1
%-----------------------------------------------------------------------
%
%-----------------------------------------------------------------------
% Question A  
%-----------------------------------------------------------------------
% No Zero Padding when size_M=10, length=20 (DFT)
% size_M=10;
% length=20;
% x=x_acf(size_M,length);
% p=abs(fftshift(fft(x)));
% diff_f=2*pi/length;
% f=-pi:diff_f:pi-diff_f;
% 
% figure(1);
% s=plot(f,p,'LineWidth',2,'color','red');
% set(s,'Marker','none');
% set(gca,'fontsize',20);
% title('No Zero-Padding: M=10, length=20,','FontSize',20);
% xlabel('Frequency (rad/sample)', 'FontSize', 20);
% ylabel('Magnitude', 'FontSize', 20);
% axis([-pi pi 0 20]);
% set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
% save_figure('../Images/Section_1/Section_1_1_A1_Zero_Padding_Effects');
% 
% % Zero Padding when size_M=10, length=256 (DFT)
% 
% size_M=10;
% length=256;
% x=x_acf(size_M,length);
% p=abs(fftshift(fft(x)));
% diff_f=2*pi/length;
% f=-pi:diff_f:pi-diff_f;
% 
% figure(2);
% s=plot(f,p,'LineWidth',2,'color','red');
% set(s, 'Marker', 'none');
% set(gca,'fontsize',20);
% title('Zero-Padding: M=10, length=256','FontSize',20);
% xlabel('Frequency (rad/sample)', 'FontSize', 20);
% ylabel('Magnitude', 'FontSize', 20);
% axis([-pi pi 0 20]);
% set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
% save_figure('../Images/Section_1/Section_1_1_A2_Zero_Padding_Effects');
% 
% % No Zero Padding when size_M=128, length=256 (DFT)
% 
% size_M=128;
% length=256;
% x=x_acf(size_M,length);
% p=abs(fftshift(fft(x)));
% diff_f=2*pi/length;
% f=-pi:diff_f:pi-diff_f;
% 
% figure(3);
% s=plot(f,p,'LineWidth',2,'color','red');
% set(s, 'Marker', 'none');
% set(gca,'fontsize',20);
% title('No Zero-Padding: M=128, length=256','FontSize',20);
% xlabel('Frequency (rad/sample)', 'FontSize', 20);
% ylabel('Magnitude', 'FontSize', 20);
% axis([-pi pi 0 200]);
% set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
% save_figure('../Images/Section_1/Section_1_1_A3_Zero_Padding_Effects');
% 
% % No Zero Padding when size_M=128, length=256 (DFT)
% 
% size_M=128;
% length=2048;
% x=x_acf(size_M,length);
% p=abs(fftshift(fft(x)));
% diff_f=2*pi/length;
% f=-pi:diff_f:pi-diff_f;
% 
% figure(4);
% s=plot(f,p,'LineWidth',2,'color','red');
% set(s, 'Marker', 'none');
% set(gca,'fontsize',20);
% title('Zero-Padding: M=128, length=2048','FontSize',20);
% xlabel('Frequency (rad/sample)', 'FontSize', 20);
% ylabel('Magnitude', 'FontSize', 20);
% axis([-pi pi 0 200]);
% set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
% save_figure('../Images/Section_1/Section_1_1_A4_Zero_Padding_Effects');

%-----------------------------------------------------------------------
% Question B
%-----------------------------------------------------------------------
M=10;
L=128;
x=x_acf(M,L);
p=fftshift(fft(x));
df=2*pi/L;
f=-pi:df:pi-df;

% FFT Real part
figure(1);
s=plot(f,real(p),'LineWidth',2,'color','red');
set(s, 'Marker', 'none');
ylim([-10 15]);
xlim([-pi pi]);
set(gca,'fontsize',20);
title('DFT Real Part','FontSize',20);
xlabel('Frequency (rad/s)', 'FontSize', 20);
ylabel('Magnitude', 'FontSize', 20);
set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
save_figure('../Images/Section_1/Section_1_1_B1_Imag_Component_Analysis');

% FFT Imag part
figure(2);
s=plot(f,imag(p),'LineWidth',2,'color','red');
xlim([-pi pi]);
set(s, 'Marker', 'none');
set(gca,'fontsize',20);
title('DFT Imaginary Part','FontSize',20);
xlabel('Frequency (rad/s)', 'FontSize', 20);
ylabel('Magnitude', 'FontSize', 20);
set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
save_figure('../Images/Section_1/Section_1_1_B2_Imag_Component_Analysis');

%-----------------------------------------------------------------------
% Question C
%-----------------------------------------------------------------------
% ACF Length
M=10;
L=128;
z=z_acf(M,L);
p=fftshift(fft(z));
df=2*pi/L;
f=-pi:df:pi-df;

% FFT real part
figure(3);
s=plot(f,real(p),'LineWidth',2,'color','red');
set(s, 'Marker', 'none');
ylim([-10 15])
xlim([-pi pi]);
set(gca,'fontsize',20);
title('Real Component of DFT','FontSize',20);
xlabel('Frequency (rad/s)', 'FontSize', 20);
ylabel('Magnitude', 'FontSize', 20);
set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
save_figure('../Images/Section_1/Section_1_1_C1_Imag_Component_Analysis');

% FFT imaginary part
figure(4);
s=plot(f,imag(p),'LineWidth',2,'color','red');
set(s, 'Marker', 'none');
ylim([-10 15])
xlim([-pi pi]);
set(gca,'fontsize',20);
title('DFT Imaginary Component','FontSize',20);
xlabel('Frequency (rad/s)', 'FontSize', 20);
ylabel('Magnitude', 'FontSize', 20);
set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
save_figure('../Images/Section_1/Section_1_1_C2_Imag_Component_Analysis');

% Asymmetrical ACF (Without Imag part)
figure(5);
s=plot(f,abs(real(p)),'LineWidth',2,'color','red');
set(s, 'Marker', 'none');
ylim([0 15])
xlim([-pi pi]);
set(gca,'fontsize',20);
title('No Imaginary Component - PSD','FontSize',20);
xlabel('Frequency (rad/s)', 'FontSize', 20);
ylabel('Magnitude', 'FontSize', 20);
set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
save_figure('../Images/Section_1/Section_1_1_C3_Imag_Component_Analysis');

% Asymmetrical ACF (With Imag part)
figure(6);
s=plot(f,abs(p),'LineWidth',2,'color','red');
set(s, 'Marker', 'none');
ylim([0 15])
xlim([-pi pi]);
set(gca,'fontsize',20);
title('Inc. Imaginary Component - PSD','FontSize',20);
xlabel('Frequency (rad/s)', 'FontSize', 20);
ylabel('Magnitude', 'FontSize', 20);
set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
save_figure('../Images/Section_1/Section_1_1_C4_Imag_Component_Analysis');

