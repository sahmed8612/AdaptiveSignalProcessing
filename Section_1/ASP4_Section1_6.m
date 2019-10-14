%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% Section 1.6
%-----------------------------------------------------------------------
%
%-----------------------------------------------------------------------
% Question A  
%-----------------------------------------------------------------------
% pcr = load('Data/PCR/PCAPCR.mat');
% 
% % Clean matrix (SVD)
% S = svd(pcr.X);
% % Matrix with Noise (SVD)
% Sig_noise = svd(pcr.Xnoise);
% error = (S - Sig_noise).^2;
% 
% bar([S, Sig_noise]);
% legend(["Original Signal", "Noisy Signal"]);
% title(sprintf("Singular Values: X & {X}_{noise}"),'FontSize',20);
% set(gca,'fontsize',20);
% ylabel("Singular Values Magnitude",'FontSize',20);
% xlabel("Singular Value",'FontSize',20);
% set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
% save_figure('../Images/Section_1/Section_1_6_A1_Robust_Regression_Singular_Values');
% 
% bar(error, 'FaceColor', 'red');
% title(sprintf("Singular Values: Squared Prediction Error"),'FontSize',20);
% set(gca,'fontsize',20);
% xlabel("Singular Value",'FontSize',20);
% ylabel("Squared Error",'FontSize',20);
% set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
% save_figure('../Images/Section_1/Section_1_6_A2_Robust_Regression_Singular_Values');

%-----------------------------------------------------------------------
% Question B  
%-----------------------------------------------------------------------
% pcr = load('Data/PCR/PCAPCR.mat');
% 
% [m, n] = size(pcr.Xnoise);
% [Unoise, Snoise, Vnoise] = svd(pcr.Xnoise);
% 
% % low-rank approximation
% r = 1:10;
% forb_norm = zeros(length(r), 1);
% power = zeros(length(r), 1);
% 
% for i = 1:length(r)
%     X_noise_hat = Unoise(1:m, 1:r(i)) * Snoise(1:r(i), 1:r(i)) * Vnoise(1:n, 1:r(i))';
%     forb_norm(i) = norm(X_noise_hat-pcr.X, 'fro');
%     sigma = diag(Snoise);
%     power(i) = sum(sum(sigma(r(i)+1:end).^2));
% 
% end
% 
% plot(r, forb_norm,'LineWidth', 2,'color','red');
% set(gca,'fontsize',20);
% title("Approximation Error: Low Rank");
% xlabel("Low Rank Approximation");
% ylabel("Appromixmation Error");
% xticks(r);
% set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
% save_figure('../Images/Section_1/Section_1_6_B1_Robust_Regression_Singular_Values');

%-----------------------------------------------------------------------
% Question C  
%-----------------------------------------------------------------------
% pcr = load('Data/PCR/PCAPCR.mat');
% 
% %OLS Weights and Errors
% OLS_B = pcr.Xnoise' * pcr.Xnoise \ pcr.Xnoise' * pcr.Y;
% OLS_Y = pcr.Xnoise * OLS_B;
% OLS_Y_test = pcr.Xtest * OLS_B;
% OLS_training_error = norm(pcr.Y - OLS_Y, 'fro');
% OLS_testing_error = norm(pcr.Ytest - OLS_Y_test, 'fro');
% 
% 
% %PCR Weights and Errors
% [m, n] = size(pcr.Xnoise);
% r = 1:10;
% PCR_training_error = zeros(length(r), 1);
% PCR_test_error = zeros(length(r), 1);
% 
% for i = 1:length(r)
%     [Unoise, Snoise, Vnoise] = svd(pcr.Xnoise);
%     [Utest, Stest, Vtest] = svd(pcr.Xtest);
%     
%     Xnoise_tilda = Unoise(1:m, 1:r(i)) * Snoise(1:r(i), 1:r(i)) * Vnoise(1:n, 1:r(i))';
%     Xtest_tilda = Utest(1:m, 1:r(i)) * Stest(1:r(i), 1:r(i)) * Vtest(1:n, 1:r(i))';
% 
%     PCR_B = Vnoise(1:n, 1:r(i)) / Snoise(1:r(i), 1:r(i)) * Unoise(1:m, 1:r(i))' * pcr.Y;
%     Y_PCR = Xnoise_tilda * PCR_B;
%     Ytest_PCR = Xtest_tilda * PCR_B;
%     % PCR Errors
%     PCR_training_error(i) = norm(pcr.Y - Y_PCR, 'fro');
%     PCR_test_error(i) = norm(pcr.Ytest - Ytest_PCR, 'fro');
%     
% end
% 
% figure;
% plot(r, PCR_training_error, "DisplayName", "PCR",'LineWidth', 2,'color','red');
% hold on;
% plot([r(1) r(end)], [OLS_training_error OLS_training_error], "LineStyle", "--", "DisplayName", "OLS",'LineWidth', 2)
% title("Training Error: OLS & PCR",'FontSize',20);
% set(gca,'fontsize',20);
% xlabel("Low Rank Approximation",'FontSize',20);
% ylabel("Training Error",'FontSize',20);
% xticks(r);
% legend("show");
% set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
% save_figure('../Images/Section_1/Section_1_6_C1_Robust_Regression_Singular_Values');
% 
% figure;
% plot(r, PCR_test_error, "DisplayName", "PCR",'LineWidth', 2,'color','red');
% hold on;
% plot([r(1) r(end)], [OLS_testing_error OLS_testing_error], "LineStyle", "--", "DisplayName", "OLS",'LineWidth', 2)
% title("Testing Error: OLS & PCR",'FontSize',20);
% set(gca,'fontsize',20);
% xlabel("Low Rank Approximation",'FontSize',20);
% ylabel("Test Error",'FontSize',20);
% xticks(r);
% legend("show");
% set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
% save_figure('../Images/Section_1/Section_1_6_C2_Robust_Regression_Singular_Values');

%-----------------------------------------------------------------------
% Question D 
%-----------------------------------------------------------------------
% pcr = load('Data/PCR/PCAPCR.mat');
% iter = 100;
% 
% % OLS Weights and errors
% OLS_B = pcr.Xnoise' * pcr.Xnoise \ pcr.Xnoise' * pcr.Y;
% OLS_error = zeros(iter, 1);
% 
% for j = 1:iter
%     [Y, Y_OLS] = regval(OLS_B);
%     OLS_error(j) = norm(Y - Y_OLS, 2);
% end
% 
% % PCR Weights and errors
% [m, n] = size(pcr.Xnoise);
% r = 1:10;
% PCR_error = zeros(iter, length(r));
% 
% for i = 1:length(r)
%     [Unoise, Snoise, Vnoise] = svd(pcr.Xnoise);
%     Xnoise_tilda = Unoise(1:m, 1:r(i)) * Snoise(1:r(i), 1:r(i)) * Vnoise(1:n, 1:r(i))';
%     PCR_B = Vnoise(1:n, 1:r(i)) / Snoise(1:r(i), 1:r(i)) * Unoise(1:m, 1:r(i))' * pcr.Y;
%     for j = 1:iter
%         [Y, Y_PCR] = regval(PCR_B);
%         PCR_error(j, i) = norm(Y - Y_PCR, 2);
%     end
%     
% end
% 
% for j = 1:iter
%     plot([r(1) r(end)], [OLS_error(j) OLS_error(j)], 'Color', 'red', 'LineWidth', 0.5);
%     hold on;
% end
% for j = 1:iter
%     plot(r, PCR_error(j, 1:end), 'Color', 'blue', 'LineWidth', 0.5);
%     hold on;
% end
% p1 = plot(r, mean(PCR_error, 1), 'Color','red', "DisplayName", "PCR");
% hold on;
% p2 = plot([r(1) r(end)], [mean(OLS_error) mean(OLS_error)], 'Color', 'blue', "LineStyle", "--", "DisplayName", "OLS");
% title("Model Error Analysis: OLS & PCR");
% set(gca,'fontsize',20);
% xlabel("Low Rank Approximation");
% ylabel("Model Prediction Errors");
% xticks(r);
% legend([p1 p2]);
% set(gca,'Box','off','TickDir','out','TickLength',[.01 .01],'XMinorTick','on','YMinorTick','on','XGrid','on','YGrid','on','XColor',[.3 .3 .3],'YColor', [.3 .3 .3],'LineWidth',1);
% save_figure('../Images/Section_1/Section_1_6_D1_Robust_Regression_Singular_Values');

