clc; clear; close all;
% -------------------------------------------------------------------------
% Purpose: Replicating figures and tables in Section 7 and Appendices D,E
%          of the paper
%          "Understanding regressions with observations collected at
%           high-frequency over long span" by Chang, Lu and Park (2024)
% -------------------------------------------------------------------------
% Outputs:
%     Figures saved in '../Figures':
%       Figure 6: Simulated Means of Robust Wald Test Statistics (Section 7)
%       Figure 7: Simulated Rejection Probabilities of Robust Wald Tests (Appendix D)
%       Figure 8: Simulated Means of Robust Wald Test Statistics and Instability
%                 of Test Results Caused by One-Day Change in Sampling Frequency (Appendix D)
%       Figure 9: Simulated Means of Robust Wald Test Statistics Under Joint Asymptotics
%                 (Appendix D)
%       Figure 10: Simulated Means of Robust Wald Test Statistics
%                 (Additional simulation in Appendix E)
%       Figure 11: Simulated Rejection Probabilities of Robust Wald Tests
%                 (Additional simulation in Appendix E)
%     Table saved in '../Tables':
%       Tables13: Simulated Size-Adjusted Power (Appendix D)
% -------------------------------------------------------------------------
% This version: September 2024 by Ye Lu (ye.lu1@sydney.edu.au)
% -------------------------------------------------------------------------
set(0,'DefaultTextInterpreter','latex'); set(0,'DefaultAxesTickLabelInterpreter','latex'); set(0,'DefaultLegendInterpreter', 'latex')
set(0,'defaultAxesXGrid','on'); set(0,'defaultAxesYGrid','on'); set(0,'defaultAxesFontSize',13); set(0,'DefaultLineLineWidth',1.5)
x0=0; y0=0; width=800; height=500; savefigfolder = '../Figures';
% -------------------------------------------------------------------------

tic
folderpath = '../MC/'; % Folder to access MC data

%% Section 7: (Main) Simulation
prefix = 'Section_7_';

%% Figure 6:
filename = 'Htest_stat_T30.csv'; Htest_stat_T30 = readmatrix(strcat(folderpath,prefix,filename));
filename = 'Gtest_stat_T30.csv'; Gtest_stat_T30 = readmatrix(strcat(folderpath,prefix,filename));
ad = 1; nw = 2; rt = 3; crt = 4; delta = (1:63)/252;
figure; set(gcf, 'position', [x0, y0, width, height]); xl=0.25; ylab = 'mean of test statistics';
subplot(2,2,1); box on; hold on; plot(delta, Htest_stat_T30(nw,:), '-.'); plot(delta, Htest_stat_T30(ad,:), '-'); hold off
legend('$H$-test: NW', '$H$-test: AD'); xlim([0 xl]); ylim([2 7]); xlabel('$\delta$'); ylabel(ylab)
subplot(2,2,2); box on; hold on; plot(delta, Htest_stat_T30(rt,:), '-.'); plot(delta, Htest_stat_T30(crt,:), '-'); hold off
legend('$H$-test: RT', '$H$-test: CRT'); xlim([0 xl]); xlabel('$\delta$'); ylabel(ylab)
subplot(2,2,3); box on; hold on; plot(delta, Gtest_stat_T30(nw,:), '-.'); plot(delta, Gtest_stat_T30(ad,:), '-'); hold off
legend('$G$-test: NW', '$G$-test: AD'); xlim([0 xl]); ylim([2 20]); xlabel('$\delta$'); ylabel(ylab)
subplot(2,2,4); box on; hold on; plot(delta, Gtest_stat_T30(rt,:), '-.'); plot(delta, Gtest_stat_T30(crt,:), '-'); hold off
legend('$G$-test: RT', '$G$-test: CRT'); xlim([0 xl]); ylim([0 160]); xlabel('$\delta$'); ylabel(ylab)
saveas(gcf, sprintf('%s/Figure6.png', savefigfolder)); % saveas(gcf, sprintf('%s/Figure6', savefigfolder), 'epsc')

%% Figure 7:
filename = 'Htest_rej_T30.csv'; Htest_rej_T30 = readmatrix(strcat(folderpath,prefix,filename));
filename = 'Htest_rej_T50.csv'; Htest_rej_T50 = readmatrix(strcat(folderpath,prefix,filename));
filename = 'Gtest_rej_T30.csv'; Gtest_rej_T30 = readmatrix(strcat(folderpath,prefix,filename));
filename = 'Gtest_rej_T50.csv'; Gtest_rej_T50 = readmatrix(strcat(folderpath,prefix,filename));
figure; set(gcf, 'position', [x0, y0, width, height]); ylab = 'rejection probabilities';
subplot(2,2,1); box on; hold on
plot(delta, Htest_rej_T30(nw,:), '-.'); plot(delta, Htest_rej_T30(ad,:), '-'); plot(delta, Htest_rej_T30(rt,:), ':', 'LineWidth', 2.3); plot(delta, Htest_rej_T30(crt,:), '--'); hold off
legend('$H$-test: NW', '$H$-test: AD', '$H$-test: RT', '$H$-test: CRT'); xlim([0 xl]); ylim([0 1]); xlabel('$\delta$'); ylabel(ylab); title('$T=30$')
subplot(2,2,2); box on; hold on
plot(delta, Htest_rej_T50(nw,:), '-.'); plot(delta, Htest_rej_T50(ad,:), '-'); plot(delta, Htest_rej_T50(rt,:), ':', 'LineWidth', 2.3); plot(delta, Htest_rej_T50(crt,:), '--'); hold off
legend('$H$-test: NW', '$H$-test: AD', '$H$-test: RT', '$H$-test: CRT'); xlim([0 xl]); ylim([0 1]); xlabel('$\delta$'); ylabel(ylab); title('$T=50$')
subplot(2,2,3); box on; hold on
plot(delta, Gtest_rej_T30(nw,:), '-.'); plot(delta, Gtest_rej_T30(ad,:), '-'); plot(delta, Gtest_rej_T30(rt,:), ':', 'LineWidth', 2.3); plot(delta, Gtest_rej_T30(crt,:), '--'); hold off
legend('$G$-test: NW', '$G$-test: AD', '$G$-test: RT', '$G$-test: CRT'); xlim([0 xl]); ylim([0 1]); xlabel('$\delta$'); ylabel(ylab); title('$T=30$')
subplot(2,2,4); box on; hold on
plot(delta, Gtest_rej_T50(nw,:), '-.'); plot(delta, Gtest_rej_T50(ad,:), '-'); plot(delta, Gtest_rej_T50(rt,:), ':', 'LineWidth', 2.3); plot(delta, Gtest_rej_T50(crt,:), '--'); hold off
legend('$G$-test: NW', '$G$-test: AD', '$G$-test: RT', '$G$-test: CRT'); xlim([0 xl]); ylim([0 1]); xlabel('$\delta$'); ylabel(ylab); title('$T=50$')
saveas(gcf, sprintf('%s/Figure7.png', savefigfolder)); % saveas(gcf, sprintf('%s/Figure7', savefigfolder), 'epsc')

%% Figure 8:
filename = 'Htest_instability_T30.csv'; Htest_instability = readmatrix(strcat(folderpath,prefix,filename));
filename = 'Gtest_instability_T30.csv'; Gtest_instability = readmatrix(strcat(folderpath,prefix,filename));
co = [0.8500 0.3250 0.0980]; cp = [0 0.4470 0.7410]; % line colors
figure; set(gcf, 'position', [x0, y0, width, height])
subplot(2,2,1); box on; hold on
plot(delta, Htest_stat_T30(ad,:), 'Color', co, 'Linestyle', '-'); plot(delta, Htest_stat_T30(crt,:), 'Color',cp, 'LineStyle','-.', 'LineWidth',1.8); hold off;
legend('$H$-test: AD','$H$-test: CRT','Location','northeast'); ylim([2 4]); xlabel('$\delta$'); title('Simulated mean of test statistics')
subplot(2,2,2); box on; hold on
plot(delta, Gtest_stat_T30(ad,:), 'Color', co, 'Linestyle', '-'); plot(delta, Gtest_stat_T30(crt,:), 'Color',cp, 'LineStyle','-.', 'LineWidth',1.8); hold off;
legend('$H$-test: AD','$H$-test: CRT','Location','northeast'); ylim([3.45 3.6]); xlabel('$\delta$'); title('Simulated mean of test statistics')
subplot(2,2,3);box on; hold on
plot(delta(2:62), Htest_instability(ad,2:62), 'Color',co); plot(delta(2:62), Htest_instability(crt,2:62), 'Color',cp, 'LineStyle','-.', 'LineWidth',1.8); hold off
legend('$H$-test: AD','$H$-test: CRT','Location','northeast'); xlim([0 xl]); ylabel('percent'); xlabel('$\delta$'); ylim([0 30]); title('Instability of test results')
subplot(2,2,4);box on; hold on
plot(delta(2:62), Gtest_instability(ad,2:62), 'Color',co); plot(delta(2:62), Gtest_instability(crt,2:62), 'Color',cp, 'LineStyle','-.', 'LineWidth',1.8); hold off
legend('$G$-test: AD','$G$-test: CRT','Location','northeast'); xlim([0 xl]); ylabel('percent'); xlabel('$\delta$'); ylim([0 10]); title('Instability of test results')
saveas(gcf, sprintf('%s/Figure8.png', savefigfolder)); % saveas(gcf, sprintf('%s/Figure8', savefigfolder), 'epsc')

%% Figure 9:
filename = 'Htest_stat_Double.csv'; Htest_stat = readmatrix(strcat(folderpath,prefix,filename));
filename = 'Gtest_stat_Double.csv'; Gtest_stat = readmatrix(strcat(folderpath,prefix,filename));
figure; set(gcf, 'position', [x0, y0, width, height]); ylab = 'mean of test statistics';
subplot(2,2,1); box on; hold on; plot(delta, Htest_stat(nw,:), '-.'); plot(delta, Htest_stat(ad,:), '-'); hold off
legend('$H$-test: NW', '$H$-test: AD'); xlim([0 xl]); ylim([0 20]); xlabel('$\delta$'); ylabel(ylab)
subplot(2,2,2); box on; hold on; plot(delta, Htest_stat(rt,:), '-.'); plot(delta, Htest_stat(crt,:), '-'); hold off
legend('$H$-test: RT', '$H$-test: CRT'); xlim([0 xl]); xlabel('$\delta$'); ylabel(ylab)
subplot(2,2,3); box on; hold on; plot(delta, Gtest_stat(nw,:), '-.'); plot(delta, Gtest_stat(ad,:), '-'); hold off
legend('$G$-test: NW', '$G$-test: AD'); xlim([0 xl]); ylim([0 20]); xlabel('$\delta$'); ylabel(ylab)
subplot(2,2,4); box on; hold on; plot(delta, Gtest_stat(rt,:), '-.'); plot(delta, Gtest_stat(crt,:), '-'); hold off
legend('$G$-test: RT', '$G$-test: CRT'); xlim([0 xl]); xlabel('$\delta$'); ylabel(ylab)
saveas(gcf, sprintf('%s/Figure9.png', savefigfolder)); % saveas(gcf, sprintf('%s/Figure9', savefigfolder), 'epsc')

%% Table 13:
power = cell(2,2,2); % two tests, two alternatives, and two bandwidths
iH =1; iG = 2; iad = 1; icrt = 2; del = [1 21 63]; df = 1; cv005 = chi2inv(0.95,df); cv001 = chi2inv(0.99,df); % nominal critical values

% -- H-test for stationary regression model --
filename = 'Htest_bet_1_ad.csv';  loadname{iad}  = strcat(folderpath,prefix,filename);
filename = 'Htest_bet_1_crt.csv'; loadname{icrt} = strcat(folderpath,prefix,filename);

% Obtain *adjusted critical values* which yield 5% size:
cv = (cv005:0.001:cv001); n = numel(cv); rejprob = zeros(n,numel(del),2); a = zeros(2,3); b = zeros(size(a)); cv_adj  = zeros(2,3);
for i = 1:2
    t = readmatrix(loadname{i});
    for j = 1:n
        rejprob(j,:,i) = mean(t>cv(j));
    end
    [a(i,:), b(i,:)] = min(abs(rejprob(:,:,i)-0.05));
end
cv_adj(iad,:)  = [cv(b(iad,1)) cv(b(iad,2)) cv(b(iad,3))]; cv_adj(icrt,:) = [cv(b(icrt,1)) cv(b(icrt,2)) cv(b(icrt,3))];

% Calculate size-adjusted power under two alternatives using the adjusted critical values:
filename = 'Htest_bet_095_ad.csv';  loadname{iad}  = strcat(folderpath,prefix,filename);
filename = 'Htest_bet_095_crt.csv'; loadname{icrt} = strcat(folderpath,prefix,filename);
for i = 1:2
    t = readmatrix(loadname{i}); power{iH,1,i} = mean(t > cv_adj(i,:));
end
filename = 'Htest_bet_090_ad.csv';  loadname{iad}  = strcat(folderpath,prefix,filename);
filename = 'Htest_bet_090_crt.csv'; loadname{icrt} = strcat(folderpath,prefix,filename);
for i = 1:2
    t = readmatrix(loadname{i}); power{iH,2,i} = mean(t > cv_adj(i,:));
end

% -- G-test for cointegrating regression model --
filename = 'Gtest_bet_1_ad.csv';  loadname{iad}  = strcat(folderpath,prefix,filename);
filename = 'Gtest_bet_1_crt.csv'; loadname{icrt} = strcat(folderpath,prefix,filename);

% Obtain *adjusted critical values* which yield 5% size:
cv = (cv005:0.001:cv001); n = numel(cv); rejprob = zeros(n,numel(del),2); a = zeros(2,3); b = zeros(size(a)); cv_adj  = zeros(2,3);
for i = 1:2
    t = readmatrix(loadname{i});
    for j = 1:n
        rejprob(j,:,i) = mean(t>cv(j));
    end
    [a(i,:), b(i,:)] = min(abs(rejprob(:,:,i)-0.05));
end
cv_adj(iad,:)  = [cv(b(iad,1)) cv(b(iad,2)) cv(b(iad,3))]; cv_adj(icrt,:) = [cv(b(icrt,1)) cv(b(icrt,2)) cv(b(icrt,3))];

% Calculate size-adjusted power under two alternatives using the adjusted critical values:
filename = 'Gtest_bet_099_ad.csv';  loadname{iad}  = strcat(folderpath,prefix,filename);
filename = 'Gtest_bet_099_crt.csv'; loadname{icrt} = strcat(folderpath,prefix,filename);
for i = 1:2
    t = readmatrix(loadname{i}); power{iG,1,i} = mean(t > cv_adj(i,:));
end
filename = 'Gtest_bet_0985_ad.csv';  loadname{iad}  = strcat(folderpath,prefix,filename);
filename = 'Gtest_bet_0985_crt.csv'; loadname{icrt} = strcat(folderpath,prefix,filename);
for i = 1:2
    t = readmatrix(loadname{i}); power{iG,2,i} = mean(t > cv_adj(i,:));
end

% -- Display Table 13 --
fileID = fopen('../Tables/Table13.txt','w'); TableDisplay(power, fileID); fclose(fileID);

%% Appendix E: Additional Simulation Results
prefix = 'Appendix_E_';

%% Figure 10:
filename = 'Htest_stat_T30.csv'; Htest_stat_T30 = readmatrix(strcat(folderpath,prefix,filename));
filename = 'Htest_stat_T50.csv'; Htest_stat_T50 = readmatrix(strcat(folderpath,prefix,filename));
filename = 'Gtest_stat_T30.csv'; Gtest_stat_T30 = readmatrix(strcat(folderpath,prefix,filename));
filename = 'Gtest_stat_T50.csv'; Gtest_stat_T50 = readmatrix(strcat(folderpath,prefix,filename));
figure; set(gcf, 'position', [x0, y0, width, height])
subplot(2,2,1); box on; hold on; plot(delta, Htest_stat_T30(nw,:), '-.'); plot(delta, Htest_stat_T30(ad,:), '-'); hold off
legend('$H$-test: NW', '$H$-test: AD'); xlim([0 xl]); ylim([2 6]); xlabel('$\delta$'); ylabel(ylab)
subplot(2,2,2); box on; hold on; plot(delta, Htest_stat_T30(rt,:), '-.'); plot(delta, Htest_stat_T30(crt,:), '-'); hold off
legend('$H$-test: RT', '$H$-test: CRT'); xlim([0 xl]); xlabel('$\delta$'); ylabel(ylab)
subplot(2,2,3); box on; hold on; plot(delta, Gtest_stat_T30(nw,:), '-.'); plot(delta, Gtest_stat_T30(ad,:), '-'); hold off
legend('$G$-test: NW', '$G$-test: AD'); xlim([0 xl]); ylim([2 3]); xlabel('$\delta$'); ylabel(ylab)
subplot(2,2,4); box on; hold on; plot(delta, Gtest_stat_T30(rt,:), '-.'); plot(delta, Gtest_stat_T30(crt,:), '-'); hold off
legend('$G$-test: RT', '$G$-test: CRT'); xlim([0 xl]); ylim([0 14]); xlabel('$\delta$'); ylabel(ylab)
saveas(gcf, sprintf('%s/Figure10.png', savefigfolder)); % saveas(gcf, sprintf('%s/Figure10', savefigfolder), 'epsc')

%% Figure 11:
filename = 'Htest_rej_T30.csv'; Htest_rej_T30 = readmatrix(strcat(folderpath,prefix,filename));
filename = 'Htest_rej_T50.csv'; Htest_rej_T50 = readmatrix(strcat(folderpath,prefix,filename));
filename = 'Gtest_rej_T30.csv'; Gtest_rej_T30 = readmatrix(strcat(folderpath,prefix,filename));
filename = 'Gtest_rej_T50.csv'; Gtest_rej_T50 = readmatrix(strcat(folderpath,prefix,filename));
figure; set(gcf, 'position', [x0, y0, width, height]); ylab = 'rejection probabilities';
subplot(2,2,1); box on; hold on
plot(delta, Htest_rej_T30(nw,:), '-.'); plot(delta, Htest_rej_T30(ad,:), '-'); plot(delta, Htest_rej_T30(rt,:), ':', 'LineWidth', 2.3); plot(delta, Htest_rej_T30(crt,:), '--'); hold off
legend('$H$-test: NW', '$H$-test: AD', '$H$-test: RT', '$H$-test: CRT'); xlim([0 xl]); ylim([0 1]); xlabel('$\delta$'); ylabel(ylab); title('$T=30$')
subplot(2,2,2); box on; hold on
plot(delta, Htest_rej_T50(nw,:), '-.'); plot(delta, Htest_rej_T50(ad,:), '-'); plot(delta, Htest_rej_T50(rt,:), ':', 'LineWidth', 2.3); plot(delta, Htest_rej_T50(crt,:), '--'); hold off
legend('$H$-test: NW', '$H$-test: AD', '$H$-test: RT', '$H$-test: CRT'); xlim([0 xl]); ylim([0 1]); xlabel('$\delta$'); ylabel(ylab); title('$T=50$')
subplot(2,2,3); box on; hold on
plot(delta, Gtest_rej_T30(nw,:), '-.'); plot(delta, Gtest_rej_T30(ad,:), '-'); plot(delta, Gtest_rej_T30(rt,:), ':', 'LineWidth', 2.3); plot(delta, Gtest_rej_T30(crt,:), '--'); hold off
legend('$G$-test: NW', '$G$-test: AD', '$G$-test: RT', '$G$-test: CRT'); xlim([0 xl]); ylim([0 1]); xlabel('$\delta$'); ylabel(ylab); title('$T=30$')
subplot(2,2,4); box on; hold on
plot(delta, Gtest_rej_T50(nw,:), '-.'); plot(delta, Gtest_rej_T50(ad,:), '-'); plot(delta, Gtest_rej_T50(rt,:), ':', 'LineWidth', 2.3); plot(delta, Gtest_rej_T50(crt,:), '--'); hold off
legend('$G$-test: NW', '$G$-test: AD', '$G$-test: RT', '$G$-test: CRT'); xlim([0 xl]); ylim([0 1]); xlabel('$\delta$'); ylabel(ylab); title('$T=50$')
saveas(gcf, sprintf('%s/Figure11.png', savefigfolder)); % saveas(gcf, sprintf('%s/Figure11', savefigfolder), 'epsc')

toc

%% Auxiliary functions for displaying Table 13:

% Auxiliary function 1
function TableDisplay(power, fileID)
iH = 1; iG = 2; l = 65; % length of the table
disp(repmat('-',1,l));
fprintf(fileID, strcat(repmat('-',1,l),'\n'));
fprintf(strcat(repmat('%12s',1,5),'\n'), 'H-test', '', 'Daily','Monthly','Quarterly');
fprintf(fileID, strcat(repmat('%12s',1,5),'\n'), 'H-test', '', 'Daily','Monthly','Quarterly'); fprintf(fileID, strcat(repmat('-',1,l),'\n'));
SubTableDisplay(power,iH,'bet0=0.95','bet0=0.90',l, fileID);

fprintf(fileID, strcat(repmat('-',1,l),'\n'));
fprintf(strcat(repmat('%12s',1,5),'\n'), 'G-test', '', 'Daily','Monthly','Quarterly');
fprintf(fileID, strcat(repmat('%12s',1,5),'\n'), 'G-test', '', 'Daily','Monthly','Quarterly'); fprintf(fileID, strcat(repmat('-',1,l),'\n'));
SubTableDisplay(power,iG,'bet0=0.990','bet0=0.985',l, fileID);
fprintf(fileID, strcat(repmat('-',1,l),'\n'));
end

% Auxiliary function 2
function SubTableDisplay(power,i,alt1,alt2,l, fileID)
iad = 1; icrt = 2;
disp(repmat('-',1,l));
data = power{i,1,iad}; % 1st alternative, AD
fprintf(strcat(repmat('%12s',1,2),repmat('%12.4f',1,3),'\n'), 'AD     ', alt1, data);
fprintf(strcat(repmat('%12s',1,2),repmat('%12.4f',1,2),'%12s','\n'), '', '',data(1)/data(3)-1,data(2)/data(3)-1,'--');
fprintf(fileID, strcat(repmat('%12s',1,2),repmat('%12.4f',1,3),'\n'), 'AD     ', alt1, data);
fprintf(fileID, strcat(repmat('%12s',1,2),repmat('%12.4f',1,2),'%12s','\n'), '', '',data(1)/data(3)-1,data(2)/data(3)-1,'--');

data = power{i,2,iad}; % 2nd alternative, AD
fprintf(strcat(repmat('%12s',1,2),repmat('%12.4f',1,3),'\n'), '', alt2, data);
fprintf(strcat(repmat('%12s',1,2),repmat('%12.4f',1,2),'%12s','\n'), '', '',data(1)/data(3)-1,data(2)/data(3)-1,'--');
fprintf(fileID, strcat(repmat('%12s',1,2),repmat('%12.4f',1,3),'\n'), '', alt2, data);
fprintf(fileID, strcat(repmat('%12s',1,2),repmat('%12.4f',1,2),'%12s','\n'), '', '',data(1)/data(3)-1,data(2)/data(3)-1,'--');

disp(repmat('-',1,l));
data = power{i,1,icrt}; % 1st alternative, CRT
fprintf(strcat(repmat('%12s',1,2),repmat('%12.4f',1,3),'\n'), 'CRT    ', alt1, data);
fprintf(strcat(repmat('%12s',1,2),repmat('%12.4f',1,2),'%12s','\n'), '', '',data(1)/data(3)-1,data(2)/data(3)-1,'--');
fprintf(fileID, strcat(repmat('%12s',1,2),repmat('%12.4f',1,3),'\n'), 'CRT    ', alt1, data);
fprintf(fileID, strcat(repmat('%12s',1,2),repmat('%12.4f',1,2),'%12s','\n'), '', '',data(1)/data(3)-1,data(2)/data(3)-1,'--');

data = power{i,2,icrt}; % 2nd alternative, CRT
fprintf(strcat(repmat('%12s',1,2),repmat('%12.4f',1,3),'\n'), '', alt2, data);
fprintf(strcat(repmat('%12s',1,2),repmat('%12.4f',1,2),'%12s','\n'), '', '',data(1)/data(3)-1,data(2)/data(3)-1,'--');
fprintf(fileID, strcat(repmat('%12s',1,2),repmat('%12.4f',1,3),'\n'), '', alt2, data);
fprintf(fileID, strcat(repmat('%12s',1,2),repmat('%12.4f',1,2),'%12s','\n'), '', '',data(1)/data(3)-1,data(2)/data(3)-1,'--');
disp(repmat('-',1,l));
end
