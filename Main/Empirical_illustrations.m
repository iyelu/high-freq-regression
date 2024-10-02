clc; clear; close all;
% -------------------------------------------------------------------------
% Purpose: Replicating figures and tables in Sections 2, 4, 5, 6, and 
%          Appendices B,C of the paper 
%          "Understanding regressions with observations collected at 
%           high-frequency over long span" by Chang, Lu and Park (2024)
% -------------------------------------------------------------------------
% Outputs:
%     Figures saved in '../Figures':
%        Figure1: Data plots for Models I-IV (Section 2)
%        Figure2: Wald tests in Models I-IV (Section 2)
%        Figure3: Estimated residual AR coefficients in Models I-IV (Section 4)
%        Figure4: Robust Wald tests in Models I-IV (Section 5)
%        Figure5: Short and long rates before and after the observed Greenspan 
%                 conundrum at quarterly, monthly, and daily frequencies (Section 6)
%     Tables saved in '../Tables':
%       Tables1-3: Stationary test results for daily, monthly, quarterly regressions 
%                   with pre and post-Greenspan conundrum samples (Section 6)
%       Tables4-6: Cointegration test results for daily, monthly, quarterly regressions 
%                   with pre and post-Greenspan conundrum samples (Appendix B)
%       Tables7-9: Stationary test results for daily, monthly, quarterly regressions 
%                   with pre and post-financial crisis samples (Appendix C.1)
%       Tables10-12: Cointegration test results for daily, monthly, quarterly 
%                   regressions with pre and post-financial crisis samples (Appendix C.2)
% -------------------------------------------------------------------------
% This version: September 2024 by Ye Lu (ye.lu1@sydney.edu.au)
% -------------------------------------------------------------------------

tic
addpath(genpath('../Functions/')); % add path to functions to be called
% -------------------------------------------------------------------------
% Setting of figures
% -------------------------------------------------------------------------
set(0,'DefaultTextInterpreter','latex'); set(0,'DefaultAxesTickLabelInterpreter','latex'); set(0,'DefaultLegendInterpreter', 'latex')
set(0,'defaultAxesXGrid','on'); set(0,'defaultAxesYGrid','on'); set(0,'defaultAxesFontSize',13); set(0,'DefaultLineLineWidth',1)
x0=0; y0=0; width=800; height=500; xl=0.25; titles = {'Model I', 'Model II', 'Model III', 'Model IV'}; 
savefigfolder = '../Figures';

% -------------------------------------------------------------------------
% Call function "data_clearning.m" in '../Functions/' to produce
%  4 datasets saved in timetables for four empirical models: Models I-IV
% -------------------------------------------------------------------------
dir_rawData = '../Data/'; % directory where raw data are saved
dir_output = './';        % directory where cleaned datasets are saved
output = data_cleaning(dir_rawData, dir_output); % 4 timetables produced
M1 = output.M1; M2 = output.M2; M3 = output.M3; M4 = output.M4;

% -------------------------------------------------------------------------
% Testing 3 null hypotheses: 
%   $H_0:\beta_0=0$ (single hypothesis for intercept)
%   $H_0:\beta_1=1$ (single hypothesis for slope coefficient)
%   $H_0:\beta_0=0 and \beta_1=1$ (joint hypothesis for both intercept and slope)
% using data sampled from varying frequencies (daily to quarterly) for
%   Non-robust test: the Wald test with non-robust variance estimator
%   Robust tests with a long-run variance estimator which adopts, respectively 
%                          Andrews bandwidth (AD)
%                          Newey-West bandwidth (NW)
%                          Discrete-time rule of thumb (RT)
%                          Continuous-time rule of thumb (CRT)
% -------------------------------------------------------------------------
truebet0 = 0; truebet1 = 1; % True parameters in the null hypotheses
delta = (1:63)/252;  % length of sampling intervals from daily (1/252) to quarterly (1/4) frequencies
nd  = numel(delta);  % number of frequencies considered
 
% Initicating matrices for saving calculated test statistics:
wd  = zeros(4,nd,3); % non-robust test statistic with non-robust variance estimator
ad  = zeros(4,nd,3); % robust test statistic with Andrews (AD) bandwidth
nw  = zeros(4,nd,3); % robust test statistic with Newey-West (NW) bandwidth
rt  = zeros(4,nd,3); % robust test statistic with discrete-time rule of thumb (RT) bandwidth
crt = zeros(4,nd,3); % robust test statistic with continuout-time rule of thumb (CRT) bandwidth
 
% Initicating matrices for saving $\hat\rho$ for Figure 3:
rhohat = zeros(4,nd);
 
% -------------------------------------------------------------------------
% Discrete regressions with data sampled at varying frequencies for  
%   4 models
%   testing 3 null hypotheses with 5 test statistics
%   samping frequency ranging from daily to quarterly
% -------------------------------------------------------------------------
for mdl = 1:4   % 4 models
    switch mdl
        case 1  % Model I
            y = M1.TB10y; x = M1.TB3m;
        case 2  % Model II
            y = M2.Eudollar; x = M2.TB3m;
        case 3  % Model III
            y = M3.Frwd3m; x = M3.Spot;
        case 4  % Model IV
            y = M4.SP500futures; x = M4.SP500;
    end
    % linear regression with discrete data y = bet0 + bet1 * x + u
    for i = 1:nd
        lhs = y(i:i:end); rhs = [ones(size(lhs)), x(i:i:end)];
        bet = rhs\lhs; uhat = lhs - rhs * bet;
        rhohat(mdl,i) = uhat(1:end-1)\uhat(2:end);
    end
    % hypothesis testing using non-robust and robust Wald tests
    for j = 1:3  % 3 null hypotheses
        switch j
            case 1
                R = [1 0]; r = truebet0;  % H0: beta0 = 0
            case 2
                R = [0 1]; r = truebet1;  % H0: beta1 = 1
            case 3
                R = eye(2); r = [truebet0; truebet1]; % H0: beta0 = 0 and beta1 = 1
        end
        for i = 1:nd  % varing sampling frequencies
            [~,Gstat,Fstat] = Gtest(i/252,y(i:i:end),x(i:i:end),1,R,r); [~,Hstat] = Htest(i/252,y(i:i:end),x(i:i:end),1,R,r); wd(mdl,i,j) = Fstat; 
            if mdl < 3  % Model I and II: H-test
                ad(mdl,i,j) = Hstat.ad91; nw(mdl,i,j) = Hstat.nw94; rt(mdl,i,j) = Hstat.rt; crt(mdl,i,j) = Hstat.crt;
            end
            if mdl > 2  % Model III and IV: G-test
                ad(mdl, i,j) = Gstat.ad91; nw(mdl, i,j) = Gstat.nw94; rt(mdl, i,j) = Gstat.rt; crt(mdl,i,j) = Gstat.crt;
            end
        end
    end
end
%% Figure 1:
figure; set(gcf, 'position', [x0,y0,width,height]);
subplot(2,2,1); box on; hold on; plot(M1.Time,M1.TB3m); plot(M1.Time,M1.TB10y); hold off
legend('3-month T-bill rate', '10-year T-bond rate'); ylabel('percent');
subplot(2,2,2); box on; hold on; plot(M2.Time,M2.TB3m); plot(M2.Time,M2.Eudollar); hold off
legend('3-month T-bill rate', '3-month Eurodollar rate'); ylabel('percent');
subplot(2,2,3); box on; hold on; plot(M3.Time, M3.Spot); plot(M3.Time, M3.Frwd3m); hold off
legend('Spot US/UK', 'Forward US/UK'); ylabel('log rate');
subplot(2,2,4); box on; hold on; plot(M4.Time, M4.SP500); plot(M4.Time, M4.SP500futures); hold off
legend('S\&P 500', 'S\&P 500 futures', 'location', 'northwest'); ylabel('log price')
saveas(gcf, sprintf('%s/Figure1.png', savefigfolder)); % saveas(gcf, sprintf('%s/Figure1', savefigfolder), 'epsc')

%% Figure 2:
figure; set(gcf, 'position', [x0, y0, width, height]);
for mdl = 1:4
    subplot(2,2,mdl); box on; plot(delta, wd(mdl,:,3), '-', 'linewidth', 1.5); xlim([0 xl]); xlabel('$\delta$'); ylabel('test statistic');
end
saveas(gcf, sprintf('%s/Figure2.png', savefigfolder)); % saveas(gcf, sprintf('%s/Figure2', savefigfolder), 'epsc')

%% Figure 3:
figure; set(gcf, 'position', [x0, y0, width, height])
for mdl = 1:4
    subplot(2,2,mdl); box on
    plot(delta, rhohat(mdl,:), '-', 'linewidth', 1.5)
    xlim([0 xl]); ylim([0 1]); xlabel('$\delta$'); ylabel('$\tilde\rho$','rotation',0);
end
saveas(gcf, sprintf('%s/Figure3.png', savefigfolder)); % saveas(gcf, sprintf('%s/Figure3', savefigfolder), 'epsc')

%% Figure 4:
figure; set(gcf, 'position', [x0, y0, width, height]); yl = [9000 5000 1000 60];
for mdl=1:4
    subplot(2,2,mdl); box on; hold on
    plot(delta, nw(mdl,:,3), '-.', 'linewidth', 1.5)
    plot(delta, ad(mdl,:,3), '--', 'linewidth', 1.5)
    plot(delta, rt(mdl,:,3),  ':', 'linewidth', 2.3)
    plot(delta, crt(mdl,:,3), '-', 'linewidth', 1.5)
    hold off
    if mdl < 3
        legend('$H$-test: NW', '$H$-test: AD', '$H$-test: RT', '$H$-test: CRT')
    else
        legend('$G$-test: NW', '$G$-test: AD', '$G$-test: RT', '$G$-test: CRT')
    end
    xlim([0 xl]); ylim([0 yl(mdl)]); xlabel('$\delta$'); ylabel('test statistic');
end
saveas(gcf, sprintf('%s/Figure4.png', savefigfolder)); % saveas(gcf, sprintf('%s/Figure4', savefigfolder), 'epsc')

%% Pre and Post-Greenspan conundrum samples *(before and after January 1, 2001)
S1 = timerange(M1.Time(1),'2000-Dec-31','closed'); S2 = timerange('2001-Jan-1', M1.Time(end),'closed'); nS = 2;
data = cell(nS,1);
% Daily subsamples
data{1}.daily.time = M1.Time(S1,:); data{1}.daily.y = M1.TB10y(S1,:); data{1}.daily.x = M1.TB3m(S1,:);
data{2}.daily.time = M1.Time(S2,:); data{2}.daily.y = M1.TB10y(S2,:); data{2}.daily.x = M1.TB3m(S2,:);
% Monthly subsamples
for S = 1:nS
    ind = [diff(day(data{S}.daily.time))<0; true]; % find end of the month index and manually add the last day
    data{S}.monthly.time = data{S}.daily.time(ind); data{S}.monthly.y = data{S}.daily.y(ind); data{S}.monthly.x = data{S}.daily.x(ind);
end
% Quarterly subsamples:
for S = 1:nS
    ind = [diff(quarter(data{S}.daily.time))~=0; true]; % find end of the quarter index and manually add the last day
    data{S}.quarterly.time = data{S}.daily.time(ind); data{S}.quarterly.y = data{S}.daily.y(ind); data{S}.quarterly.x = data{S}.daily.x(ind);
end

%% Figure 5: 
set(0,'DefaultLineLineWidth',1.2); x0=0; y0=0; width=800; height=600;
figure; set(gcf, 'position', [x0,y0,width,height]);
for S = 1:nS
    subplot(3,2,S)
    box on; hold on; plot(data{S}.quarterly.time, data{S}.quarterly.x); plot(data{S}.quarterly.time, data{S}.quarterly.y); hold off
    ylabel('percent'); if S ==1; ylim([0 20]); else; ylim([0 8]); end; legend('Short rate', 'Long rate')
    subplot(3,2,S+2)
    box on; hold on; plot(data{S}.monthly.time, data{S}.monthly.x); plot(data{S}.monthly.time, data{S}.monthly.y); hold off
    ylabel('percent'); if S ==1; ylim([0 20]); else; ylim([0 8]); end; legend('Short rate', 'Long rate')
    subplot(3,2,S+4)
    box on; hold on; plot(data{S}.daily.time, data{S}.daily.x); plot(data{S}.daily.time, data{S}.daily.y); hold off
    ylabel('percent'); if S ==1; ylim([0 20]); else; ylim([0 8]); end; legend('Short rate', 'Long rate')
end
saveas(gcf, sprintf('%s/Figure5.png', savefigfolder)); % saveas(gcf, sprintf('%s/Figure5', savefigfolder), 'epsc')

% True parameters in three null hypotheses considered:
truealp = 0; truebet = 1;
Hypo = {'alp=0', 'bet=1', 'alp=0, bet=1'}; nh = length(Hypo);                         % to be printed on screen
Hypo_txt = {'$\mbox H_0^\alpha$', '$\mbox H_0^\beta$', '$\mbox H_0^{\alpha,\beta}$'}; % to be printed to txt files
% Length of sampling interval for quarterly, monthly and daily data, respectively 
delta_list = [1/4 21/252 1/252];

filenames_H = {'../Tables/Table1.txt', '../Tables/Table2.txt', '../Tables/Table3.txt'};
filenames_G = {'../Tables/Table4.txt', '../Tables/Table5.txt', '../Tables/Table6.txt'};
Caps_H = {strcat(repmat(':',1,17),' TABLE 1: Quarterly regression, H-test ',repmat(':',1,17)),...
    strcat(repmat(':',1,17),' TABLE 2: Monthly regression, H-test ',repmat(':',1,17)),...
    strcat(repmat(':',1,18),' TABLE 3: Daily regression, H-test ',repmat(':',1,18))};
Caps_G = {strcat(repmat(':',1,17),' TABLE 4: Quarterly regression, G-test ',repmat(':',1,17)),...
    strcat(repmat(':',1,17),' TABLE 5: Monthly regression, G-test ',repmat(':',1,17)),...
    strcat(repmat(':',1,18),' TABLE 6: Daily regression, G-test ',repmat(':',1,18))};

data{1} = M1(S1,:); data{2} = M1(S2,:);

for d = 1:3  % three sampling frequencies
    delta = delta_list(d); wd  = zeros(nh,nS);   % non-robust Wald test
    G_ad  = zeros(3,2); G_nw  = zeros(3,2); G_rt  = zeros(3,2); G_crt = zeros(3,2); % robust G-test
    H_ad  = zeros(3,2); H_nw  = zeros(3,2); H_rt  = zeros(3,2); H_crt = zeros(3,2); % robust H-test

    for S = 1:nS % loop for subsamples
        switch d
            case 1 % quarterly
                ind = [diff(quarter(data{S}.Time))~=0; true]; y = data{S}.TB10y(ind); x = data{S}.TB3m(ind);
            case 2 % monthly
                ind = [diff(day(data{S}.Time))<0; true]; y = data{S}.TB10y(ind); x = data{S}.TB3m(ind);
            case 3 % daily
                y = data{S}.TB10y; x = data{S}.TB3m;
        end
        for j = 1:nh  % loop for null hypotheses: H0_alp, H0_bet, H0_alp_bet
            switch j
                case 1
                    R = [1 0]; r = truealp;
                case 2
                    R = [0 1]; r = truebet;
                case 3
                    R = eye(2); r = [truealp; truebet];
            end
            [~,Gstat,Fstat] = Gtest(delta,y,x,1,R,r); [~,Hstat] = Htest(delta,y,x,1,R,r); wd(j,S)   = Fstat;
            G_ad(j,S) = Gstat.ad91; G_nw(j,S) = Gstat.nw94; G_rt(j,S) = Gstat.rt; G_crt(j,S) = Gstat.crt;
            H_ad(j,S) = Hstat.ad91; H_nw(j,S) = Hstat.nw94; H_rt(j,S) = Hstat.rt; H_crt(j,S) = Hstat.crt;
        end
    end
    %% Tables 1-3:
    disp(Caps_H{d}); fileID = fopen(filenames_H{d},'w'); TableDisplay(data, Hypo, Hypo_txt, wd, H_ad, H_crt, H_nw, H_rt, fileID); fclose(fileID);
    %% Tables 4-6: 
    disp(Caps_G{d}); fileID = fopen(filenames_G{d},'w'); TableDisplay(data, Hypo, Hypo_txt, wd, G_ad, G_crt, G_nw, G_rt, fileID); fclose(fileID);
end

%% Pre and post financial crisis samples *(before and after January 1, 2008)
S1 = timerange(M1.Time(1),'2007-Dec-31','closed'); S2 = timerange('2008-Jan-1', M1.Time(end),'closed'); nS = 2; data{1} = M1(S1,:); data{2} = M1(S2,:);
filenames_H = {'../Tables/Table7.txt', '../Tables/Table8.txt', '../Tables/Table9.txt'};
filenames_G = {'../Tables/Table10.txt', '../Tables/Table11.txt', '../Tables/Table12.txt'};
Caps_H = {strcat(repmat(':',1,17),' TABLE 7: Quarterly regression, H-test ',repmat(':',1,17)),...
    strcat(repmat(':',1,17),' TABLE 8: Monthly regression, H-test ',repmat(':',1,17)),...
    strcat(repmat(':',1,18),' TABLE 9: Daily regression, H-test ',repmat(':',1,18))};
Caps_G = {strcat(repmat(':',1,17),' TABLE 10: Quarterly regression, G-test ',repmat(':',1,17)),...
    strcat(repmat(':',1,17),' TABLE 11: Monthly regression, G-test ',repmat(':',1,17)),...
    strcat(repmat(':',1,18),' TABLE 12: Daily regression, G-test ',repmat(':',1,18))};
l = 72; % length of the table

for d = 1:3  % three sampling frequencies
    delta  = delta_list(d); wd  = zeros(nh,nS);   % non-robust Wald test
    G_ad  = zeros(3,2); G_nw  = zeros(3,2); G_rt  = zeros(3,2); G_crt = zeros(3,2); % robust G-test
    H_ad  = zeros(3,2); H_nw  = zeros(3,2); H_rt  = zeros(3,2); H_crt = zeros(3,2); % robust H-test

    for S = 1:nS % loop for subsamples
        switch d
            case 1  % quarterly
                ind = [diff(quarter(data{S}.Time))~=0; true]; y = data{S}.TB10y(ind); x = data{S}.TB3m(ind);
            case 2  % monthly
                ind = [diff(day(data{S}.Time))<0; true]; y = data{S}.TB10y(ind); x = data{S}.TB3m(ind);
            case 3  % daily
                y = data{S}.TB10y; x = data{S}.TB3m;
        end
        for j = 1:nh  % loop for null hypotheses: H0_alp, H0_bet, H0_alp_bet
            switch j
                case 1
                    R = [1 0]; r = truealp;
                case 2
                    R = [0 1]; r = truebet;
                case 3
                    R = eye(2); r = [truealp; truebet];
            end
            [~,Gstat,Fstat] = Gtest(delta,y,x,1,R,r); [~,Hstat] = Htest(delta,y,x,1,R,r); wd(j,S)   = Fstat;
            G_ad(j,S) = Gstat.ad91; G_nw(j,S) = Gstat.nw94; G_rt(j,S) = Gstat.rt; G_crt(j,S) = Gstat.crt;
            H_ad(j,S) = Hstat.ad91; H_nw(j,S) = Hstat.nw94; H_rt(j,S) = Hstat.rt; H_crt(j,S) = Hstat.crt;
        end
    end
    %% Tables 7-9:
    disp(Caps_H{d}); fileID = fopen(filenames_H{d},'w'); TableDisplay(data, Hypo, Hypo_txt, wd, H_ad, H_crt, H_nw, H_rt, fileID); fclose(fileID);
    %% Tables 10-12:
    disp(Caps_G{d}); fileID = fopen(filenames_G{d},'w'); TableDisplay(data, Hypo, Hypo_txt, wd, G_ad, G_crt, G_nw, G_rt, fileID); fclose(fileID);
end

rmpath(genpath('./Functions/')); % remove path to auxiliary functions
toc

%% Auxiliary function for displaying tables

function TableDisplay(data, Hypo, Hypo_txt, wd, ad, crt, nw, rt, fileID)
l = 72; % length of the table
% -- First subsample --
disp(repmat('-',1,l)); fprintf('First subsample: %d-%d\n',year(data{1}.Time(1)), year(data{1}.Time(end)));
S = 1; disp(repmat('-',1,l));
fprintf(strcat(repmat('%12s',1,6),'\n'), 'H0', 'Wald', 'AD','CRT','NW','RT'); disp(repmat('-',1,l));
for j = 1:3  % 3 hypotheses
    % print on screen
    fprintf(strcat('%12s',repmat('%12.2f',1,5),'\n'), Hypo{j},wd(j,S),ad(j,S),crt(j,S),nw(j,S),rt(j,S));
    % print to txt file
    fprintf(fileID,strcat('%25s',repmat(' & $%8.2f$',1,5),'\\\\ \n'), Hypo_txt{j},wd(j,S),ad(j,S),crt(j,S),nw(j,S),rt(j,S));
    if j == 2  % print p-values
        df = 1;
        % print on screen
        fprintf(strcat('%12s%12s',repmat('%12.4f',1,4),'\n'),...
            '', '', 1-chi2cdf(ad(j,S),df), 1-chi2cdf(crt(j,S),df), 1-chi2cdf(nw(j,S),df), 1-chi2cdf(rt(j,S),df));
        % print to txt file
        fprintf(fileID, strcat('%25s%12s',repmat(' & $[%6.4f]$',1,4),'\\\\ \n'),...
            '', '&', 1-chi2cdf(ad(j,S),df), 1-chi2cdf(crt(j,S),df), 1-chi2cdf(nw(j,S),df), 1-chi2cdf(rt(j,S),df));
    end
end
% -- Second subsample --
disp(repmat('-',1,l)); fprintf('Second subsample: %d-%d\n',year(data{2}.Time(1)), year(data{2}.Time(end)));
S = 2; disp(repmat('-',1,l));
fprintf(strcat(repmat('%12s',1,6),'\n'), 'H0', 'Wald', 'AD','CRT','NW','RT'); disp(repmat('-',1,l));
fprintf(fileID,'\n'); % print a ling break in txt file
for j = 1:3 % 3 hypotheses
    % print on screen
    fprintf(strcat('%12s',repmat('%12.2f',1,5),'\n'), Hypo{j},wd(j,S),ad(j,S),crt(j,S),nw(j,S),rt(j,S));
    % print to txt file
    fprintf(fileID,strcat('%25s',repmat(' & $%8.2f$',1,5),'\\\\ \n'), Hypo_txt{j},wd(j,S),ad(j,S),crt(j,S),nw(j,S),rt(j,S));
end
disp(repmat('-',1,l)); fprintf(' \n \n ');
end