clear; clc; close all
% -------------------------------------------------------------------------
% Purpose: Simulation for stationary and cointegration models presented
%          in Section 7: Simulation of Chang, Lu and Park (2024)
% -------------------------------------------------------------------------
%          Results will be used to generate 
%              Figure 6 in Section 7 and Figures 7,8,9 in Appendix D
% -------------------------------------------------------------------------
% Saved results in '../' include:
%          (1) Mean test statistics (63 frequencies, 4 bandwidths)
%              'Section_7_Htest_stat_T30.csv'
%              'Section_7_Htest_stat_T50.csv'
%              'Section_7_Htest_stat_Double.csv'  
%              'Section_7_Gtest_stat_T30.csv'
%              'Section_7_Gtest_stat_T50.csv'
%              'Section_7_Gtest_stat_Double.csv'  
%          (2) Rejection probabilities (63 frequencies, 4 bandwidths)
%              'Section_7_Htest_rej_T30.csv'
%              'Section_7_Htest_rej_T50.csv'
%              'Section_7_Htest_rej_Double.csv'  
%              'Section_7_Gtest_rej_T30.csv'
%              'Section_7_Gtest_rej_T50.csv'
%              'Section_7_Gtest_rej_Double.csv' 
%          (3) Instability of the test results (61 frequencies, 4 bandwidths)
%              'Section_7_Htest_instability_T30.csv'
%              'Section_7_Htest_instability_T50.csv'
%              'Section_7_Gtest_instability_T30.csv'
%              'Section_7_Gtest_instability_T50.csv'
% -------------------------------------------------------------------------
% Call functions:
%   ../../Functions/MC/OU_est.m:  estimate parameters in an OU model
%   ../../Functions/MC/data_sim_main.m: simulate data according to MC models in Section 7
%   ../../Functions/test/myols.m: OLS estimation
% -------------------------------------------------------------------------
% Runtime: 6.42 hours on MacBook Air with M2 chip
% -------------------------------------------------------------------------
% This version: September 2024 by Ye Lu (ye.lu1@sydney.edu.au)
% -------------------------------------------------------------------------

savefolder = '../'; prefix = 'Section_7'; tspan = {'T30','T50','Double'};
addpath(genpath('../../Functions/'));

niter = 3000;         % number of simulation iterations
Tlist = [30 50];      % time span considered (in years)
nT    = length(Tlist);
del   = (1:63)/252;   % length of sampling intervals considered (in years)
ndel  = length(del);

% -------------------------------------------------------------------------
% Null hypothesis
% -------------------------------------------------------------------------
truealp = 0; truebet = 1; R = eye(2); r = [truealp; truebet];
coeff   = [truealp truebet]; % assume H0 is true in DGP

% -------------------------------------------------------------------------
% Set true parameters in the STATIONARY regression model
%     (X and U: stationary Ornstein-Uhlenbeck processes)
% by fitting the stationary simulation model to empirical Model II
% -------------------------------------------------------------------------
M2 = readtimetable('../../Data/M2_TBEU.csv');
y  = M2.Eudollar; x = M2.TB3m; [~,~,~,~,~,~,u] = myols(y,x,1);
[para.stat.sigx, para.stat.kapx] = OU_est(x, 1/252);  % parameters in process X
[para.stat.sigu, para.stat.kapu] = OU_est(u, 1/252);  % parameters in process U

% -------------------------------------------------------------------------
% Set true parameters in the COINTEGRATING regression model
% by fitting the nonstationary simulation model to empirical Model III
% -------------------------------------------------------------------------
M3 = readtimetable('../../Data/M3_logFX.csv'); %
y = M3.Frwd3m; x = M3.Spot; [~,~,~,~,~,~,u] = myols(y,x,1);
d = diff(x); n = numel(x); T3 = numel(x)/252;
para.coin.sigx = sqrt(d'*d/T3);                      % parameters in process X
[para.coin.sigu, para.coin.kapu] = OU_est(u, 1/252); % parameters in process U

tic
% -------------------------------------------------------------------------
% Stationary regression with H-tests
% -------------------------------------------------------------------------
mdl = 'stat';
ad  = zeros(niter,ndel,nT); nw = zeros(size(ad)); rt = zeros(size(ad)); crt = zeros(size(ad));
ad2 = zeros(niter,ndel); nw2 = zeros(size(ad2)); rt2 = zeros(size(ad2)); crt2 = zeros(size(ad2));

for m = 1 : niter
    fprintf('%s, iter=%d\n', mdl, m); rng(m)
    for k = 1:nT
        T = Tlist(k); [y, x] = data_sim_main(T, del(1), para, mdl, coeff); % simulate daily data
        for i = 1:ndel
            delta = del(i);
            % Test statitics (time span fixed)
            [~,Hstat] = Htest(delta,y(i:i:end),x(i:i:end),1,R,r);
            ad(m,i,k) = Hstat.ad91; nw(m,i,k) = Hstat.nw94;
            rt(m,i,k) = Hstat.rt;  crt(m,i,k) = Hstat.crt;
            if T == 50
                % Test statistics (time span varies with delta)
                n = round(3 * delta^(-1/2)/delta);
                [~, Hstat2] = Htest(delta,y(i:i:i*n),x(i:i:i*n),1,R,r);
                ad2(m,i)    = Hstat2.ad91; nw2(m,i) = Hstat2.nw94;
                rt2(m,i)    = Hstat2.rt; crt2(m,i)  = Hstat2.crt;
            end
        end
    end
end
% -- Parse and save the simulation results --
tst = 'H'; savename = cell(8,1);
for i = 1:3
    savename{i}   = sprintf('%s%s_%stest_stat_%s.csv',savefolder,prefix,tst,tspan{i});
    savename{i+3} = sprintf('%s%s_%stest_rej_%s.csv', savefolder,prefix,tst,tspan{i});
    if i ~=3
        savename{i+6} = sprintf('%s%s_%stest_instability_%s.csv',savefolder,prefix,tst,tspan{i});
    end
end
SaveResults(ad, ad2, nw, nw2, rt, rt2, crt, crt2, niter, savename)

% -------------------------------------------------------------------------
% Cointegrating regression with G-tests
% -------------------------------------------------------------------------
mdl = 'coin';
ad  = zeros(niter,ndel,nT); nw = zeros(size(ad)); rt = zeros(size(ad)); crt = zeros(size(ad));
ad2 = zeros(niter,ndel); nw2 = zeros(size(ad2)); rt2 = zeros(size(ad2)); crt2 = zeros(size(ad2));

for m = 1 : niter
    fprintf('%s, iter=%d\n', mdl, m); rng(m)
    for k = 1:nT
        T = Tlist(k); [y, x] = data_sim_main(T, del(1), para, mdl, coeff); % simulate daily data
        for i = 1:ndel
            delta = del(i);
            % Test statitics (time span fixed)
            [~,Gstat] = Gtest(delta,y(i:i:end),x(i:i:end),1,R,r);
            ad(m,i,k) = Gstat.ad91; nw(m,i,k) = Gstat.nw94;
            rt(m,i,k) = Gstat.rt;  crt(m,i,k) = Gstat.crt;
            if T == 50
                % Test statistics (time span varies with delta)
                n = round(3 * delta^(-1/2)/delta);
                [~, Gstat2] = Gtest(delta,y(i:i:i*n),x(i:i:i*n),1,R,r);
                ad2(m,i) = Gstat2.ad91; nw2(m,i) = Gstat2.nw94;
                rt2(m,i) = Gstat2.rt; crt2(m,i)  = Gstat2.crt;
            end
        end
    end
end
% -- PARSE and save the simulation results --
tst = 'G'; 
for i = 1:3
    savename{i}   = sprintf('%s%s_%stest_stat_%s.csv',savefolder,prefix,tst,tspan{i});
    savename{i+3} = sprintf('%s%s_%stest_rej_%s.csv', savefolder,prefix,tst,tspan{i});
    if i ~=3
        savename{i+6} = sprintf('%s%s_%stest_instability_%s.csv',savefolder,prefix,tst,tspan{i});
    end
end
SaveResults(ad, ad2, nw, nw2, rt, rt2, crt, crt2, niter, savename)

toc
rmpath(genpath('../../Functions/'));

%% Parse and save the simulation results

function SaveResults(ad, ad2, nw, nw2, rt, rt2, crt, crt2, niter, savename)
% -------------------------------------------------------------------------
% Purpose: This auxiliary function helps save
% (i)   Mean statistics
% (ii)  Rejection probabilities
% (iii) Percentages of instability of test restuls
%  for
% (a) 4 bandwidths
% (b) 63 frequencies, and
% (c) 3 different time spans (T=30, T=50, T under double asymptotics)
% -------------------------------------------------------------------------
ndel  = 63; % number of sampling frequencies

% objects to be saved
test_stat   = zeros(4,ndel,3);
rej_prob    = zeros(4,ndel,3);
instability = zeros(4,ndel,2);

% indices for bandwidths
iad = 1; inw = 2; irt = 3; icrt = 4;

% 5% critical value of the test for joint hypothesis of 2 restrictions
df = 2; cv = chi2inv(0.95,df);

% -- results for T=30 and T=50 (t = 1,2) --
for t = 1:2
    % mean test statistics
    test_stat(iad,:,t)  = mean(ad(:,:,t));
    test_stat(inw,:,t)  = mean(nw(:,:,t));
    test_stat(irt,:,t)  = mean(rt(:,:,t));
    test_stat(icrt,:,t) = mean(crt(:,:,t));

    % rejection probabilities
    rej_prob(iad,:,t)  = mean(ad(:,:,t)>cv);
    rej_prob(inw,:,t)  = mean(nw(:,:,t)>cv);
    rej_prob(irt,:,t)  = mean(rt(:,:,t)>cv);
    rej_prob(icrt,:,t) = mean(crt(:,:,t)>cv);

    % instability percentages (subject to 1-day change in sampling frequency)
    rej{iad} = ad(:,:,t)>cv; rej{icrt} = crt(:,:,t)>cv;
    rej{irt} = rt(:,:,t)>cv; rej{inw}  = nw(:,:,t)>cv;
    D = zeros(niter,63,4);
    for i = 2:62
        for k = 1:4
            D(:,i,k)=abs(rej{k}(:,i-1)-rej{k}(:,i))+abs(rej{k}(:,i)-rej{k}(:,i+1))>0;
        end
    end
    instability(iad,:,t)  = 100 * mean(D(:,:,iad));
    instability(inw,:,t)  = 100 * mean(D(:,:,inw));
    instability(irt,:,t)  = 100 * mean(D(:,:,irt));
    instability(icrt,:,t) = 100 * mean(D(:,:,icrt));
end

% -- results for double asymptotics (t=3) --
test_stat(iad,:,3)  = mean(ad2);
test_stat(inw,:,3)  = mean(nw2);
test_stat(irt,:,3)  = mean(rt2);
test_stat(icrt,:,3) = mean(crt2);

rej_prob(iad,:,3)  = mean(ad2 > cv);
rej_prob(inw,:,3)  = mean(nw2 > cv);
rej_prob(irt,:,3)  = mean(rt2 > cv);
rej_prob(icrt,:,3) = mean(crt2 > cv);

% -- write to csv files --
for t  = 1:3
    writematrix(test_stat(:,:,t), savename{t})
    writematrix(rej_prob(:,:,t), savename{t+3})
    if t~=3
        writematrix(instability(:,:,t), savename{t+6})
    end
end

end