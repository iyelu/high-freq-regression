clc;clear;close all
% -------------------------------------------------------------------------
% Purpose: Simulation for stationary and cointegration models presented
%          in Appendix E: Additional Simulation Results of Chang, Lu and
%          Park (2024)
% -------------------------------------------------------------------------
%          Results will be used to generate Figures 10,11 in Appendix E
% -------------------------------------------------------------------------
% Saved results in '../' include:
%          (1) Mean test statistics (63 frequencies, 4 bandwidths)
%              'Appendix_E_Htest_stat_T30.csv'
%              'Appendix_E_Htest_stat_T50.csv'
%              'Appendix_E_Gtest_stat_T30.csv'
%              'Appendix_E_Gtest_stat_T50.csv'
%          (2) Rejection probabilities (63 frequencies, 4 bandwidths)
%              'Appendix_E_Htest_rej_T30.csv'
%              'Appendix_E_Htest_rej_T50.csv'
%              'Appendix_E_Gtest_rej_T30.csv'
%              'Appendix_E_Gtest_rej_T50.csv'
% -------------------------------------------------------------------------
% Call functions:
%   ../../Functions/MC/OU_est.m: estimate parameters in OU model 
%   ../../Functions/MC/SR_mle.m: estimate parameters in Feller's square root model   
%   ../../Functions/data_sim_additional.m: simulate data according to MC models in Appendix E
%   ../../Functions/test/myols.m
% -------------------------------------------------------------------------
% Runtime: 8.31 hours on iMac
% -------------------------------------------------------------------------
% This version: September 2024 by Ye Lu (ye.lu1@sydney.edu.au)
% -------------------------------------------------------------------------

savefolder = '../'; prefix = 'Appendix_E'; tspan = {'T30','T50'};
addpath(genpath('../../Functions/'));

niter = 3000;        % number of simulation iterations
Tlist = [30 50];     % time span considered (in years)
nT    = length(Tlist);
del   = (1:63)/252;  % length of sampling intervals considered (in years)
ndel  = length(del);

% ---------------------------------------------------------------------
% Null hypothesis
% ---------------------------------------------------------------------
truealp = 0; truebet = 1; R = eye(2); r = [truealp; truebet];
coeff   = [truealp truebet]; % assume H0 is true in DGP

% ---------------------------------------------------------------------
% Set true parameters in the STATIONARY regression model
%     (X: Feller's SR process | U: mean-zero Ornstein-Uhlenbeck process)
% by fitting the stationary simulation model to empirical Model II
% ---------------------------------------------------------------------
M2 = readtimetable('../../Data/M2_TBEU.csv');

% select data before 2007 to avoid 3-month T-bill rate hitting zero-lower-bound
S  = timerange(M2.Time(1), '31-Dec-2007', 'closed');  data = M2(S,:);
y = data.Eudollar; x = data.TB3m; [~,~,~,~,~,~,u] = myols(y,x,1);

% estimate parameters in SR and OU processes
[para.stat.kapx, para.stat.mux, para.stat.sigx] = SR_mle(x,1/252);  
[para.stat.sigu, para.stat.kapu]                = OU_est(u,1/252); 

% correlation between two Brownian motions of X and U
top = diff(x)'*diff(u); bot = (1/252)*para.stat.sigx*para.stat.sigu*sum(x); para.stat.rho = top/bot;

% ---------------------------------------------------------------------
% Set true parameters in the COINTEGRATING regression model
% ---------------------------------------------------------------------
%-- Model for X: Heston stochastic volatility model ------------------
% ---------------------------------------------------------------------
rf = 0.02; d = 0.015; % constant risk-free rate and S&P 500 dividend yield
kapx = 3; mux = 0.1; sigx = 0.25; rho = -0.8; lambda = 4;

% specifiy the drift function
C = [rf-d; kapx*mux];                    % constant term 
D = [0, lambda*(1-rho^2)-0.5; 0, -kapx]; % slope term 
para.coin.Xmu = @(x)(C + D * x);         % drift function

% specify the diffusion function (component-wise)
para.coin.Xsig11 = @(x)(sqrt((1-rho^2)*x(2)));
para.coin.Xsig12 = @(x)(rho*sqrt(x(2))); 
para.coin.Xsig21 = @(x)(0);
para.coin.Xsig22 = @(x)(sigx*sqrt(x(2)));

% initial values in the simulation
para.coin.Xinitial = [log(100); mux];

% ---------------------------------------------------------------------
% -- Model for U: GARCH stochastic volatility model -------------------
% ---------------------------------------------------------------------
kapu = 20; kapv = 2; muv = 0.0002; sigv = 0.008;

% specifiy the drift function
C = [0; kapv*muv];                % constant term 
D = [-kapu, 0; 0, -kapv];         % slope term 
para.coin.Umu = @(x)(C + D * x);  % drift function

% specify the diffusion function (component-wise)
para.coin.Usig11 = @(x)(sqrt(x(2)));
para.coin.Usig12 = @(x)(0); 
para.coin.Usig21 = @(x)(0);
para.coin.Usig22 = @(x)(sigv*sqrt(x(2)));

% initial values in the simulation
para.coin.Uinitial = [0; muv];

tic
% ---------------------------------------------------------------------
% Stationary regression with H-tests
% ---------------------------------------------------------------------
mdl = 'stat'; ad = zeros(niter,ndel,nT); nw = zeros(size(ad)); rt = zeros(size(ad)); crt = zeros(size(ad));
for m = 1 : niter
    fprintf('%s, iter=%d\n', mdl, m); rng(m)
    for k = 1:nT
        T = Tlist(k);
        [y, x] = data_sim_additional(T, del(1), para, mdl, coeff);  % simulate daily data
        for i = 1:ndel
            delta = del(i);
            [~,Hstat] = Htest(delta,y(i:i:end),x(i:i:end),1,R,r);
            ad(m,i,k) = Hstat.ad91; nw(m,i,k) = Hstat.nw94;
            rt(m,i,k) = Hstat.rt;  crt(m,i,k) = Hstat.crt;
        end
    end
end
% -- PARSE and save the simulation results --
tst = 'H'; savename = cell(4,1);
for i = 1:2
    savename{i}   = sprintf('%s%s_%stest_stat_%s.csv',savefolder,prefix,tst,tspan{i});
    savename{i+2} = sprintf('%s%s_%stest_rej_%s.csv', savefolder,prefix,tst,tspan{i});
end
SaveResults(ad, nw, rt, crt, savename)

% ---------------------------------------------------------------------
% Cointegrating regression with G-tests
% ---------------------------------------------------------------------
mdl = 'coin'; ad = zeros(niter,ndel,nT); nw = zeros(size(ad)); rt = zeros(size(ad)); crt = zeros(size(ad));
for m = 1 : niter
    fprintf('%s, iter=%d\n', mdl, m); rng(m)
    for k = 1:nT
        T = Tlist(k);
        [y, x] = data_sim_additional(T, del(1), para, mdl, coeff); % simulate daily data
        for i = 1:ndel
            delta = del(i);
            [~,Gstat] = Gtest(delta,y(i:i:end),x(i:i:end),1,R,r);
            ad(m,i,k) = Gstat.ad91; nw(m,i,k) = Gstat.nw94;
            rt(m,i,k) = Gstat.rt;  crt(m,i,k) = Gstat.crt;

        end
    end
end
% -- PARSE and save the simulation results --
tst = 'G'; savename = cell(4,1);
for i = 1:2
    savename{i}   = sprintf('%s%s_%stest_stat_%s.csv',savefolder,prefix,tst,tspan{i});
    savename{i+2} = sprintf('%s%s_%stest_rej_%s.csv', savefolder,prefix,tst,tspan{i});
end
SaveResults(ad, nw, rt, crt, savename)

toc

rmpath(genpath('../Functions/'));

%% Parse and save the simulation results

function SaveResults(ad, nw, rt, crt, savename)
% -------------------------------------------------------------------------
% Purpose: This auxiliary function helps save
% (i)   Mean statistics
% (ii)  Rejection probabilities
%  for
% (a) 4 bandwidths
% (b) 63 frequencies, and
% (c) 2 different time spans (T=30, T=50)
% -------------------------------------------------------------------------
ndel  = 63; % number of sampling frequencies

% objects to be saved
test_stat   = zeros(4,ndel,3);
rej_prob    = zeros(4,ndel,3);

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
end

% -- write to csv files --
for t  = 1:2
    writematrix(test_stat(:,:,t), savename{t})
    writematrix(rej_prob(:,:,t), savename{t+2})
end

end

