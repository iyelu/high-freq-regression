clc;clear;close all
% -------------------------------------------------------------------------
% Purpose: Simulation for the POWER analysis of 
%          stationary and cointegration models presented
%          in Section 7: Simulation of Chang, Lu and Park (2024)
% -------------------------------------------------------------------------
%          Results will be used to generate Table 13 in Appendix D
% -------------------------------------------------------------------------
% Saved results in '../' include 
%          (1) Test statistics (5000 MC, 3 frequencies) under H0: bet = 1
%              'Section_7_Htest_bet_1_ad.csv'
%              'Section_7_Htest_bet_1_crt.csv' 
%              'Section_7_Gtest_bet_1_ad.csv'
%              'Section_7_Gtest_bet_1_crt.csv' 
%          (2) Test statistics (5000 MC, 3 frequencies) under two alternatives
%              'Section_7_Htest_bet_095_ad.csv'   H1: bet = 0.95
%              'Section_7_Htest_bet_095_crt.csv'  H1: bet = 0.95
%              'Section_7_Htest_bet_090_ad.csv'   H1: bet = 0.90
%              'Section_7_Htest_bet_090_crt.csv'  H1: bet = 0.90
%              'Section_7_Gtest_bet_099_ad.csv'   H1: bet = 0.990
%              'Section_7_Gtest_bet_099_crt.csv'  H1: bet = 0.990
%              'Section_7_Gtest_bet_0985_ad.csv'  H1: bet = 0.985
%              'Section_7_Gtest_bet_0985_crt.csv' H1: bet = 0.985
% -------------------------------------------------------------------------
% Call functions:
%   ../../Functions/MC/data_sim_main.m: simulate data according to MC models in Section 7
% -------------------------------------------------------------------------
% Runtime: 8.22 hours on MacBook Air with M2 Chip
% -------------------------------------------------------------------------
% This version: September 2024 by Ye Lu (ye.lu1@sydney.edu.au)
% -------------------------------------------------------------------------

savefolder = '../'; prefix = 'Section_7';
addpath(genpath('../../Functions/'));

niter = 5000; T = 50; % number of MC iterations and time span
del   = [1 21 63];      % daily, monthly, quarterly frequencies

% two regression models with two tests
mdls = {'stat','coin'}; tsts = {'H','G'};

% Parameters in stationary and cointegrating regression models detailed in Section 7
para.stat.kapx=0.1020; para.stat.sigx=1.5513; para.stat.kapu=6.9011; para.stat.sigu=2.7565;
para.coin.sigx=0.0998; para.coin.kapu=1.5718; para.coin.sigu=0.0097;

% Single null hypothesis: bet = 1
R = [0 1]; r = 1;

% three true beta values: one for the null and two for the alternatives
iH = 1; iG = 2;
b{iH} = [1 0.95 0.9];   % true betas for stationary reg. with H-test
b{iG} = [1 0.99 0.985]; % true betas for cointegrating reg. with G-test
bp{iH} = {'1', '095', '090'}; bp{iG} = {'1', '099', '0985'};

tic
for i = 1:2  % two models/tests
    mdl = mdls{i}; tst = tsts{i}; b0 = b{i}; bp0 = bp{i};
    for j = 1:3 % three true beta values (null and two alteranatives)
        truealp = 0; truebet = b0(j); coeff = [truealp truebet];

        % names of the csv files saving the simulated test statistics
        adsave  = sprintf('%s%s_%stest_bet_%s_ad.csv',savefolder,prefix,tst,bp0{j});
        crtsave = sprintf('%s%s_%stest_bet_%s_crt.csv',savefolder,prefix,tst,bp0{j});

        % save test statistics for 3 frequencies
        Wad  = zeros(niter,3); Wcrt = zeros(niter,3);

        for m = 1 : niter
            fprintf('%s, hypo = %d, iter=%d\n', mdl, j, m); rng(m)
            [y, x] = data_sim_main(T, del(1)/252, para, mdl, coeff);  % simulate daily data

            for d = 1:3
                delta = del(d)/252;
                switch mdl
                    case 'coin'
                        [~, Gstat] = Gtest(delta,y(del(d):del(d):end),x(del(d):del(d):end),1,R,r);
                        Wad(m,d) = Gstat.ad91; Wcrt(m,d)  = Gstat.crt;
                    case 'stat'
                        [~, Hstat] = Htest(delta,y(del(d):del(d):end),x(del(d):del(d):end),1,R,r);
                        Wad(m,d) = Hstat.ad91; Wcrt(m,d)  = Hstat.crt;
                end
            end
        end
        % write to csv files
        writematrix(Wad, adsave); writematrix(Wcrt, crtsave)
    end
end
toc

rmpath(genpath('../../Functions/'));
