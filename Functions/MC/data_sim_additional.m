function [Yt, Xt] = data_sim_additional(T, delta, para, mdl, coeff)
% -------------------------------------------------------------------------
% Purpose: data simulation for stationary and cointegration regressions
%          in Appendix E: additional simulation results
% -------------------------------------------------------------------------
% Input:
%   T:     length of time span (unit: year)
%   delta: legnth of sampling interval  (unit: year)
%   para:  structure variable containing true parameters in simulation
%          para.stat: parameters in stationary regression model
%          para.coin: parameters in cointegrating regression model
%   mdl:   'stat': stationary regression | 'coin': cointegrating regression
%   coeff: true intercept and slope coefficients in the regression model
% -------------------------------------------------------------------------
% Output:
%   Yt, Xt: one simulated daily discrete sample (sample size: n) of
%           the dependent variable and independent variable
%           (saved in n by 1 vectors)
% -------------------------------------------------------------------------
% Call function:
%   ../Functions/MC/SR_sim.m:           simulate Feller's SR process
%   ../Functions/MC/OU_sim.m:           simulate Ornstein-Uhlenbeck process
%   ../Functions/MC/Duffusion2D_sim. m: simulate bivariate diffusion 
% -------------------------------------------------------------------------
n   = floor(T/delta); % number of discrete observations in one sample 
alp = coeff(1);       % true intercept 
bet = coeff(2);       % true slope coefficient

switch mdl
    case 'stat'
        % -- Stationary regression ------------------------------------
        %       Yt = Xt + Ut
        %    where Xt is a Feller's square root (SR) process:
        %            dXt = kap_x * (mu_x - Xt) * dt + sig_x dVt
        %    and Ut is a OU process:  
        %            dUt = -kap_u * Ut dt + sig_u dWt
        %    where V and W are standard BMs with correlation rho
        % -------------------------------------------------------------
        pX  = [para.stat.mux, para.stat.kapx, para.stat.sigx];
        pU  = [0,             para.stat.kapu, para.stat.sigu];
        rho = para.stat.rho;  

        % -- simulate correlated shocks of X and U --------------------
        nosim = 1;
        sh = randn(2 * nosim, n);
        co = [1, rho; rho, 1]; % corr of shocks of dVt and dWt
        a  = chol(co);
        b  = kron(eye(nosim), a');
        sh = b * sh;           % normalized shocks of dWt and dVt (2*nosim by n)
        sh_x = sh(1:2:end, :)';
        sh_u = sh(2:2:end, :)';

        % -- simulate processes X, U and Y ----------------------------
        Xt = SR_sim(delta, pX, sh_x);
        Ut = OU_sim(delta, pU, sh_u);
        Yt = alp + bet * Xt + Ut;

    case 'coin'
        % -- Cointegration type of regression -------------------------
        %       Yt = Xt + Ut
        %    where Xt follows the Heston (1993) SV model
        %    and Ut has the GARCH stochastic volatility
        %    both Xt and Ut are bivariate diffusion processes
        % -------------------------------------------------------------       

        output = Diffusion2D_sim(delta, n, para.coin.Xmu, para.coin.Xsig11, ...
            para.coin.Xsig12,para.coin.Xsig21,para.coin.Xsig22,para.coin.Xinitial);
        Xt = output(:,1);

        output = Diffusion2D_sim(delta, n, para.coin.Umu, para.coin.Usig11, ...
            para.coin.Usig12,para.coin.Usig21,para.coin.Usig22,para.coin.Uinitial);
        Ut = output(:,1);
        
        Yt = alp + bet * Xt + Ut;
end
end




