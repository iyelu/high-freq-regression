function [Yt, Xt] = data_sim_main(T, delta, para, mdl, coeff)
% ---------------------------------------------------------------------
% Purpose: data simulation for main stationary and cointegration models
%          in Section 7: Simulation
% ---------------------------------------------------------------------
% Input:
%   T:     length of time span (unit: year)
%   delta: legnth of sampling interval  (unit: year)
%   para:  structure variable containing true parameters in simulation
%          para.stat: parameters in stationary regression model
%          para.coin: parameters in cointegrating regression model
%   mdl:   'stat': stationary regression | 'coin': cointegrating regression
%   coeff: true intercept and slope coefficients in the regression model
% ---------------------------------------------------------------------
% Output:
%   Yt, Xt: one simulated daily discrete sample (sample size: n) of
%           the dependent variable and independent variable
%           (saved in n by 1 vectors)
% ---------------------------------------------------------------------
% Call function:
%   OU_sim.m: simulate OU process
%   BM_sim.m: simulate Brownian motion
% ---------------------------------------------------------------------
n   = floor(T/delta); % number of discrete observations in one sample 
alp = coeff(1);       % true intercept 
bet = coeff(2);       % true slope coefficient

switch mdl
    case 'stat'
        % -- Stationary regression ------------------
        %       Yt = Xt + Ut
        %    where Xt and Ut are both OU processes
        %      dXt = - kapx * Ut dt + sigx dVt
        %      dUt = - kapu * Ut dt + sigu dWt
        %    with V and W two independent standard BM
        % -------------------------------------------
        pX = [0, para.stat.kapx, para.stat.sigx];
        pU = [0, para.stat.kapu, para.stat.sigu];
        sh = randn(n, 2);  % independent shocks for X and U
        Xt = OU_sim(delta, pX, sh(:,1));
        Ut = OU_sim(delta, pU, sh(:,2));
        Yt = alp + bet * Xt + Ut;

    case 'coin'
        % -- Cointegrating regression --------------
        %       Yt = Xt + Ut
        %    where 
        %       Xt is a Brownian Motion with vol sigx
        %       Ut is a stationary OU process:
        %          dUt = - kapu * Ut dt + sigu dWt
        % ------------------------------------------
        pU = [0, para.coin.kapu, para.coin.sigu];
        sh = randn(n, 2);  % independent shocks for X and U
        Ut = OU_sim(delta, pU, sh(:,1));
        Xt = BM_sim(delta, para.coin.sigx, sh(:,2));
        Yt = alp + bet * Xt + Ut;
end

end

