function [kappa, mu, sigma] = SR_mle(V, delta)

% -----------------------------------------------------------------------------
%  Feller's square root (SR) process
%      dVt = kappa * (mu - Vt) + sigma * sqrt(Vt) * dWt    
% -----------------------------------------------------------------------------
% Input: 
%    V: discrete data to fit the SR model 
%    delta: length of sampling interval
% -----------------------------------------------------------------------------
% Output: 
%    kappa, mu, sigma: MLE of parameters      
% -----------------------------------------------------------------------------
% Call functions: 
%    SR_ols.m      (OLS estimation of SR model) 
%    SR_loglik.m   (log-likelihood function of SR model) 
% -----------------------------------------------------------------------------

% Initial values for MLE : OLS estimates 
[kappa0, mu0, sigma0] = SR_ols(V, delta);
theta0  = [kappa0, mu0, sigma0];

% MLE 
options = optimset('LargeScale', 'off', 'MaxIter', 500, 'MaxFunEvals', 500, 'TolFun', 1e-4, 'TolX', 1e-4, 'TolCon', 1e-4);
func    = @(theta) SR_loglik(theta, V, delta);

[theta_mle, Fval, Exitflag] = fminsearch(@(theta) func(theta), theta0, options);

kappa = theta_mle(1); mu = theta_mle(2); sigma = theta_mle(3);

end
