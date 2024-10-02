function [kappahat, muhat, sigmahat] = SR_ols(V, delta)

% -----------------------------------------------------------------------------
%  Feller's square root (SR) process
%      dVt = kappa * (mu - Vt) + sigma * sqrt(Vt) * dWt    
% -----------------------------------------------------------------------------
% Input: 
%    V: discrete data to fit the SR model 
%    delta: length of sampling interval
% -----------------------------------------------------------------------------
% Output: 
%    kappa, mu, sigma: OLS estimates of parameters

%    OLS is based on the discretized regression model: 

%    [V(i) - V(i - 1)] / sqrt[delta * V(i - 1)] = 
%       [mu * sqrt(delta) / sqrt(V(i - 1)) - sqrt(delta * V(i - 1)] * kappa 
%       + sigma * epsilon(i),   epsilon(i) ~ N(0, 1) 
% -----------------------------------------------------------------------------

if isrow(V)
  V = V';
end

Vlag = V(1 : end - 1);

muhat = mean(V); % OLS estimate for mu 

% regressand (y) and regressors (x) of the discretized model
y     = diff(V) ./ sqrt(delta * Vlag);   
x     = muhat * sqrt(delta) ./ sqrt(Vlag) - sqrt(delta * Vlag);

kappahat = x \ y; % OLS estimate in the discretized regression model
resid    = y - kappahat * x;
sigmahat = sqrt(var(resid, 1));

end
