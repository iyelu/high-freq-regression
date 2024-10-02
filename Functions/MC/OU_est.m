function [sighat, kaphat, muhat] = OU_est(data, delta)

% ---------------------------------------------------------------------
% Purpose: Fit the Ornstein-Uhlenbeck model to a discrete sample
% ---------------------------------------------------------------------
% Input:
%   data:  a discrete sample of OU process save in a vector
%   delta: length of sampling interval (unit: year) of the data
%          e.g. delta = 1/252 for daily data and 1/4 for quarterly data
% ---------------------------------------------------------------------
% Output: 
%   sighat: estimated volatility parameter
%   kaphat: estimated mean-reversion parameter
%   muhat:  estimated mean parameter
% ---------------------------------------------------------------------

n = numel(data); % discrete sample size
T = n * delta;   % length of time span in years
d = diff(data);  
qv = d'*d;       % estimated total quadratic variation 
sighatsq = qv/T;

sighat = sqrt(sighatsq);
kaphat = 0.5 * sighatsq/var(data,1);
muhat  = mean(data);
