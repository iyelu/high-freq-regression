function output = SR_sim(delta, param, sh)
% ------------------------------------------
% Purpose: Simulation of a Feller's square root diffusion model
%                dXt = mu(Xt) dt + sig(Xt) dW_t
%          where
%                mu(x)=kap_x*(mu_x-x), sig(x)=sig_x*sqrt(x)
%          using Milstein Approximation
% ------------------------------------------
% Input:
%   T, delta: time span and sampling interval
%   param:    model parameters
%   sh:       Brownian motion shocks
% ------------------------------------------

% number of observations in the simulated process
n = size(sh,1);

% parameters
mu_x  = param(1);
kap_x = param(2);
sig_x = param(3);

% Define the functions used in Milstein Approximation
mu = @(x)(kap_x * (mu_x - x));
sigma = @(x)(sig_x * sqrt(x));
sigmaprime = @(x)(sig_x/(2 * sqrt(x)));

x = zeros(n, 1);

% Initial draw from time invariant distribution Gamma(alpha,beta)
alpha = 2 * kap_x * mu_x / sig_x^2;
beta = 2 * kap_x / sig_x^2;
x(1) = gamrnd(alpha, 1/beta);

for i = 2 : n
    xx = x(i - 1);
    x(i) = xx + delta * mu(xx) ...
              + sqrt(delta) * sigma(xx) * sh(i) ...
              + delta * (sigma(xx) * sigmaprime(xx) / 2) * (sh(i)^2 - 1);
end

output = x;

end
