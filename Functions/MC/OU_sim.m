function [x] = OU_sim(delta, param, shock)

% -- Simulate a discrete sample of Ornstein-Uhlenbeck process ---
  
mu    = param(1);
kappa = param(2);
sigma = param(3);

n = numel(shock);
x = zeros(n, 1);
% Initial distribution: Normal(mu, sigma^2/ (2 * kappa))
x(1) = mu + sqrt(sigma^2 / (2 * kappa)) * shock(1);
for i = 2 : n
    x(i) = (mu + exp(- kappa * delta) * (x(i-1, :) - mu)) + ...
        sigma * sqrt((1 - exp(- 2 * kappa * delta))/ (2 * kappa))  * shock(i);
end
end
