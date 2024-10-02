function x = sim_BM(delta, eta, shock) 

% -- Simulate a discrete sample of Brownian motion ---

  % -------------------------------------
  % delta: length of sampling interval 
  % eta:   volatility of the BM
  % shock: BM shocks
  % -------------------------------------
  
  n = numel(shock); % number of discrete-observations
  x = zeros(n, 1);
  x(2:end) = cumsum(eta * sqrt(delta) * shock(2:end));

end
