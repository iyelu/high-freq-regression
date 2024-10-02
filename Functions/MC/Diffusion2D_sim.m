function output = Diffusion2D_sim(delta, n, mu, sig11, sig12, sig21, sig22, X0)
% -------------------------------------------------------------------------
% Purpose: Simulation of bivariate diffusion process
%                dXt = mu(Xt) dt + sig(Xt) dW_t
%          where
%                mu(x)  is a 2 by 1 vector valued function 
%                sig(x) is a 2 by 2 matrix valued function
%
% Simulation method: Euler Approximation
% ------------------------------------------
% Input:
%   delta: length of sampling interval
%   n:     number of observations in the simulated sample\
%   mu:    function handle (drift function)
%   sig11, sig12, sig21, sig22: function handle (diffusion functions)
%   X0:    2 by 1 vector of initial values of the component processes
% -------------------------------------------------------------------------

% Total number of observations to simulate initially
nobs = floor(n * 1.1); % simulate 10% more observations for initial burn

% Euler approximation with 30 sub-intervals per sampling interval
delta_euler = delta/30;
n_euler     = nobs * 30; 

% Initialization 
x = zeros(n_euler, 2); x(1,1) = X0(1); x(1,2) = X0(2); 

% Simulation by Euler approximation
sh = normrnd(0,1,2,n_euler);
for i = 2 : n_euler
    xx = x(i - 1, :)';
    sigma = [sig11(xx), sig12(xx); sig21(xx), sig22(xx)];
    y = xx + delta_euler * mu(xx) + sqrt(delta_euler) * sigma * sh(:,i); 
    x(i, :) = y';
end

% Leave 29 observations out of 30, now have the desired frequency
x = x(30:30:end,:);

% Burn the initial observations to make total n observations - 
output = x(end-n+1:end,:);

end