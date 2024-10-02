function [V, bn] = lrvnw94(data, kernel, demean, w, an)
% ------------------------------------------------------------
% PURPOSE:  Discrete time long-run variance (LRV) estimation
%           using Newey and West (1994) AUTOMATIC bandwidth
% ------------------------------------------------------------
% INPUTS:        
%        DATA   - n by k vector of time series data
%        KERNEL - Kernel function used in LRV estimation  
%                 'QS' stands for Quadratic Spectral kernel  
%                 'BT' stands for Bartlett kernel  
%                 'PZ' stands for Parzen kernel  
%        W      - weighting vector for elements in multivarate data
%        DEMEAN - Logical true of false (0 or 1) indicating whether 
%                 the mean should be substracted when computing 
%                 the covariance
% ------------------------------------------------------------
% OUTPUTS: 
%         V  - Longrun variance estimate
%         BN - Length of bandwidth used in the estimator
% ------------------------------------------------------------
% FUNCTIONS: 
%            kernel_quadratic, kernel_parzen, kernel_bartlett
% ------------------------------------------------------------
% Ye Lu (ye.lu1@sydney.edu.au)
% Last modified: 05/07/2018 

% ------------------------------------------------------------
% Input Checking 
% ------------------------------------------------------------

n = size(data, 1);            

if nargin == 1
  kernel = 'PZ';
  demean = true;
  w = 1;
elseif nargin == 2
  demean = true;
  w = 1;
elseif nargin == 3
  w = 1;
elseif nargin == 4 
  %do nothing 
elseif nargin == 5 
  %do nothing 
else
  error('The number of inputs should be 1-5.')
end

if ndims(data)>2
    error('DATA must be a N by K matrix of data.')
end
if isempty(demean) 
  demean = true;
end 
if ~ismember(demean,[0 1]) 
    error('DEMEAN must be either logical true or false.')
end

% w should be a column vector or a scalar
if size(w,1)<size(w,2)
  w = w';
end 

% If DEMEAN is true, demean the data 
if demean 
  data = data-repmat(mean(data), n, 1);
end 

% ------------------------------------------------------------
% Kernel constants:
% - r:  characteristic exponent of the kernel function  
% - cr: constant in the optimal bandwidth 
%       See equations (6.2) in Andrews (1991)
% - p: useful for the inner bandwidth
% ------------------------------------------------------------

switch kernel
  case 'QS' % quadratic spectral kernel
    r  = 2;
    cr = 1.3221;
    p  = 2/25;
  case 'PZ' % Parzen kernel
    r  = 2;
    cr = 2.6614;
    p  = 4/25;
  case 'BT' % bartlett kernel 
    r  = 1;
    cr = 1.1447;
    p  = 2/9;
end

% Inner bandwidth 
if nargin < 5
  an = round(4*(n/100)^p);
end 

% ------------------------------------------------------------
% Bandwidth:
% - For both AD91 and NW94 bandwidths, the formula is 
% - band = cr * (theta_r^2 * n)^(1/(2r + 1))
% - while theta_r is estimated different in AD and NW
% ------------------------------------------------------------

% an = round(4*(n/100)^p);                      % "inner bandwidth"
top = fq(data, r, an);  % f(q) in Andrews (1991)
bot = fq(data, 0, an);  % f in Andrews (1991)
a = diag(top).^2;
b = diag(bot).^2;
thetarsq = (w'*a)/(w'*b); % a(q) in Andrews (1991) 
bn       = cr*(thetarsq*n)^(1/(2*r+1));         

% ------------------------------------------------------------
% Kernal weights 
% ------------------------------------------------------------
domain = (0:n-1)/bn;                
switch kernel
case 'QS' 
    [kw, ~] = kernel_quadratic(domain); 
  case 'PZ'  
    [kw, ~] = kernel_parzen(domain);    
  case 'BT' 
    [kw, ~] = kernel_bartlett(domain);  
end

% ------------------------------------------------------------
% Compute LRV estimator 
% ------------------------------------------------------------
V = data'*data/n; % variance estimate
V = V * kw(1);
for i = 1:n-1
  Gammai = (data(i+1:end, :)' * data(1:end-i,:))/n;
  GplusGprime = Gammai + Gammai';
  V = V + kw(i+1) * GplusGprime;
end 

return 

% ------------------------------------------------------------
% -- Compute f^(q) in the optimal bandwidth as in Andrews (1991) where q is r here 
%    with truncation parameter an
%
%    f^(q) = sum_j |j|^q * Gamma(j), j = -infty to infty 
%          = Gamma(0) + sum_(j=1)^infty [Gamma(j) + Gamma(j)'] * j^q
% 
%    q = 0, 1, 2 
%    When q = 0, f^(q) is the long-run variance
% ------------------------------------------------------------
function y = fq(data, q, an)
  if size(data,1)<size(data,2) 
    data = data';
  end 

  n = size(data, 1);
  y = data'*data/n;  % Gamma(0)
  for j = 1:(an-1) 
    Gammaj = data(j+1:end, :)' * data(1:end-j, :)/n;
    GplusGprime = Gammaj + Gammaj';
    y = y + j^q*GplusGprime; 
  end 
return 
