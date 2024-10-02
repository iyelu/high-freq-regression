function V = lrvnw87(data, kernel, demean, bn)
% ------------------------------------------------------------
% PURPOSE:  Discrete time long-run variance (LRV) estimation
%           using Newey-West NONAUTOMATIC bandwidth
% ------------------------------------------------------------
% INPUTS:        
%        DATA   - n by k vector of time series data
%        KERNEL - Kernel function used in LRV estimation  
%                 'QS' stands for Quadratic Spectral kernel  
%                 'BT' stands for Bartlett kernel  
%                 'PZ' stands for Parzen kernel  
%        DEMEAN - Logical true of false (0 or 1) indicating whether 
%                 the mean should be substracted when computing 
%                 the covariance
%        BN     - Non-negative integer specifying the bandwidth 
%                 (or the number of lags). If empty or not included, 
%                 BN = min(floor(4*(n/100)^(2/9)), n) is used.
% ------------------------------------------------------------
% OUTPUTS: 
%         V  - k by k longrun variance matrix
% ------------------------------------------------------------
% FUNCTIONS: 
%            kernel_quadratic, kernel_parzen, kernel_bartlett
% ------------------------------------------------------------
% Ye Lu (ye.lu1@sydney.edu.au)
% Last modified: 29/06/2018 

% ------------------------------------------------------------
% Input Checking 
% ------------------------------------------------------------
n = size(data, 1);            

if nargin == 1
  kernel = 'BT';
  demean = true;
  bn     = min(floor(4*(n/100)^(2/9)), n);
end 

switch kernel
  case 'QS' 
    p  = 2/25;
  case 'PZ'  
    p  = 4/25;
  case 'BT'  
    p  = 2/9;
end

if nargin == 2
  demean = true;
  bn     = min(floor(4*(n/100)^p), n);
elseif nargin == 3
  bn     = min(floor(4*(n/100)^p), n);
end

if ndims(data)>2
    error('DATA must be a T by K matrix of data.')
end
if isempty(demean) 
  demain = true;
end 
if ~ismember(demean,[0 1]) 
    error('DEMEAN must be either logical true or false.')
end

% If DEMEAN is true, demean the data 
if demean 
  data = data-repmat(mean(data), n, 1);
end 

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
V = data' * data / n; % variance estimate
V = V * kw(1);
for i = 1:n-1
  Gammai = (data(i+1:end, :)' * data(1:end-i, :))/n;
  GplusGprime = Gammai + Gammai';
  V = V + kw(i+1) * GplusGprime;
end 



