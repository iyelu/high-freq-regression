function [V, bn] = lrvad91(data, kernel, demean, w)
% ------------------------------------------------------------
% PURPOSE:  Discrete time long-run variance (LRV) estimation
%           using Andrews (1991) AUTOMATIC bandwidth
% ------------------------------------------------------------
% INPUTS:        
%        DATA   - n by k vector of time series data
%        W      - k dimentional weighting vector  
%        KERNEL - Kernel function used in LRV estimation  
%                 'QS' stands for Quadratic Spectral kernel  
%                 'BT' stands for Bartlett kernel  
%                 'PZ' stands for Parzen kernel  
%        W      - weithing vector for elements in multivariate data
%        DEMEAN - Logical true of false (0 or 1) indicating whether 
%                 the mean should be substracted when computing 
%                 the covariance
% ------------------------------------------------------------
% OUTPUTS: 
%         V  - k by k longrun variance estimate
%         BN - Length of bandwidth used in the estimator
% ------------------------------------------------------------
% FUNCTIONS: 
%            kernel_quadratic, kernel_parzen, kernel_bartlett
% ------------------------------------------------------------
% Ye Lu (ye.lu1@sydney.edu.au)
% Last modified: 03/07/2018 

% ------------------------------------------------------------
% Input Checking 
% ------------------------------------------------------------
n = size(data, 1);            
k = size(data, 2);

if nargin == 1
  kernel = 'PZ';
  demean = true;
  w=1;
elseif nargin == 2
  demean = true;
  w=1;
elseif nargin == 3
  w=1;
end

if ndims(data)>2
    error('DATA must be a N by K matrix of data.')
end
if isempty(demean) 
  demain = true;
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
% ------------------------------------------------------------

switch kernel
  case 'QS' % quadratic spectral kernel
    r  = 2;
    cr = 1.3221;
  case 'PZ' % Parzen kernel
    r  = 2;
    cr = 2.6614;
  case 'BT' % bartlett kernel 
    r  = 1;
    cr = 1.1447;
end

% ------------------------------------------------------------
% Determine bandwidth bn:
% - bn = cr * (theta_r^2 * n)^(1/(2r + 1))
% - while theta_r is estimated by assuming data is AR(1)
% ------------------------------------------------------------

% multivariate case of Andrews' bandwidth

top = 0; bot = 0;
for i = 1 : k
  y = data(2: end, i);
  x = data(1: end-1, i);
  rho   = x\y;
  sigsq = data(:,i)'*data(:,i) / n;
  if r == 1
    top = top+w(i)*4*rho^2*sigsq^2/(1-rho)^6/(1+rho)^2;
    bot = bot+w(i)*sigsq^2/(1-rho)^4;
  elseif r == 2
    top = top+w(i)*4*rho^2*sigsq^2/(1-rho)^8;
    bot = bot+w(i)*sigsq^2/(1-rho)^4;
  end
end
thetarsq = top/bot;
bn  = cr*(thetarsq*n)^(1/(2*r+1));

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
  Gammai = (data(i+1:end, :)' * data(1:end-i,:))/n;
  GplusGprime = Gammai + Gammai';
  V = V + kw(i+1) * GplusGprime;
end 

