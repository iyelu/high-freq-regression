function [b, s2, vcv, tstat, R2, Rbar, uhat] = myols(y,x,c,b0,white,nw87,andrews)
% Linear regression estimation with robust standard errors.
%
% USAGE:
%   [B,TSTAT,S2,VCVNW,R2,RBAR,YHAT] = olsnw(Y,X,C,NWLAGS)
%
% INPUTS:
%   Y       - T by 1 vector of dependent data
%   X       - T by K vector of independent data
%   C       - [OPTIONAL] 1 or 0 to indicate whether a constant should be included (1: include constant). The default value is 1.
%   B0      - true parameters
%   WHITE   - 1 or 0 logical input to include/exclude White robust cov estimate (only take into account heteroskedasticity but not serial correlation)
%   NW87    - 1 or 0 logical input to include/exclude Newey-West (1987) robust cov estimate (nonautomatic bandwidth)
%   ANDREWS - 1 or 0 logical input to include/exclude Andrews (1991) robust cov estimate (automatic HAC estimate)
%
% OUTPUTS:
%   B        - A K(+1 is C=1) vector of parameters. If a constant is included, it is the first parameter
%   S2   - Estimated error variance of the regression 
%   VCV  - Variance-covariance matrix of the estimated parameters
%     VCV.NONROBUST 
%     VCV.WHITE: White robust (only heteroskedasticity, no serial correlation)
%     VCV.NW87:  Newey-West robust (nonautomatic bandwith) 
%     VCV.AD91:  Andrews (1991) robust (automatic bandwidth
%   TSTAT - K(+1) vector of t-statistics computed using different standard errors
%     TSTAT.NONROBUST 
%     TSTAT.WHITE 
%     TSTST.NW87
%     TSTAT.AD91
%   R2   - R-squared of the regression.  Centered if C=1
%   RBAR - Adjusted R-squared. Centered if C=1
%   UHAT - Residuals
%
% COMMENTS:
%   The model estimated is Y = X * B + u where Var(u)=sigma^2.
%
% EXAMPLES:
%   Regression with automatic BW selection
%       b = ols(y,x)
%   Regression without a constant
%       b = ols(y,x,0)
%   Regression with White standard errors
%       b = ols(y,x,1,0)
% ---------------------
% Copyright: Ye Lu
% Date: 03/07/2018


% ---------------------
% Input Checking
% ---------------------

% Check y
if size(y,1)<size(y,2)
    y=y';
end
if size(y,2)~=1
    error('Y must be a column vector')
end
T = size(y,1);

% Check X
if ~isempty(x)
    if size(x,1)~=T
        error('X must have the same number of rows as Y.')
    end
    if size(x,2)>T
        error('The number of columns of X must be grater than or equal to T')
    end
    if rank(x)<size(x,2)
        error('X is rank deficient')
    end
end

K=size(x,2);

% Check c
if isempty(c)
    c=1;
end

if ~ismember(c,[0 1])
    error('C must be either 0 or 1.')
end

% Check c,X
if c==0 && size(x,2)==0
    error('The model must include a constant or at least one X')
end
if c==1
    x=[ones(T,1) x];
    K=K+1;
end 

% Check number of argument in the function 
switch nargin
    case 2
        c=1; b0=zeros(K,1); 
    case 3
        b0=zeros(K,1); white=0; nw87=0; andrews=0;
    case 4
        white=0; nw87=0; andrews=0;
    case 5
        nw87=0; andrews=0;
    case 6
        andrews=0;
    case 7
        % do nothing
    otherwise
        error('2 to 7 inputs only')
end

% Check b0
if size(b0,1)<size(b0,2)
    b0=b0';
end

% ---------------------
% Main program for OLS  
% ---------------------

% Compute beta, fitted values, and residuals

b    = x \ y;
yhat = x * b;
uhat = y - yhat;

% Estimated residual variance and (nonrobust) estimated parameter covariance

s2  = uhat' * uhat / T;
vcv.nonrobust = inv(x' * x) * s2;
tstat.nonrobust = (b-b0) ./ sqrt(diag(vcv.nonrobust));

% R2 and Rbar (adjusted R2)

if c==1 % Use centered if a constant is included
    ytilde=y-mean(y);
    R2=1 - (uhat'*uhat)/(ytilde'*ytilde);
    Rbar=1 - (uhat'*uhat)/(ytilde'*ytilde) * (T-1)/(T-K);
else % Use non centered versions
    R2=1 - (uhat'*uhat)/(y'*y);
    Rbar=1 - (uhat'*uhat)/(y'*y) * (T-1)/(T-K);
end

% ---------------------
% Robust covariance of estimated parameters (betahat) 
% ---------------------

% 1. White robust covariance  

if white
 xu  = x .* repmat(uhat, 1, K);
 XpX = inv(x' * x);
 vcv.white   = XpX * (xu' * xu) * XpX;
 tstat.white = (b-b0) ./ sqrt(diag(vcv.white));
end

% 2. Newey-West (1987) robust covariance (nonautomatric bandwith) 

if nw87
 xu  = x .* repmat(uhat, 1, K);
 lrv = lrvnw87(xu, 'PZ', 0);
 vcv.nw87   = inv(x'*x/T) * lrv * inv(x'*x/T) / T;
 tstat.nw87 = (b-b0) ./ sqrt(diag(vcv.nw));
end

% 3. Andrews (1991) robust covariance (automatric bandwith) 

if andrews
 xu = x .* repmat(uhat, 1, K);
 w = zeros(K,1); 
 w(end) = 1;
 lrv = lrvad91(xu, 'PZ', 0, w);
 vcv.ad91   = inv(x'*x/T) * lrv * inv(x'*x/T) / T;
 tstat.ad91 = (b-b0) ./ sqrt(diag(vcv.ad91));
end 



