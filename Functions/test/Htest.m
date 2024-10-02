function [b,Hstat,bw,lrv] = Htest(delta,y,x,c,R,r)
% --
% Purpose: Robust H test 
% --
% INPUTS:
%   DELTA   - sampling interval of the simulated data
%   Y       - n by 1 vector of dependent data
%   X       - n by K vector of independent data
%   C       - [OPTIONAL] 1 or 0 to indicate whether a constant should be included (1: include constant). The default value is 1.
%   R       - Restriction matrix 
%   r       - H0: R*b=r
% --
% OUTPUTS:
%   B       - A K(+1 is C=1) vector of parameters. If a constant is included, it is the first parameter
%   HSTAT   - K(+1) vector of H statistic computed using different bandwidths
%             HSTAT.AD91 | HSTAT.NW94 | HSTAT.RT | HSTAT.CRT
%   BW      - Structural variable of bandwidth Choices of different schemes 
%             BW.AD91 | BW.NW94 | BW.RT | BW.CRT
%   LRV     - Structural variable of LRV estimates using different bandwidths
%             LRV.AD91 | LRV.NW94 | LRV.RT | LRV.CRT
% ---------------------
% Ye Lu | This version: Aug 2020

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
n = size(y,1);

% Check X
if ~isempty(x)
  if size(x,1)~=n
      error('X must have the same number of rows as Y.')
  end
  if size(x,2)>n
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
    x=[ones(n,1) x];
    K=K+1;
end 

% Check r 
if size(r,1)<size(r,2)
    r=r';
end

% ---------------------
% Main program for OLS  
% ---------------------

% Compute beta, fitted values, and residuals
b    = x \ y;
yhat = x * b;
uhat = y - yhat;

w = [0.5 0.5];
xu  = x .* repmat(uhat, 1, K);

% -- Andrews 1991 --
[lrv.ad91, bw.ad91] = lrvad91(xu, 'PZ', 0, w);
vcv  = (x'*x)^(-1)*(n*lrv.ad91)*(x'*x)^(-1);
Hstat.ad91 = (R*b-r)'*(R*vcv*R')^(-1)*(R*b-r);

% -- Newey West 1994 --
[lrv.nw94, bw.nw94]  = lrvnw94(xu, 'PZ', 0, w);
vcv  = (x'*x)^(-1)*(n*lrv.nw94)*(x'*x)^(-1);
Hstat.nw94 = (R*b-r)'*(R*vcv*R')^(-1)*(R*b-r);

% -- Discrete time rule of thumb --
bw.rt = n^(1/5);   % "wrong" rule of thumb
lrv.rt = lrvnw87(xu, 'PZ', 0, bw.rt);
vcv  = (x'*x)^(-1)*(n*lrv.rt)*(x'*x)^(-1);
Hstat.rt = (R*b-r)'*(R*vcv*R')^(-1)*(R*b-r);

% -- Continuous time rule of thumb --
bw.crt = min(n^(1/5)/(delta^(4/5)),n);   % "right" rule of thumb
lrv.crt = lrvnw87(xu, 'PZ', 0, 1.5*bw.crt);
vcv  = (x'*x)^(-1)*(n*lrv.crt)*(x'*x)^(-1);
Hstat.crt = (R*b-r)'*(R*vcv*R')^(-1)*(R*b-r);

