function [b,Gstat,Wald,bw,lrv] = Gtest(delta,y,x,c,R,r)
% --
% Purpose: Robust G test and Nonrobust Wald test
% --
% INPUTS:
%   DELTA   - Sampling interval of the simulated data 
%   Y       - N by 1 vector of dependent data
%   X       - N by K vector of independent data
%   C       - [OPTIONAL] 1 or 0 to indicate whether a constant should be included (1: include constant). The default value is 1.
%   R       - Restriction matrix 
%   r       - H0: R*b=r
% --
% OUTPUTS:
%   B        - A K(+1 is C=1) vector of parameters. If a constant is included, it is the first parameter
%   GSTAT    - K(+1) vector of G-statistics computed using different bandwidths
%             GSTAT.AD91 | GSTAT.NW94 | GSTAT.RT | GSTAT.CRT
%   WALD    - Nonrobust Wald statistic
%   BW      - Structural variable of bandwidth Choices of different schemes 
%             BW.AD91 | BW.NW94 | BW.RT | BW.CRT
%   LRV     - Structural variable of LRV estimates using different bandwidths
%             LRV.AD91 | LRV.NW94 | LRV.RT | LRV.CRT
% ---------------------
% By Ye Lu, Aug 2020

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
        error('The number of columns of X must be grater than or equal to n')
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

% --------------
% Main program  
% --------------

% Compute beta, fitted values, and residuals
b    = x \ y;
yhat = x * b;
uhat = y - yhat;

% -- ad: Andrews 1991 --
[lrv.ad91, bw.ad91] = lrvad91(uhat, 'PZ', 0);
vcv  = inv(x'*x)*lrv.ad91;
Gstat.ad91 = (R*b-r)'*(R*vcv*R')^(-1)*(R*b-r);

% -- nw: Newey West 1994 --
[lrv.nw94, bw.nw94] = lrvnw94(uhat, 'PZ', 0);
vcv  = inv(x'*x)*lrv.nw94;
Gstat.nw94 = (R*b-r)'*(R*vcv*R')^(-1)*(R*b-r);

% -- rt: discrete time rule of thumb --
%    rt = n^1/5
bw.rt = n^(1/5);   
lrv.rt = lrvnw87(uhat,'PZ',0,bw.rt);
vcv  = inv(x'*x)*lrv.rt;
Gstat.rt = (R*b-r)'*(R*vcv*R')^(-1)*(R*b-r);

% -- crt: continuous time rule of thumb --
%    crt = T^1/5 / delta
bw.crt  = min(n^(1/5)/(delta^(4/5)), n);   
lrv.crt = lrvnw87(uhat, 'PZ', 0, 1.5*bw.crt);
vcv  = inv(x'*x) * lrv.crt;
Gstat.crt = (R*b-r)'*(R*vcv*R')^(-1)*(R*b-r);

% -- wd: Wald test --
s2   = uhat'*uhat/numel(y);
vcv  = inv(x'*x)*s2;
Wald = (R*b-r)'*(R*vcv*R')^(-1)*(R*b-r);



