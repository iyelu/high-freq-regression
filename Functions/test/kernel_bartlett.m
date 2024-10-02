function [wBT, wBTRn] = kernel_bartlett(x)

% input 
% x: domain of the kernel density function

% output
% wBT: range of the kernel density function
% wBTRn: range of the renormalized kernel density function

w = zeros(size(x));
cBT = 2/3; 
BT = (abs(x) <= 1); 
BTRn = (abs(cBT * x) <= 1); 
wBT = w;
wBT(BT) = 1 - abs(x(BT));
wBTRn = w;
wBTRn(BTRn) = 1 - abs(cBT * x(BTRn));
