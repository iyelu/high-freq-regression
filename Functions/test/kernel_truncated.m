function [wBT] = kernel_uniform(x)

% input 
% x: domain of the kernel density function

% output
% wBT: range of the kernel density function

w = zeros(size(x));
BT = (abs(x) <= 1); 
wBT = w;
wBT(BT) = 1;
