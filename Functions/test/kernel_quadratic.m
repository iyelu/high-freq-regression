function [wQS, wQSRn] = kernel_quadratic(x)

% input 
% x: domain of the kernel density function

% output
% wQS: range of the kernel density function
% wQSRn: range of the renormalized kernel density function

w = zeros(size(x));
argQS = 6*pi*x/5;
w1 = 3./(argQS.^2);
w2 = (sin(argQS)./argQS)-cos(argQS);
wQS = w1.*w2;
wQS(x == 0) = 1;
wQSRn = wQS; % Renormalization constant = 1