function [wPZ, wPZRn] = kernel_parzen(x)

% input 
% x: domain of the kernel density function

% output
% wPZ: range of the kernel density function
% wPZRn: range of the renormalized kernel density function

w = zeros(size(x));
cPZ = 0.539285;
PZ1 = (abs(x) >= 0) & (abs(x) <= 1/2); % range1 of Parzen kernel
PZ2 = (abs(x) >= 1/2) & (abs(x) <= 1); % range2 of Parzen kernel
PZ1Rn = (abs(cPZ*x) >= 0) & (abs(cPZ*x) <= 1/2);
PZ2Rn = (abs(cPZ*x) >= 1/2) & (abs(cPZ*x) <= 1);
wPZ = w;
wPZ(PZ1) = 1-6*x(PZ1).^2+6*abs(x(PZ1)).^3;
wPZ(PZ2) = 2*(1-abs(x(PZ2))).^3;
wPZRn = w;
wPZRn(PZ1Rn) = 1-6*(cPZ*x(PZ1Rn)).^2 ...
    + 6*abs(cPZ*x(PZ1Rn)).^3;
wPZRn(PZ2Rn) = 2*(1-abs(cPZ*x(PZ2Rn))).^3;
