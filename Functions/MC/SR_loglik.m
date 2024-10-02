function lnL = SR_loglik(params, data, delta)
% =========================================================================
% PURPOSE : Log-likelihood objective function (multiplied by -1) for the 
%           Feller's square root (SR) process using MATLAB ncx2pdf function.
% =========================================================================
% Input: 
%           params = model paramsameters (kappa, mu, sigma)
%           data   = discrete sample of SR process
%           delta  = length of sampling interval of the data
% =========================================================================
% Output:
%           lnL = Objective function value 
% =========================================================================

    DataF       = data(2 : end);
    DataL       = data(1 : end - 1);

    kappa       = params(1);
    mu          = params(2);
    sigma       = params(3);

    c  = 2 * kappa / (sigma ^ 2 * (1 - exp(-kappa * delta)));
    q  = 2 * kappa * mu / sigma ^ 2 - 1;
    u  = c * exp(-kappa * delta) * DataL;
    v  = c * DataF;
    nc = 2 * u;
    df = 2 * q + 2; 
    s  = 2 * v;
    
    gpdf = ncx2pdf(s, df, nc);
    ppdf = 2 * c * gpdf;
    lnL  = sum(-log(ppdf));
end
