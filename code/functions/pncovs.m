function p = pncovs(rho, n, ncovs, tail)

%PNCOVS(rho,n,ncovs) turns back the probability for Pearson's linear
% correlation for testing the hypothesis of no correlation. This function
% has been adapted from MATLAB built-in function 'pvalPearson' (type 
% "open('internal.stats.corrPearson')" for details).
%
%       rho     - vector or matrix of correlation coefficients
%       n       - number of observations
%       ncovs   - number of covariates used to derive residualized
%                 variables
%       tail    - 'both', 'right', or 'left'
%
%This function is supposed to be used on correlations coefficients that have
%been derived from Pearson / Spearman correlations between residualized
%variables. The function allows to specify the number of covariates that have
%been used to derive residuals. Degrees of freedom are adjusted accordingly.
%The resulting p-values correspond to those of partial correlations.
%This might be helpful in scenarios that involve permutation testing where
%the effects of covariates are partialled out before data permutation.

t = rho.*sqrt((n-2-ncovs)./(1-rho.^2)); % +/- Inf where rho == 1

switch tail
    case 'both'
        p = 2*tcdf(-abs(t),n-2-ncovs);
    case 'right'
        p = tcdf(-t,n-2-ncovs);
    case 'left'
        p = tcdf(t,n-2-ncovs);
end

end
