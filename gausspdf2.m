function [q, normConst] = gausspdf2(x, Delta)
% GAUSSPDF  Computes the PDF for a Gaussian
%   Q = GAUSSPDF2(X, DELTA) computes the (multivariate) zero-mean Gaussian
%   distribution with precision DELTA at the points in X. Three possible
%   cases
%    i) X is an (nx x N) matrix of N column vectors, DELTA is an (nx x nx)
%    precision matrix, common to all vectors.
%    ii) X is an arbitrary matrix with scalars, DELTA is a scalar precision
%    of a monovariate Gaussian density.
%    iii) X is an (1 x N) matrix of N scalars, DELTA is an (1 x N) matrix
%    of N precisions of N different monovariate Gaussian densities.
% 
%   [Q, NORMCONST] = GAUSSPDF2(X, DELTA) also returns the normalising
%   constant of the PDF.
%
%   Fredrik Lindsten, 2009-10-19 / 2011-02-15
%   lindsten@isy.liu.se


[nx, sz] = size(Delta); % The rows of delta must equal the dimension
if(sz == nx)
    % Delta - nx x nx, same covariance for all data points
    normConst = sqrt((2*pi)^nx/det(Delta));
    if(nx == 1)
        % x may be a (k-dim) matrix, i.e. M x N x ... scalars
        q = exp(-(1/2)*Delta*x.^2);
    else
        % Must have x - nx x N
        % x'*Rinv*x == sum(x.*(Delta*x), 1)
        q = exp(-(1/2)*sum(x.*(Delta*x), 1));
    end
    q = q/normConst;
elseif(nx == 1)
    % Delta - 1 x N, different covariances for all points - monovariate
    % x     - M x N, matrix of scalars, each colum has the same covariance
    normConst = sqrt(2*pi./Delta);    
    q = exp(-(1/2)* bsxfun(@times, x.^2, Delta));
    q = bsxfun(@rdivide, q, normConst);
%     M = size(x,1);
%     q = exp(-(1/2)* x.^2 .* repmat(Delta, M, 1));
%     q = q./ repmat(normConst, M, 1);    
else
    error('DELTA incorrect size');
end



