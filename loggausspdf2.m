function q= loggausspdf2(x, Delta)
% GAUSSPDF  Computes the logarithm of a Gaussian PDF
%   Q = LOGGAUSSPDF2(X, DELTA) computes the logarithm of the (multivariate)
%   zero-mean Gaussian distribution with precision DELTA at the points in
%   X. Three possible cases
%
%    i) X is an (nx x N) matrix of N column vectors, DELTA is an (nx x nx)
%    precision matrix, common to all vectors.
%
%    ii) X is an arbitrary matrix with scalars, DELTA is a scalar precision
%    of a monovariate Gaussian density.
%
%    iii) X is an (M x N) matrix of scalars, DELTA is an (1 x N) matrix
%    of N precisions of N different monovariate Gaussian densities.
% 
%   [Q, NORMCONST] = LOGGAUSSPDF2(X, DELTA) also returns the normalising
%   constant of the PDF.
%
%   Fredrik Lindsten, 2009-10-19 / 2011-02-15
%   lindsten@isy.liu.se


[nx, sz] = size(Delta); % The rows of delta must equal the dimension
if(sz == nx)
    % Delta - nx x nx, same covariance for all data points
    lZ = nx/2*log(2*pi) - 1/2*log(det(Delta));
    if(nx == 1)
        % x may be a (k-dim) matrix, i.e. M x N x ... scalars
        q = -(1/2)*Delta*x.^2;
    else
        % Must have x - nx x N
        % x'*Rinv*x == sum(x.*(Delta*x), 1)
        q = -(1/2)*sum(x.*(Delta*x), 1);
    end
    q = q - lZ;
elseif(nx == 1)
    % Delta - 1 x N, different covariances for all points - monovariate
    % x     - M x N, matrix of scalars, each colum has the same covariance
    lZ = 1/2*log(2*pi) - 1/2*log(Delta);    
    q = -(1/2)* bsxfun(@times, x.^2, Delta);
    q = bsxfun(@minus, q, lZ);    
    % M = size(x,1);
    % q = -(1/2)* x.^2 .* repmat(Delta, M, 1);
    % q = q - repmat(log(normConst), M, 1);
else
    error('DELTA incorrect size');
end



