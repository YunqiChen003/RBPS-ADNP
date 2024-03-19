function G=KS_asyn_cor(M,m)
fprintf('KS_asyn_cor\n');
% m: input parameters of dynamic system;
% M: KF output;
tic;
% RTS-smoothing
G.xhats(:,m.ss.T) = M.xhat(:,m.ss.T);
G.Xkserror(:,m.ss.T)=m.x(:,m.ss.T)-G.xhats(:,m.ss.T);
G.Ps(:,:,m.ss.T) = M.P(:,:,m.ss.T);
for n = (m.ss.T-1):-1:1
    G.G1(:,:,n) = M.Prexx(:,:,n)/M.P(:,:,n+1) - M.Prex1y(:,:,n)*M.Kx2y(:,:,n+1)'/M.P(:,:,n+1);
    G.G2(:,:,n) = -M.Prexx(:,:,n)/M.P(:,:,n+1)*M.Kx2y(:,:,n+1) + M.Prex1y(:,:,n)*(M.Kx2y(:,:,n+1)'/M.P(:,:,n+1)*M.Kx2y(:,:,n+1) + inv(M.S(:,:,n+1)));
    G.xhats(:,n) = M.xhat(:,n) + G.G1(:,:,n)*(G.xhats(:,n+1)-M.xpre(:,n+1)) + G.G2(:,:,n)*M.innovation(:,n+1);  % smoothing state
    G.Xkserror(:,n)=m.x(:,n)-G.xhats(:,n);
    G.Pbx_xy(:,:,n) = M.P(:,:,n) - G.G1(:,:,n)*M.Pre(:,:,n+1)*G.G1(:,:,n)' - G.G2(:,:,n)*M.Prex2y(:,:,n+1)'*G.G1(:,:,n)'-...
        G.G1(:,:,n)*M.Prex2y(:,:,n+1)*G.G2(:,:,n)' - G.G2(:,:,n)*M.S(:,:,n+1)*G.G2(:,:,n)'; % backward covariance  var(x_{n}|x_{n+1},Y_{n+1});
    G.Ps(:,:,n) = G.G1(:,:,n)*G.Ps(:,:,n+1)*G.G1(:,:,n)' + G.Pbx_xy(:,:,n)  ; % smoothing covariance
    G.Psxx(:,:,n) = G.G1(:,:,n)*G.Ps(:,:,n+1); % two-step smoothing covariance cov[x_{n},x_{n+1}|Y_{1:T}];
end
G.time=toc;


