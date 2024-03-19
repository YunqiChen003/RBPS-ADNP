function G=KF_LG_2D_asyn_cor(m)
fprintf('KF_LG_2D_asyn_cor\n');
% m: input parameters of dynamic system;
tic;
% filtering
G.xhat(:,1) = m.ss.m0;
G.P(:,:,1) = m.ss.P0;
G.Xerror(:,1)=m.x(:,1)-G.xhat(:,1);
for n=2:m.ss.T
    G.xpre(:,n) = m.ss.A*G.xhat(:,n-1);          % one step prediction
    G.Prexx(:,:,n-1) = G.P(:,:,n-1)*m.ss.A';    % covariance between x_{n-1} and x_{n}
    G.Pre(:,:,n) = m.ss.A*G.Prexx(:,:,n-1)+m.ss.Q;     % prediction covariance
    G.Prex1y(:,:,n-1)=G.Prexx(:,:,n-1)*m.ss.H';     % covariance between x_{n-1} and y_{n};
    G.Prex2y(:,:,n)=G.Pre(:,:,n)*m.ss.H'+m.ss.M;          % covariance between x_{n} and y_{n};
    
    G.S(:,:,n)=m.ss.H*G.Pre(:,:,n)*m.ss.H' + m.ss.M'*m.ss.H' + m.ss.H*m.ss.M + m.ss.R;
    G.Kx2y(:,:,n) = G.Prex2y(:,:,n)/G.S(:,:,n);
    G.yhat(:,n)=m.ss.H*G.xpre(:,n);
    G.innovation(:,n)=m.y(:,n) - G.yhat(:,n);
    G.xhat(:,n) = G.xpre(:,n) + G.Kx2y(:,:,n)*G.innovation(:,n);   % filtering state
    G.Xerror(:,n)=m.x(:,n)-G.xhat(:,n);
    G.P(:,:,n) = G.Pre(:,:,n) - G.Kx2y(:,:,n)*G.S(:,:,n)*G.Kx2y(:,:,n)';             % filtering covariance   
end
G.time=toc;