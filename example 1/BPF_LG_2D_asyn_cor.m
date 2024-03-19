function e=BPF_LG_2D_asyn_cor(m)
fprintf('BPF_LG_2D_asyn_cor\n');
xparticle=zeros(m.ss.dimx,m.ss.N,m.ss.T);
nIdx=zeros(m.ss.N,m.ss.T-1);
W=zeros(m.ss.N,m.ss.T);
Xhat=zeros(m.ss.dimx,m.ss.T);
Xerror=zeros(m.ss.dimx,m.ss.T);
NumResampling=0;
%-------------------------------------------------------------------
% Particle filter
%-------------------------------------------------------------------
tic;
for t=1:m.ss.T
    if t==1
        xparticle(:,:,t)=mvnrnd(m.ss.m0,m.ss.P0,m.ss.N)';
        W(:,t)=ones(m.ss.N,1)./m.ss.N;
    else
        ESS=ceil(1/sum(W(:,t-1).^2));
        if ESS<=m.ss.Nthres
            nIdx(:,t-1)=resampling(W(:,t-1));% Resampling;
            q=ones(m.ss.N,1)./m.ss.N;
            NumResampling=NumResampling+1;
        else
            nIdx(:,t-1)=(1:m.ss.N)';
            q=W(:,t-1);
        end
        xtsub1=xparticle(:,nIdx(:,t-1),t-1);
        xparticle(:,:,t)=m.ss.A*xtsub1+mvnrnd(zeros(m.ss.dimx,1),m.ss.Q,m.ss.N)';
        innovation=m.y(:,t)-(m.ss.H+m.ss.C)*xparticle(:,:,t)+m.ss.C*m.ss.A*xtsub1;
        logly=loggausspdf2(innovation,inv(m.ss.R0))';
        logly=logly-max(logly);
        W(:,t)=q.*exp(logly);
        W(:,t)=W(:,t)/sum(W(:,t));
    end
    Xhat(:,t)=sum(repmat(W(:,t)',m.ss.dimx,1).*xparticle(:,:,t),2);
    Xerror(:,t)=m.x(:,t)-Xhat(:,t);
end
e.time=toc;
e.xparticle=xparticle;
e.nIdx=nIdx;
e.W=W;
e.Xhat=Xhat;
e.Xerror=Xerror;
e.NumResampling=NumResampling;
