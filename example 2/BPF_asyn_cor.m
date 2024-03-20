function e=BPF_asyn_cor(m)
fprintf('BPF_asyn_cor\n');
xparticle=zeros(m.ss.dimx,m.ss.N,m.ss.T);
nIdx=zeros(m.ss.N,m.ss.T-1);
W=zeros(m.ss.N,m.ss.T);
Xhat=zeros(m.ss.dimx,m.ss.T);
Xerror=zeros(m.ss.dimx,m.ss.T);
NumResampling=0;
tic;
for t=1:m.ss.T
    if t==1
        %生成初始时刻的样本；
        xparticle(:,:,t)=(mvnrnd(m.ss.m0',m.ss.P0,m.ss.N))';
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
        mx=sys_eq(xtsub1,m);
        xparticle(:,:,t)=mx+sqrtm(m.ss.Qxx)*randn(m.ss.dimx,m.ss.N);
        muy=mea_eq(xparticle(:,:,t))+m.ss.Qxy'/m.ss.Qxx*(xparticle(:,:,t)-mx);
        logly=loggausspdf2(m.y(:,t)-muy,inv(m.ss.R0))';
        logly=logly-max(logly);
        W(:,t)=q.*exp(logly); 
        W(:,t)=W(:,t)./sum(W(:,t));
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
