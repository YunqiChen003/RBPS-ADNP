function e=BRBPF_LG_2D_ind(m)
fprintf('BRBPF_LG_2D_ind\n');
xparticle=zeros(m.ss.dimxn,m.ss.N,m.ss.T);
xlfilter=zeros(m.ss.dimxl,m.ss.N,m.ss.T);
Plfilter=zeros(m.ss.dimxl,m.ss.dimxl,m.ss.T);
nIdx=zeros(m.ss.N,m.ss.T-1);
W=zeros(m.ss.N,m.ss.T);
Xnhat=zeros(1,m.ss.T);
Xlhat=zeros(1,m.ss.T);
Xhat=zeros(2,m.ss.T);
Xerror=zeros(2,m.ss.T);
NumResampling=0;
tic;
for t=1:m.ss.T
    if t==1
        xparticle(:,:,t)=mvnrnd(m.ss.mn0,m.ss.Pn0,m.ss.N)';
        xlfilter(:,:,t)=repmat(m.ss.ml0,1,m.ss.N);
        Plfilter(:,:,t)=m.ss.Pl0;
        W(:,t)=ones(m.ss.N,1)./m.ss.N;
    else
        %% PF time update
        ESS=ceil(1/sum(W(:,t-1).^2));
        if ESS<=m.ss.Nthres
            nIdx(:,t-1)=resampling(W(:,t-1));% Resampling;
            q=ones(m.ss.N,1)./m.ss.N;
            NumResampling=NumResampling+1;
        else
            nIdx(:,t-1)=(1:m.ss.N)';
            q=W(:,t-1);
        end
        Sigma_nn=m.ss.A(1,2)*Plfilter(:,:,t-1)*m.ss.A(1,2)'+m.ss.Qnn;
        alpha=m.ss.A(1,1)*xparticle(:,nIdx(:,t-1),t-1)+m.ss.A(1,2)*xlfilter(:,nIdx(:,t-1),t-1);
        xparticle(:,:,t)=alpha+sqrt(Sigma_nn)*randn(1,m.ss.N);
        %% KF time update
        % KF dynamical measurement update;
        L=Plfilter(:,:,t-1)*m.ss.A(1,2)'/Sigma_nn;
        Xlfilterstar=xlfilter(:,nIdx(:,t-1),t-1)+L*(xparticle(:,:,t)-alpha);
        Plfilterstar=Plfilter(:,:,t-1)-L*Sigma_nn*L';
        % KF time update;
        flbar=m.ss.Qnl'/m.ss.Qnn*(xparticle(:,:,t)-m.ss.A(1,1)*xparticle(:,nIdx(:,t-1),t-1));
        Xlpre=flbar+m.ss.Albar*Xlfilterstar;
        Plpre=m.ss.Albar*Plfilterstar*m.ss.Albar'+m.ss.Qllbar;
        %% KF measurement update
        M=m.ss.H(1,2)*Plpre*m.ss.H(1,2)'+m.ss.R;
        K=Plpre*m.ss.H(1,2)'/M;
        Plfilter(:,:,t)=Plpre-K*M*K';
        muy=m.ss.H(1,1)*xparticle(:,:,t)+m.ss.H(1,2)*Xlpre;
        inovation=m.y(:,t)-muy;
        xlfilter(:,:,t)=Xlpre+K*(inovation);
        %% PF measurement update
        logly=loggausspdf2(inovation,inv(M))';
        logly=logly-max(logly);
        W(:,t)=q.*exp(logly);
        W(:,t)=W(:,t)./sum(W(:,t));
    end
    %% Compute filtering estimates;
    Xnhat(t)=sum(repmat(W(:,t)',m.ss.dimxn,1).*xparticle(:,:,t),2);
    Xlhat(t)=sum(repmat(W(:,t)',m.ss.dimxl,1).*xlfilter(:,:,t),2);
    Xhat(:,t)=[Xnhat(t);Xlhat(t)];
    Xerror(:,t)=m.x(:,t)-Xhat(:,t);
end
e.time=toc;
e.xparticle=xparticle;
e.xlfilter=xlfilter;
e.Plfilter=Plfilter;
e.nIdx=nIdx;
e.W=W;
e.Xhat=Xhat;
e.Xerror=Xerror;
e.NumResampling=NumResampling;
