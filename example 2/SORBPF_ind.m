function e=SORBPF_ind(m)
fprintf('SORBPF_ind\n');
xparticle=zeros(m.ss.dimxn,m.ss.N,m.ss.T);
xlfilter=zeros(m.ss.dimxl,m.ss.N,m.ss.T);
Plfilter=zeros(m.ss.dimxl,m.ss.dimxl,m.ss.T);
nIdx=zeros(m.ss.N,m.ss.T-1);
W=zeros(m.ss.N,m.ss.T);
Xnhat=zeros(m.ss.dimxn,m.ss.T);
Xlhat=zeros(m.ss.dimxl,m.ss.T);
Xhat=zeros(m.ss.dimx,m.ss.T);
Xerror=zeros(m.ss.dimx,m.ss.T);
NumResampling=0;
tic;
for t=1:m.ss.T
    if t==1
        xparticle(:,:,t)=m.ss.mn0+chol(m.ss.Pn0)'*randn(m.ss.dimxn,m.ss.N);
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
        Sigma_nn=m.ss.An*Plfilter(:,:,t-1)*m.ss.An'+m.ss.Qnn;
        alpha=xparticle(:,nIdx(:,t-1),t-1)+m.ss.An*xlfilter(:,nIdx(:,t-1),t-1);
%         xparticle(:,:,t)=alpha+chol(Sigma_nn)'*randn(m.ss.dimxn,m.ss.N);
%% PF time update;
        % Calculate suboptimal importance density function;
        J=zeros(2,2,m.ss.N);
        J(1,1,:)=alpha(1,:)./sqrt(alpha(1,:).^2+alpha(2,:).^2);
        J(1,2,:)=alpha(2,:)./sqrt(alpha(1,:).^2+alpha(2,:).^2);
        J(2,1,:)=-alpha(2,:)./(alpha(1,:).^2+alpha(2,:).^2);
        J(2,2,:)=alpha(1,:)./(alpha(1,:).^2+alpha(2,:).^2);
        Hso=J;
        myso=mea_eq(alpha);
        G0=multiprod(Hso,Sigma_nn);
        sigmay=bsxfun(@plus,multiprod(G0,multitransp(Hso)),m.ss.R);
        sigmayinv=multi_3d_Matrix_inv(sigmay);
        G=multiprod(multiprod(Sigma_nn,multitransp(Hso)),sigmayinv);
        mnso=alpha+squeeze(multiprod(G,reshape(m.y(:,t)-myso,m.ss.dimy,1,m.ss.N)));
        Ppso=bsxfun(@minus,Sigma_nn,multiprod(G,G0));
        xparticle(:,:,t)=mvnrnd(mnso',Ppso,m.ss.N)';
        %% KF time update
        % KF dynamical measurement update;
        L=Plfilter(:,:,t-1)*m.ss.An'/Sigma_nn;    
        Xlfilterstar=xlfilter(:,nIdx(:,t-1),t-1)+L*(xparticle(:,:,t)-alpha);
        Plfilterstar=Plfilter(:,:,t-1)-L*Sigma_nn*L';
        % KF time update;
        flbar=m.ss.gamanl'*(xparticle(:,:,t)-xparticle(:,nIdx(:,t-1),t-1));
        Xlpre=flbar+m.ss.Albar*Xlfilterstar;
        Plpre=m.ss.Albar*Plfilterstar*m.ss.Albar'+m.ss.Qllbar;
        %% KF measurement update
        M=m.ss.H*Plpre*m.ss.H'+m.ss.R;
        K=Plpre*m.ss.H'/M;
        innovation=m.y(:,t)-mea_eq(xparticle(:,:,t));
        Plfilter(:,:,t)=Plpre-K*M*K';
        xlfilter(:,:,t)=Xlpre+K*innovation;
        %% PF measurement update
        logly=loggausspdf2(xparticle(:,:,t)-alpha,inv(Sigma_nn))'+loggausspdf2(innovation,inv(M))'-log(mvnpdf(xparticle(:,:,t)',mnso',Ppso))+log(q);
        logly=logly-max(logly);
        W(:,t)=exp(logly);
        W(:,t)=W(:,t)./sum(W(:,t));
    end
    %% Compute filtering estimates;
        Xnhat(:,t)=sum(repmat(W(:,t)',m.ss.dimxn,1).*xparticle(:,:,t),2);
        Xlhat(:,t)=sum(repmat(W(:,t)',m.ss.dimxl,1).*xlfilter(:,:,t),2);
        Xhat(:,t)=[Xnhat(:,t);Xlhat(:,t)];
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
