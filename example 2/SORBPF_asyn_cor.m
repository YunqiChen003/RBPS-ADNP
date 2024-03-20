function e=SORBPF_asyn_cor(m)
fprintf('SORBPF_asyn_cor\n');
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
        xparticle(:,:,t)=m.ss.mn0+sqrtm(m.ss.Pn0)*randn(m.ss.dimxn,m.ss.N);
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
        %% KF time update
        % KF dynamical measurement update;
        L=Plfilter(:,:,t-1)*m.ss.An'/Sigma_nn;
        Plfilterstar=Plfilter(:,:,t-1)-L*Sigma_nn*L';
        % KF time update;
        Plpre=m.ss.Albar*Plfilterstar*m.ss.Albar'+m.ss.Qllbar;
        %% KF measurement update
        M=m.ss.Hbar1*Plfilterstar*m.ss.Hbar1'+m.ss.Hbar2*m.ss.Albar*Plfilterstar*m.ss.Hbar1'+...
            m.ss.Hbar1*Plfilterstar*m.ss.Albar'*m.ss.Hbar2'+m.ss.Hbar2*Plpre*m.ss.Hbar2'+m.ss.R0;%PF measurement update;
        K=(m.ss.Albar*Plfilterstar*m.ss.Hbar1'+Plpre*m.ss.Hbar2')/M;%KF measurement update;
        Plfilter(:,:,t)=Plpre-K*M*K';
        %% PF time update;
        % Calculate suboptimal importance density function;
        J=zeros(2,2,m.ss.N);
        J(1,1,:)=alpha(1,:)./sqrt(alpha(1,:).^2+alpha(2,:).^2);
        J(1,2,:)=alpha(2,:)./sqrt(alpha(1,:).^2+alpha(2,:).^2);
        J(2,1,:)=-alpha(2,:)./(alpha(1,:).^2+alpha(2,:).^2);
        J(2,2,:)=alpha(1,:)./(alpha(1,:).^2+alpha(2,:).^2);
        Hso=bsxfun(@plus,m.ss.D1*L+m.ss.D3,J);
        bso=mea_eq(alpha)-squeeze(multiprod(J,reshape(alpha,m.ss.dimxn,1,m.ss.N)))+m.ss.D1*(xlfilter(:,nIdx(:,t-1),t-1)-L*alpha)-m.ss.D3*xparticle(:,nIdx(:,t-1),t-1);
        myso=squeeze(multiprod(Hso,reshape(alpha,m.ss.dimxn,1,m.ss.N)))+bso;
        G0=multiprod(Hso,Sigma_nn);
        sigmay=bsxfun(@plus,multiprod(G0,multitransp(Hso)),M);
        sigmayinv=multi_3d_Matrix_inv(sigmay);
        G=multiprod(multiprod(Sigma_nn,multitransp(Hso)),sigmayinv);
        mnso=alpha+squeeze(multiprod(G,reshape(m.y(:,t)-myso,m.ss.dimy,1,m.ss.N)));
        Ppso=bsxfun(@minus,Sigma_nn,multiprod(G,G0));
        xparticle(:,:,t)=mvnrnd(mnso',Ppso,m.ss.N)';
        %% KF time update
        % KF dynamical measurement update;
        Xlfilterstar=xlfilter(:,nIdx(:,t-1),t-1)+L*(xparticle(:,:,t)-alpha);
        flbar=m.ss.gamanl'*(xparticle(:,:,t)-xparticle(:,nIdx(:,t-1),t-1));
        % KF time update;
        Xlpre=flbar+m.ss.Albar*Xlfilterstar;
        %% KF measurement update
        muy=mea_eq(xparticle(:,:,t))+m.ss.gamany'*(xparticle(:,:,t)-xparticle(:,nIdx(:,t-1),t-1))-m.ss.gamalybar'*flbar+...
            m.ss.Hbar1*Xlfilterstar+m.ss.Hbar2*Xlpre;
        inovation=m.y(:,t)-muy;
        xlfilter(:,:,t)=Xlpre+K*(inovation);
%         q1=gausspdf2(xparticle(:,:,t)-alpha,inv(Sigma_nn))';
%         q2=gausspdf2(inovation,inv(M))';
%         q3=mvnpdf(xparticle(:,:,t)',mnso',Ppso);
        %% PF measurement update
%         W(:,t)=q.*q1.*q2./q3;
        logly=loggausspdf2(xparticle(:,:,t)-alpha,inv(Sigma_nn))'+loggausspdf2(inovation,inv(M))'-log(mvnpdf(xparticle(:,:,t)',mnso',Ppso))+log(q);
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

