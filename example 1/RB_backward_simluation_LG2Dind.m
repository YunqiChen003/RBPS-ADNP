function G=RB_backward_simluation_LG2Dind(e,m)
% e: output of Rao-Blackwellzed particle filters;
% m: model parametes;
fprintf('RB_backward_simluation_LG2Dind\n');
tic;
% Backward simulator
ind=sample_ancestors(e.W(:,m.ss.T),m.ss.Ns);
G.x_bwd(:,:,m.ss.T) = e.xparticle(:,ind,m.ss.T);
G.xnhat(:,m.ss.T)=mean(G.x_bwd(:,:,m.ss.T),2);
for t = (m.ss.T-1) : (-1) : 1
    %% Backward information filter update;
    if t==m.ss.T-1
        Omega_hat=m.ss.H(1,2)'/m.ss.R*m.ss.H(1,2);
    else
        Omega_hat=G.Omega{t+1}+m.ss.H(1,2)'/m.ss.R*m.ss.H(1,2);
    end
    %% Backward information filter prediction;
    halfQ=chol(m.ss.Qllbar)';
    M=halfQ'*Omega_hat*halfQ+eye(m.ss.dimxl);
    Psi=m.ss.Albar'*Omega_hat*m.ss.Albar-m.ss.Albar'*Omega_hat*halfQ/M*halfQ'*Omega_hat*m.ss.Albar;
    Omega=Psi+m.ss.An'/m.ss.Qnn*m.ss.An;
    %% Compute predictive PDF;
    halfPfilter=chol(e.Plfilter(:,:,t))';
    Upsilon=halfPfilter'*Omega*halfPfilter+eye(m.ss.dimxl);
    for j = 1:m.ss.Ns
        %% Backward information filter update;
        xn1minusfn=G.x_bwd(:,j,t+1)-e.xparticle(:,:,t);
        flbar=m.ss.gamanl'*xn1minusfn;
        innovation=m.y(:,t+1)-m.ss.H(1,1)*G.x_bwd(:,j,t+1);
        if t==m.ss.T-1
            lambda_hat_cand=m.ss.H(1,2)'/m.ss.R*innovation;
        else
            lambda_hat_cand=G.lambda{t+1}(:,j)+m.ss.H(1,2)'/m.ss.R*innovation;
        end
        %% Backward information filter prediction;
        m1=lambda_hat_cand-Omega_hat*flbar;
        halfQm=halfQ'*m1;
        tau=sum(flbar.*(Omega_hat*flbar),1)-2*sum(lambda_hat_cand.*flbar,1)-sum(halfQm.*(M\halfQm),1);
        eta=(m.ss.Albar'-m.ss.Albar'*Omega_hat*halfQ/M*halfQ')*m1;
        logZ=-(1/2)*(tau+sum(xn1minusfn.*(m.ss.Qnn\xn1minusfn),1));
        lambda_cand=eta+m.ss.An'/m.ss.Qnn*xn1minusfn;
        %% Compute predictive PDF;
        halfPlOx=halfPfilter'*(lambda_cand-Omega*e.xlfilter(:,:,t));
        kappa=sum(e.xlfilter(:,:,t).*(Omega*e.xlfilter(:,:,t)),1)-2*sum(lambda_cand.*e.xlfilter(:,:,t),1)-sum(halfPlOx.*(Upsilon\halfPlOx),1);
        %% Compute backward sampling weights;
        logws=logZ-1/2*kappa;
        logws=logws-max(logws);
        ws=e.W(:,t).*(exp(logws))';
        ws=ws./sum(ws);
        %% Backward sampling;
        index=sample_ancestors(ws,1);
        G.x_bwd(:,j,t) = e.xparticle(:,index,t);
        G.lambda{t}(:,j)=lambda_cand(:,index);
    end
    G.Omega{t}=Omega;
    G.xnhat(:,t)=mean(G.x_bwd(:,:,t),2);
end
G.time=toc;

