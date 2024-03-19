function G=RB_backward_simluation_LG2Dasyn_cor(e,m)
% e: output of Rao-Blackwellzed particle filters;
% m: model parametes;
fprintf('RB_backward_simluation_asyn_cor\n');
tic;
% Backward simulator
ind=sample_ancestors(e.W(:,m.ss.T),m.ss.Ns);
G.x_bwd(:,:,m.ss.T) = e.xparticle(:,ind,m.ss.T);
G.xnhat(:,m.ss.T)=mean(G.x_bwd(:,:,m.ss.T),2);

for t = (m.ss.T-1) : (-1) : 1
    %% Backward information filter update;
    Omega_hat_11=m.ss.Hbar1'/m.ss.R0*m.ss.Hbar1;
    Omega_hat_12=m.ss.Hbar1'/m.ss.R0*m.ss.Hbar2;
    if t==m.ss.T-1
        Omega_hat_22=m.ss.Hbar2'/m.ss.R0*m.ss.Hbar2;
    else
        Omega_hat_22=G.Omega{t+1}+m.ss.Hbar2'/m.ss.R0*m.ss.Hbar2;
    end
    %% Backward information filter prediction;
    halfQ=chol(m.ss.Qllbar)';
    M=halfQ'*Omega_hat_22*halfQ+eye(m.ss.dimxl);
    Phi=Omega_hat_12+m.ss.Albar'*Omega_hat_22;
    Psi=Omega_hat_11+Omega_hat_12*m.ss.Albar+m.ss.Albar'*Omega_hat_12'+m.ss.Albar'*Omega_hat_22*m.ss.Albar-Phi*halfQ/M*halfQ'*Phi';
    Omega=Psi+m.ss.An'/m.ss.Qnn*m.ss.An;
    %% Compute predictive PDF;
    halfPfilter=chol(e.Plfilter(:,:,t))';
    Upsilon=halfPfilter'*Omega*halfPfilter+eye(m.ss.dimxl);
    for j=1:m.ss.Ns
        %% Backward information filter update;
        xn1minusfn=G.x_bwd(:,j,t+1)-e.xparticle(:,:,t);
        flbar=m.ss.gamanl'*xn1minusfn;
        hbar=G.x_bwd(:,j,t+1)+m.ss.gamany'*(xn1minusfn)-m.ss.gamalybar'*flbar;
        innovation=m.y(:,t+1)-hbar;
        logC=-(1/2)*sum(innovation.*(m.ss.R0\innovation),1);%x'*Rinv*x == sum(x.*(Delta*x), 1)
        lambda_hat_1_cand=m.ss.Hbar1'/m.ss.R0*innovation;
        if t==m.ss.T-1
            lambda_hat_2_cand=m.ss.Hbar2'/m.ss.R0*innovation;
        else
            lambda_hat_2_cand=G.lambda{t+1}(:,j)+m.ss.Hbar2'/m.ss.R0*innovation;
        end
        %% Backward information filter prediction;
        m1=lambda_hat_2_cand-Omega_hat_22*flbar;
        halfQm=halfQ'*m1;
        tau=sum(flbar.*(Omega_hat_22*flbar),1)-2*sum(lambda_hat_2_cand.*flbar,1)-sum(halfQm.*(M\halfQm),1);
        eta=(m.ss.Albar'-Phi*halfQ/M*halfQ')*m1+lambda_hat_1_cand-Omega_hat_12*flbar;
        logZ=logC-(1/2)*(tau+sum(xn1minusfn.*(m.ss.Qnn\xn1minusfn),1));
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

