function G=MHRBBSi_asyn_cor(e,m)
% e: output of Rao-Blackwellzed particle filters;
% m: model parametes;
fprintf('MHRBBSi_asyn_cor\n');
tic;
% Backward simulator
G.index{m.ss.T}=sample_ancestors(e.W(:,m.ss.T),m.ss.Ns);
G.x_bwd(:,:,m.ss.T) = e.xparticle(:,G.index{m.ss.T},m.ss.T);
G.xnhat(:,m.ss.T)=mean(G.x_bwd(:,:,m.ss.T),2);
G.lambda{m.ss.T}=zeros(m.ss.dimxl,m.ss.Ns);G.Omega{m.ss.T}=zeros(m.ss.dimxl,m.ss.dimxl);
for t = (m.ss.T-1) : (-1) : 1
    %% MH sampling setup;
    % Backward information update;
    G.Omega_hat_11{t+1}=m.ss.Hbar1'/m.ss.R0*m.ss.Hbar1;
    G.Omega_hat_12{t+1}=m.ss.Hbar1'/m.ss.R0*m.ss.Hbar2;
    G.Omega_hat_22{t+1}=G.Omega{t+1}+m.ss.Hbar2'/m.ss.R0*m.ss.Hbar2;
    G.index{t}=e.nIdx(G.index{t+1},t);
    G.x_bwd(:,:,t)=e.xparticle(:,G.index{t},t);
    xn1minusfn=G.x_bwd(:,:,t+1)-G.x_bwd(:,:,t);
    flbar=m.ss.gamanl'*xn1minusfn;
    hbar=mea_eq(G.x_bwd(:,:,t+1))+m.ss.gamany'*(xn1minusfn)-m.ss.gamalybar'*flbar;
    innovation=m.y(:,t+1)-hbar;
    G.logC{t+1}=-(1/2)*sum(innovation.*(m.ss.R0\innovation),1);%x'*Rinv*x == sum(x.*(Delta*x), 1)
    G.lambda_hat_1{t+1}=m.ss.Hbar1'/m.ss.R0*innovation;
    G.lambda_hat_2{t+1}=G.lambda{t+1}+m.ss.Hbar2'/m.ss.R0*innovation;
    % Backward information filter prediction;
    halfQ=chol(m.ss.Qllbar)';
    M=halfQ'*G.Omega_hat_22{t+1}*halfQ+eye(m.ss.dimxl);
    Phi=G.Omega_hat_12{t+1}+m.ss.Albar'*G.Omega_hat_22{t+1};
    Psi=G.Omega_hat_11{t+1}+G.Omega_hat_12{t+1}*m.ss.Albar+m.ss.Albar'*G.Omega_hat_12{t+1}'+m.ss.Albar'*G.Omega_hat_22{t+1}*m.ss.Albar-Phi*halfQ/M*halfQ'*Phi';
    G.Omega{t}=Psi+m.ss.An'/m.ss.Qnn*m.ss.An;
    m1=G.lambda_hat_2{t+1}-G.Omega_hat_22{t+1}*flbar;
    eta=(m.ss.Albar'-Phi*halfQ/M*halfQ')*m1+G.lambda_hat_1{t+1}-G.Omega_hat_12{t+1}*flbar;
    halfQm=halfQ'*m1;
    tau=sum(flbar.*(G.Omega_hat_22{t+1}*flbar),1)-2*sum(G.lambda_hat_2{t+1}.*flbar,1)-sum(halfQm.*(M\halfQm),1);
    G.logZ{t}=G.logC{t+1}-(1/2)*(tau+sum(xn1minusfn.*(m.ss.Qnn\xn1minusfn),1));
    G.lambda{t}=eta+m.ss.An'/m.ss.Qnn*xn1minusfn;
    % Compute predictive PDF;
    halfPfilter=chol(e.Plfilter(:,:,t))';
    Upsilon=halfPfilter'*G.Omega{t}*halfPfilter+eye(m.ss.dimxl);
    halfPlOx=halfPfilter'*(G.lambda{t}-G.Omega{t}*e.xlfilter(:,G.index{t},t));
    kappa=sum(e.xlfilter(:,G.index{t},t).*(G.Omega{t}*e.xlfilter(:,G.index{t},t)),1)-2*sum(G.lambda{t}.*e.xlfilter(:,G.index{t},t),1)-sum(halfPlOx.*(Upsilon\halfPlOx),1);
    G.logpb{t}=G.logZ{t}-1/2*kappa;
    %% MH sampling run;
    J=1;
    while J<=m.ss.rep
        indexproposal=(sample_ancestors(e.W(:,t),m.ss.Ns));
        xproposal=e.xparticle(:,indexproposal,t);
        % Backward information update;
        xn1minusfn_psl=G.x_bwd(:,:,t+1)-xproposal;
        flbar_psl=m.ss.gamanl'*xn1minusfn_psl;
        hbar_psl=mea_eq(G.x_bwd(:,:,t+1))+m.ss.gamany'*(xn1minusfn_psl)-m.ss.gamalybar'*flbar_psl;
        innovation_psl=m.y(:,t+1)-hbar_psl;
        logC_psl=-(1/2)*sum(innovation_psl.*(m.ss.R0\innovation_psl),1);%x'*Rinv*x == sum(x.*(Delta*x), 1)
        lambda_hat_1_psl=m.ss.Hbar1'/m.ss.R0*innovation_psl;
        lambda_hat_2_psl=G.lambda{t+1}+m.ss.Hbar2'/m.ss.R0*innovation_psl;
        % Backward information filter prediction;
        m1_psl=lambda_hat_2_psl-G.Omega_hat_22{t+1}*flbar_psl;
        eta_psl=(m.ss.Albar'-Phi*halfQ/M*halfQ')*m1_psl+lambda_hat_1_psl-G.Omega_hat_12{t+1}*flbar_psl;
        halfQm_psl=halfQ'*m1_psl;
        tau_psl=sum(flbar_psl.*(G.Omega_hat_22{t+1}*flbar_psl),1)-2*sum(lambda_hat_2_psl.*flbar_psl,1)-sum(halfQm_psl.*(M\halfQm_psl),1);
        logZ_psl=logC_psl-(1/2)*(tau_psl+sum(xn1minusfn_psl.*(m.ss.Qnn\xn1minusfn_psl),1));
        lambda_psl=eta_psl+m.ss.An'/m.ss.Qnn*xn1minusfn_psl;
        % Compute predictive PDF;
        halfPlOx_psl=halfPfilter'*(lambda_psl-G.Omega{t}*e.xlfilter(:,indexproposal,t));
        kappa_psl=sum(e.xlfilter(:,indexproposal,t).*(G.Omega{t}*e.xlfilter(:,indexproposal,t)),1)-2*sum(lambda_psl.*e.xlfilter(:,indexproposal,t),1)-sum(halfPlOx_psl.*(Upsilon\halfPlOx_psl),1);
        logpb_psl=logZ_psl-1/2*kappa_psl;
        % Compute MH logratio;
        logratio= logpb_psl-G.logpb{t};
        % Accept (Update) the MH samples;
        u=rand(1,m.ss.Ns);
        [~,c]=find(log(u)<logratio);
        G.index{t}(c)=indexproposal(c);
        G.logpb{t}(:,c)=logpb_psl(:,c);
        G.x_bwd(:,c,t)=xproposal(:,c);
        G.lambda{t}(:,c)=lambda_psl(:,c);
        G.logZ{t}(:,c)=logZ_psl(:,c);
        J=J+1;
    end
    G.xnhat(:,t)=mean(G.x_bwd(:,:,t),2);
end
G.time=toc;

