function e=Conditional_MF_smoothing_asyn_cor(G,m)
% G: output of Rao-blackwellized backward simluation for nonlinear sub-state;
% m: model parameters;
fprintf('Conditional_MF_smoothing_asyn_cor\n');
tic;
%% Forward conditional KF;
xlfilter=zeros(m.ss.dimxl,m.ss.Ns,m.ss.T);
Plfilter=zeros(m.ss.dimxl,m.ss.dimxl,m.ss.T);
for t=1:m.ss.T
    if t==1
        xlfilter(:,:,t)=repmat(m.ss.ml0,1,m.ss.Ns);
        Plfilter(:,:,t)=m.ss.Pl0;
    else
        Sigma_nn=m.ss.A(1,2)*Plfilter(:,:,t-1)*m.ss.A(1,2)'+m.ss.Qnn;
        alpha=m.ss.A(1,1)*G.x_bwd(:,:,t-1)+m.ss.A(1,2)*xlfilter(:,:,t-1);
        %% KF time update
        % KF dynamical measurement update;
        L=Plfilter(:,:,t-1)*m.ss.A(1,2)'/Sigma_nn;    
        Xlfilterstar=xlfilter(:,:,t-1)+L*(G.x_bwd(:,:,t)-alpha);
        Plfilterstar=Plfilter(:,:,t-1)-L*Sigma_nn*L';
        % KF time update;
        flbar=m.ss.gamaln*(G.x_bwd(:,:,t)-m.ss.A(1,1)*G.x_bwd(:,:,t-1));
        Xlpre=flbar+m.ss.Albar*Xlfilterstar;
        Plpre=m.ss.Albar*Plfilterstar*m.ss.Albar'+m.ss.Qllbar;
        %% KF measurement update
       Skesai=[Plfilterstar,Plfilterstar*m.ss.Albar';m.ss.Albar*Plfilterstar,Plpre];%PF measurement update;
       %S=Hbar(1,1)*Pkudstar*Hbar(1,1)'+Hbar(1,2)*Pkudstar*Akbar'*Hbar(1,1)'+Hbar(1,1)*Pkudstar*Akbar'*Hbar(1,2)+Hbar(1,2)*Pkpre*Hbar(1,2)'+Qyybar;%PF measurement update;
       M=m.ss.Hbar*Skesai*m.ss.Hbar'+m.ss.R0;%PF measurement update;
       K=Skesai(2,:)*m.ss.Hbar'/M;%KF measurement update;
       Plfilter(:,:,t)=Plpre-K*M*K';
       muy=m.ss.H(1,1)*G.x_bwd(:,:,t)+m.ss.gamany'*(G.x_bwd(:,:,t)-m.ss.A(1,1)*G.x_bwd(:,:,t-1))-m.ss.gamalybar'*flbar+...
           m.ss.Hbar*[Xlfilterstar;Xlpre];
       inovation=m.y(:,t)-muy;
       xlfilter(:,:,t)=Xlpre+K*(inovation);
    end
end
%% Conditional MF smoothing;
for t=m.ss.T:-1:1
    if t==m.ss.T
        Plsmoother(:,:,t)=Plfilter(:,:,t);
        xlsmoother(:,:,t)=xlfilter(:,:,t);
    else
        Plsmoother(:,:,t)=inv(inv(Plfilter(:,:,t))+G.Omega{t});
        for j=1:m.ss.Ns
            xlsmoother(:,j,t)=Plsmoother(:,:,t)*(Plfilter(:,:,t)\xlfilter(:,j,t)+G.lambda{t}(:,j));
        end
    end
    e.xlhat(:,t)=mean(xlsmoother(:,:,t),2);
end
e.time=toc;
e.xlsmoother=xlsmoother;
e.Plsmoother=Plsmoother;

