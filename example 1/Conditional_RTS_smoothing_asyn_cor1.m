function e=Conditional_RTS_smoothing_asyn_cor1(G,m)
% G: output of Rao-blackwellized backward simluation for nonlinear sub-state;
% m: model parameters;
fprintf('Conditional_RTS_smoothing_asyn_cor1\n');
tic;
%% Forward conditional KF;
xlfilter=zeros(m.ss.dimxl,m.ss.Ns,m.ss.T);
Plfilter=zeros(m.ss.dimxl,m.ss.dimxl,m.ss.T);
Xlfilterstar=zeros(m.ss.dimxl,m.ss.Ns,m.ss.T);
Plfilterstar=zeros(m.ss.dimxl,m.ss.dimxl,m.ss.T);
Xlpre=zeros(m.ss.dimxl,m.ss.Ns,m.ss.T);
inovation=zeros(m.ss.dimy,m.ss.Ns,m.ss.T);
J1=cell(1,m.ss.T-1);
J2=cell(1,m.ss.T-1);
Delta=cell(1,m.ss.T-1);
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
        Xlfilterstar(:,:,t-1)=xlfilter(:,:,t-1)+L*(G.x_bwd(:,:,t)-alpha);
        Plfilterstar(:,:,t-1)=Plfilter(:,:,t-1)-L*Sigma_nn*L';
        % KF time update;
        flbar=m.ss.gamaln*(G.x_bwd(:,:,t)-m.ss.A(1,1)*G.x_bwd(:,:,t-1));
        Xlpre(:,:,t)=flbar+m.ss.Albar*Xlfilterstar(:,:,t-1);
        Plpre=m.ss.Albar*Plfilterstar(:,:,t-1)*m.ss.Albar'+m.ss.Qllbar;
        %% KF measurement update
       Skesai=[Plfilterstar(:,:,t-1),Plfilterstar(:,:,t-1)*m.ss.Albar';m.ss.Albar*Plfilterstar(:,:,t-1),Plpre];%PF measurement update;
       %S=Hbar(1,1)*Pkudstar*Hbar(1,1)'+Hbar(1,2)*Pkudstar*Akbar'*Hbar(1,1)'+Hbar(1,1)*Pkudstar*Akbar'*Hbar(1,2)+Hbar(1,2)*Pkpre*Hbar(1,2)'+Qyybar;%PF measurement update;
       M=m.ss.Hbar*Skesai*m.ss.Hbar'+m.ss.R0;%PF measurement update;
       K=Skesai(2,:)*m.ss.Hbar'/M;%KF measurement update;
       Plfilter(:,:,t)=Plpre-K*M*K';
       Sigmaly2=Plfilterstar(:,:,t-1)*m.ss.Hbar1'+Plfilterstar(:,:,t-1)*m.ss.Albar'*m.ss.Hbar2';
       J1{t-1}=(Plfilterstar(:,:,t-1)*m.ss.Albar'-Sigmaly2*K')/Plfilter(:,:,t);
       J2{t-1}=(-Plfilterstar(:,:,t-1)*m.ss.Albar'+Sigmaly2*K')/Plfilter(:,:,t)*K+Sigmaly2/M;
       Delta{t-1}=J1{t-1}*m.ss.Albar*Plfilterstar(:,:,t-1)+J2{t-1}*Sigmaly2';
       muy=m.ss.H(1,1)*G.x_bwd(:,:,t)+m.ss.gamany'*(G.x_bwd(:,:,t)-m.ss.A(1,1)*G.x_bwd(:,:,t-1))-m.ss.gamalybar'*flbar+...
           m.ss.Hbar*[Xlfilterstar(:,:,t-1);Xlpre(:,:,t)];
       inovation(:,:,t)=m.y(:,t)-muy;
       xlfilter(:,:,t)=Xlpre(:,:,t)+K*inovation(:,:,t);
    end
end
%% Conditional RTS smoothing;
for t=m.ss.T:-1:1
    if t==m.ss.T
        Plsmoother(:,:,t)=Plfilter(:,:,t);
        xlsmoother(:,:,t)=xlfilter(:,:,t);
    else
        Plsmoother(:,:,t)=Plfilterstar(:,:,t)-Delta{t}+J1{t}*Plsmoother(:,:,t+1)*J1{t}';
        for j=1:m.ss.Ns
            xlsmoother(:,j,t)=Xlfilterstar(:,j,t)+J1{t}*(xlsmoother(:,j,t+1)-Xlpre(:,j,t+1))+J2{t}*inovation(:,j,t+1);
        end
    end
    e.xlhat(:,t)=mean(xlsmoother(:,:,t),2);
end
e.time=toc;
e.xlsmoother=xlsmoother;
e.Plsmoother=Plsmoother;

