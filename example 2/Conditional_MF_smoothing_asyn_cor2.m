function e=Conditional_MF_smoothing_asyn_cor2(G,m)
% G: output of Rao-blackwellized backward simluation for nonlinear sub-state;
% m: model parameters;
fprintf('Conditional_MF_smoothing_asyn_cor2\n');
tic;
xlfilter=zeros(m.ss.dimxl,m.ss.Ns,m.ss.T);
Plfilter=zeros(m.ss.dimxl,m.ss.dimxl,m.ss.T);
for t=1:m.ss.T
    if t==1
        xlfilter(:,:,t)=repmat(m.ss.ml0,1,m.ss.Ns);
        Plfilter(:,:,t)=m.ss.Pl0;
    else
        Sigma_nn=m.ss.An*Plfilter(:,:,t-1)*m.ss.An'+m.ss.Qnn;
        alpha=G.x_bwd(:,:,t-1)+m.ss.An*xlfilter(:,:,t-1);
        %% KF time update
        % KF dynamical measurement update;
        L=Plfilter(:,:,t-1)*m.ss.An'/Sigma_nn;    
        Xlfilterstar=xlfilter(:,:,t-1)+L*(G.x_bwd(:,:,t)-alpha);
        Plfilterstar=Plfilter(:,:,t-1)-L*Sigma_nn*L';
        % KF time update;
        flbar=m.ss.gamanl'*(G.x_bwd(:,:,t)-G.x_bwd(:,:,t-1));
        Xlpre=flbar+m.ss.Albar*Xlfilterstar;
        Plpre=m.ss.Albar*Plfilterstar*m.ss.Albar'+m.ss.Qllbar;
        %% KF measurement update
       M=m.ss.Hbar1*Plfilterstar*m.ss.Hbar1'+m.ss.Hbar2*m.ss.Albar*Plfilterstar*m.ss.Hbar1'+...
            m.ss.Hbar1*Plfilterstar*m.ss.Albar'*m.ss.Hbar2'+m.ss.Hbar2*Plpre*m.ss.Hbar2'+m.ss.R0;%PF measurement update;
       K=(m.ss.Albar*Plfilterstar*m.ss.Hbar1'+Plpre*m.ss.Hbar2')/M;%KF measurement update;
       Plfilter(:,:,t)=Plpre-K*M*K';
       muy=mea_eq(G.x_bwd(:,:,t))+m.ss.gamany'*(G.x_bwd(:,:,t)-G.x_bwd(:,:,t-1))-m.ss.gamalybar'*flbar+...
           m.ss.Hbar1*Xlfilterstar+m.ss.Hbar2*Xlpre;
       inovation=m.y(:,t)-muy;
       xlfilter(:,:,t)=Xlpre+K*(inovation);
    end
end
xlsmoother=zeros(m.ss.dimxl,m.ss.Ns,m.ss.T);
for t=m.ss.T:-1:1
    if t==m.ss.T
        Plsmoother(:,:,t)=Plfilter(:,:,t);
        xlsmoother(:,:,t)=xlfilter(:,:,t);
    else
        Plinv=inv(Plfilter(:,:,t));
        Plsmoother(:,:,t)=inv(Plinv+G.Omega{t});
        for j=1:m.ss.Ns
            xlsmoother(:,j,t)=Plsmoother(:,:,t)*(Plinv*xlfilter(:,j,t)+G.lambda{t}(:,j));
        end
    end
    e.xlhat(:,t)=mean(xlsmoother(:,:,t),2);
end
e.time=toc;
e.xlsmoother=xlsmoother;

