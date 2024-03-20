function e=Conditional_MF_smoothing_ind2(G,m)
% G: output of Rao-blackwellized backward simluation for nonlinear sub-state;
% m: model parameters;
fprintf('Conditional_MF_smoothing_ind2\n');
xlfilter=zeros(m.ss.dimxl,m.ss.Ns,m.ss.T);
Plfilter=zeros(m.ss.dimxl,m.ss.dimxl,m.ss.T);
tic;
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
        M=m.ss.H*Plpre*m.ss.H'+m.ss.R;
        K=Plpre*m.ss.H'/M;
        innovation=m.y(:,t)-mea_eq(G.x_bwd(:,:,t));
        Plfilter(:,:,t)=Plpre-K*M*K';
        xlfilter(:,:,t)=Xlpre+K*innovation;
    end 
end
xlsmoother=zeros(m.ss.dimxl,m.ss.Ns,m.ss.T);
tic;
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
