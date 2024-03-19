function G = ffbsi_asyn_cor(M,m)
fprintf('ffbsi_asyn_cor\n');
% Forward filter/backward simulator;
tic;
% Backward simulator
ind=sample_ancestors(M.W(:,m.ss.T),m.ss.Ns);
x_bwd(:,:,m.ss.T) = M.xparticle(:,ind(1:m.ss.Ns),m.ss.T);
G.xhat(:,m.ss.T)=mean(x_bwd(:,:,m.ss.T),2);
for t = (m.ss.T-1) : (-1) : 1
        for j = 1:m.ss.Ns
            fxt =m.ss.A*M.xparticle(:,:,t);
            p=gausspdf2((repmat(x_bwd(:,j,t+1),1,m.ss.N)-fxt),inv(m.ss.Q)).*gausspdf2(repmat(m.y(:,t+1)-(m.ss.H +m.ss.C)*x_bwd(:,j,t+1),1,m.ss.N)+m.ss.C*fxt,inv(m.ss.R0));
            w_sm = p'.*M.W(:,t);
            w_sm = w_sm/sum(w_sm);
            ind = find(rand(1) < cumsum(w_sm),1,'first');
            x_bwd(:,j,t) = M.xparticle(:,ind,t);
        end
        G.xhat(:,t)=mean(x_bwd(:,:,t),2);
end
G.time=toc;
G.x_bwd=x_bwd;
G.Xerror=m.x-G.xhat;

