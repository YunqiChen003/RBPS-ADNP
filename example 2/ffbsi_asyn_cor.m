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
            fxt=sys_eq(M.xparticle(:,:,t),m);
            xminusfxt=x_bwd(:,j,t+1)-fxt;
            innovation=m.y(:,t+1)-mea_eq(x_bwd(:,j,t+1))-m.ss.Qxy'/m.ss.Qxx*xminusfxt;
            p=gausspdf2(xminusfxt,inv(m.ss.Qxx)).*gausspdf2(innovation,inv(m.ss.R0));
            w_sm = p'.*M.W(:,t);
            w_sm = w_sm/sum(w_sm);
            ind = find(rand(1) < cumsum(w_sm),1,'first');
            x_bwd(:,j,t) = M.xparticle(:,ind,t);
        end
        G.xhat(:,t)=mean(x_bwd(:,:,t),2);
end
G.Xerror=m.x-G.xhat;
G.time=toc;
G.x_bwd=x_bwd;

