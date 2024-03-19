function idx = sample_ancestors(W,N)

% W = W/sum(W);
% u = 1/N*rand;
% 
% idx = zeros(N,1);
% q = 0;
% n = 0;
% for i = 1:N
%     while q < u
%         n = n+1;
%         q = q + W(n);
%     end
%     idx(i) = n;
%     u = u + 1/N;
% end
W = W/sum(W);
bins = cumsum(W);
[~,~,idx] = histcounts(rand(N,1), [0 ; bins]);



