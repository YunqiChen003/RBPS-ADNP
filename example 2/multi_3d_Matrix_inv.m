function Y=multi_3d_Matrix_inv(X)
% Calculate the inv by page of three-dimensional arrays;
% X: nx*nx*K;
% Y(:,:,i)=inv(X(:,:,1));
i=size(X,1);
k=size(X,3);
M=mat2cell(X,i,i,ones(1,k));
K=cellfun(@inv,M,'UniformOutput',false);
Y=cell2mat(K);
end