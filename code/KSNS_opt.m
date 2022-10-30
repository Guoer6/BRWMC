function S = KSNS_opt(X,neighbor_num,sim,gpu_flag)
% neighbor_num = single(neighbor_num);
if gpu_flag==1
    X = gpuArray(single(X));
    sim = gpuArray(single(sim));
end
    
feature_matrix = X;
nearst_neighbor_matrix = KSNS_neighbors(sim,neighbor_num,gpu_flag); %和KNN差不多
S = jisuanW_opt(feature_matrix,nearst_neighbor_matrix,gpu_flag);
if gpu_flag==1
    S =  double(gather(S));
end

end

function nearst_neighbor_matrix=KSNS_neighbors(sim,neighbor_num,gpu_flag)
%%%%nearst_neighbor_matrix：  represent Neighbor matrix
  N = size(sim,1);
  D = sim+diag(inf*ones(1,N));  %对角线上的1置为无穷大
  [as, si]=sort(D,2,'ascend');   %按行升序排列
  if gpu_flag==1
      nearst_neighbor_matrix=gpuArray.zeros(N,N);
  else
      nearst_neighbor_matrix = zeros(N,N);  %创建邻居矩阵
  end 
  index=si(:,1:neighbor_num);   %取前K个邻居
  for i=1:N
      nearst_neighbor_matrix(i,index(i,:))=1;   %最近邻的位置置为1，其余为0  
  end
end

function [W,objv] = jisuanW_opt(feature_matrix,nearst_neighbor_matrix,gpu_flag)
lata1 = 1;  lata2 = 1;
X=feature_matrix';  % each column of X represents a sample, and each behavior is an indicator
[~,N] = size(X);    % N represents the number of samples
C = nearst_neighbor_matrix';
rand('state',1);
W = rand(N,N);  %创建一个N维的随机矩阵
if gpu_flag == 1
    W = single(W);
end
W = W- diag(diag(W));   %去掉对角线元素
W = W./repmat(sum(W),N,1); %归一化
G  = jisuan_Kel(X);
G(isnan(G))=0;
G = G/max(G(:));
WC1 = W'*G*W-2*W*G+G;
WC = sum(diag(WC1))/2;
wucha = WC + norm(W.*(1-C),'fro')^2*lata1/2 +  norm(W,'fro')^2*lata2/2;
objv = wucha;
jingdu = 0.0001;
error = jingdu*(1+lata1+lata2);   %Iteration error threshold
we = 1;      %Initial value of error
gen=1;
while  gen<100 && we>error
    %gen
    FZ = G+lata1*C.*W;
    FM = G*W+lata1*W+lata2*W;    
    W = FZ./(FM+eps).*W;  
    WC1 = W'*G*W-2*W*G+G;
    WC = sum(diag(WC1))/2;
    wucha = WC + norm(W.*(1-C),'fro')^2*lata1/2 +  norm(W,'fro')^2*lata2/2;    
    we = abs(wucha-objv(end));
    objv = [objv,wucha];
    gen = gen+1;
end
W=W'; 
W = matrix_normalize(W,gpu_flag);
end

function K  =jisuan_Kel(X)
%X Columns represent samples, and rows represent features
X = X';
sA = (sum(X.^2, 2));
sB = (sum(X.^2, 2));
K = exp(bsxfun(@minus,bsxfun(@minus,2*X*X', sA), sB')/mean(sA));
end

function W = matrix_normalize(W,gpu_flag)
K = 10;
W(isnan(W))=0;
W(1:(size(W,1)+1):end)=0;

if gpu_flag==1
    GW = gpuArray(single(W));
else
    GW = W;
end
for round=1:K
    SW = sum(GW,2);
    ind = find(SW>0);
    SW(ind) = 1./sqrt(SW(ind));
    D1 = diag(SW);
    GW=D1*GW*D1;
end
if gpu_flag==1
    W = gather(GW);
else
    W = GW;
end

end