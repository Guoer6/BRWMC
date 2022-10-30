function [W]=SNF(Wall,K,t,ALPHA)

if nargin < 2
    K = 4;
end
if nargin < 3
    t = 5;
end

if nargin < 4
    ALPHA = 1;
end

C = length(Wall);
[m,n]=size(Wall{1});
for i = 1 : C
    Wall{i} = Wall{i}./repmat(sum(Wall{i},2),1,n);%归一化
    Wall{i} = (Wall{i} + Wall{i}')/2;%使矩阵变得对称，对应SNF步骤1
end

for i = 1 : C
    newW{i} = FindDominateSet(Wall{i},round(K));%K近邻(newW=S)
end

Wsum = zeros(m,n);
for i = 1 : C
    Wsum = Wsum + Wall{i};%Wsum为两初始矩阵之和
end
for ITER=1:t%迭代
    for i = 1 : C
        Wall0{i}=newW{i}*(Wsum - Wall{i})*newW{i}'/(C-1);%(Wsum - Wall{i})=P,newW=S
    end
    for i = 1 : C
        Wall{i} = BOnormalized(Wall0{i},ALPHA);
    end
    Wsum = zeros(m,n);
    for i = 1 : C
        Wsum = Wsum + Wall{i};
    end

    W = Wsum/C;
    W = W./repmat(sum(W,2),1,n);
    W = (W +W'+eye(n))/2;

    if ITER == 1 || (ITER < t)
        oldW = W;
    else
        break;
    end   
    
end

end

function W = BOnormalized(W,ALPHA)
if nargin < 2 % nargin是用来判断输入变量个数的函数
    ALPHA = 1;
end
W = W+ALPHA*eye(length(W));%加上小倍数的单位矩阵以确保半正定性
W = (W +W')/2;%让W变得对称
end

function newW = FindDominateSet(W,K)%寻找K近邻
[m,n]=size(W);
[YW,IW1] = sort(W,2,'descend');%对矩阵W的每一行进行降序排序
clear YW;
newW=zeros(m,n);
temp=repmat((1:n)',1,K);
I1=(IW1(:,1:K)-1)*m+temp;
newW(I1(:))=W(I1(:));
newW=newW./repmat(sum(newW,2),1,n);
clear IW1;
clear IW2;
clear temp;
end

