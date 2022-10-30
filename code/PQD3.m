function [kd,km] = PQD3(sdistd,sdistm,nd,nm,c,q)
% I = interaction ;
% c = 0.01; q = 0.1;0.8931


% sdist=pdist(I,'euclidean');
%(exp(-log(q)/c) - exp(-log(q)/sqrt(sdist+c*c)))/(1-q))
%invqkernel
for i = 1:nm
    for j = 1:nm
%         sdist(i,j)=norm(I(:,i)-I(:,j));
%         km(i,j)= (exp(-log(q)/c)-exp(-log(q)/sqrt(sdist(i,j)+c*c)))/(1-q);
        km(i,j)= (exp(-log(q)/c) - exp(-log(q)/sqrt(sdistm(i,j)+c*c)))/(1-q);
    end
end

for i = 1:nd
    for j = 1:nd
%         sdist(i,j) = norm(I(i,:)-I(j,:));
%         sdist = mapminmax(sdist,0,1);
%         kd(i,j)= (exp(-log(q)/c)-exp(-log(q)/sqrt(sdist(i,j)+c*c)))/(1-q);
        kd(i,j)= (exp(-log(q)/c) - exp(-log(q)/sqrt(sdistd(i,j)+c*c)))/(1-q);
    end
end