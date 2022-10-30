function [DD,LL]=networkfusion(DS,DG,LF,LG)
[nl,~]=size(LF);
[nd,~]=size(DS);
DD=zeros(nd,nd);
LL=zeros(nl,nl);
for i=1:nd
    for j=1:nd
        if DS(i,j)~=0
            DD(i,j)=(DS(i,j)+DG(i,j))/2;
        else
            DD(i,j)=DG(i,j);
        end
    end
end
for i=1:nl
    for j=1:nl
        if LF(i,j)~=0
            LL(i,j)=(LF(i,j)+LG(i,j))/2;
        else
            LL(i,j)=LG(i,j);
        end
    end
end

end