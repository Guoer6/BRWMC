
function reas=reconstructdi2miasd(KD,as,alpha)
% as已知关联矩阵 KD  39*39

nm=size(as,1); % 行数 39
nd=size(as,2); % 列数 292
sumas=sum(as');% 归一化 

% avas = zeros(1,nd);


 for i=1:nm
 if(sumas(1,i)~=0)
    avas(i,:)=as(i,:)/sumas(1,i);
  else avas(i,:)=zeros(1,nd);  %
    end
 end
%  avas=avas';   %  292*39
 
reas=zeros(nm,nd); %  39*292

for i=1:nm    %39
     for j=1:nd  % 292
         
         for k=1:nm
               if i~=k
         reas(i,j)=reas(i,j)+alpha*KD(i,k)*avas(k,j); 
               end
           end
    end


end
      reas=as+reas;

      