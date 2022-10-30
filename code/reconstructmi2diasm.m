
function reas=reconstructmi2diasm(KM,as,alpha)


nm=size(as,1);
nd=size(as,2);
sumas=sum(as);

 for i=1:nd
 if(sumas(1,i)~=0)
    avas(:,i)=as(:,i)/sumas(1,i);   %¡–πÈ“ªªØ
  else avas(:,i)=zeros(nm,1); 
    end
 end



 
reas=zeros(nm,nd);
 
for i=1:nm    
     for j=1:nd
         
         for k=1:nm
               %if i~=k
          reas(i,j)=reas(i,j)+alpha*KM(i,k)*avas(k,j); 
              %end
           end
     end
end

% KM=KM-diag(diag(KM));
% reas=avas*KM;

    reas=as+reas;
end

      