function [interaction]=confusion(W1,W2,beta)

[n,m]=size(W1);
interaction=zeros(n,m);

% for i=1:n
%     for j=1:m
%         if W1(i:j)==0
%             interaction(i,j)=W2(i,j);
%         else
%             interaction(i,j)=W1(i,j);
%         end
%     end
% end

interaction=beta*W1+(1-beta)*W2;

end