function [Rl]=BR(A,Wrr,Wdd,L1,L2,alpha,beta)
% beta=0;
normWrr = normFun1(Wrr);
normWdd = normFun1(Wdd);
R0=A;
Rt=R0;

% alpha=0.9;
% beta=0.1;

sl=sum(A,1);
sc=sum(A,2);
[nd,nl]=size(A);
for j = 1:nl
	text1(1,j)=norm(A(:,j));
end
z1=find(text1==0);
for i = 1:nd
	text2(i,1)=norm(A(i,:));
end
z2=find(text2==0);

% alpha=0.1;
% knownerr = reconstructmi2diasm(Wrr,A,alpha);
% knownerd = reconstructdi2miasd(Wdd,A,alpha);
% [SD,SL]=NCPHLDA(knownerr, knownerd ,Wdd,Wrr);

% [SLP]=NCPHLDA(A,A,Wdd,Wrr);
% [SLP]=LP(A,Wdd,Wrr,nl,nd,0.6);

%bi-random walk on the heterogeneous network
t=0;
delta=1;
%bi-random walk on the heterogeneous network
% for t=1:max(L1,L2)
    
%     ftl = 0;
%     ftr = 0;
%    beta=0;
    %random walk on the lncRNA similarity network
    while  (delta > 1e-6)
%         nRtleft =(1- alpha) * normWrr * (beta*SLP+(1-beta)*Rt) + alpha*R0;
%         nRtleft =(1- gama) * normWrr * Rt + gama * (beta*SD+(1-beta)*R0);
        nRtleft =(1- alpha) * normWrr *Rt + alpha*R0;
        delta =abs(sum(sum((abs(nRtleft)-abs(Rt)))));
        Rt = nRtleft;
        t=t+1;
%         ftl = 1;
    end
    %random walk on the disease similarity network
    delta=1;
    while  (delta > 1e-6)
%         nRtright = (1-alpha) * (beta*SLP+(1-beta)*Rt) * normWdd + alpha * R0;
%         nRtright = (1-gama) * Rt * normWdd + gama * (beta*SL+(1-beta)*R0);
        nRtright = (1-alpha) * Rt * normWdd + alpha*R0;
        delta =abs(sum(sum((abs(nRtright)-abs(Rt)))));
        Rt = nRtright;
        t=t+1;
%         ftr = 1;
    end
    %Rt: predictive association scores between each lncRNA-disease pair
    Rl = nRtleft;
    Rd = nRtright;
    Rt =  (nRtleft + nRtright)/2;
end
