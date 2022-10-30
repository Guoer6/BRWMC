%tju cs for bioinformatics 
clear
R=zeros(10,10);
K=20;
for t=1:10
    t
% for l1=1:10
%     l1
%     for l2=1:10
%         l2
% lncRNA_disease_Y = load('known_lncRNA_disease_interaction.txt');    %lcRNA¼²²¡¹ØÁª

load interaction;
lncRNA_disease_Y = interaction;
load disSemSim;
DS = disSemSim;

y_train=lncRNA_disease_Y;
[II,JJ] = find(lncRNA_disease_Y == 1);
[nm,nd]=size(lncRNA_disease_Y);
% DS = load('diseasesimilarity.txt');
LF = lncRNAfunsim(DS,lncRNA_disease_Y);

Pre_value =lncRNA_disease_Y;
weight=[];
for i = 1:length(II)    
    i
    y_train(II(i),JJ(i)) = 0;
%     y_train1=WKNKN( y_train,LF,DS, 5, 0.7);
%     [y_train1] = BR(y_train,LF,DS,5,5,0.1,0.9);
    y_train1=y_train;
    [DG,LG] = gaussiansimilarity(y_train1,nd,nm);
    
    K_COM1=SNF({DS,DG},K,t);
    K_COM2=SNF({LF,LG},K,t);

%     y_train1=WKNKN( y_train,LF,DS, 5, 0.7 );
%     [F_1] = BR(y_train1,K_COM2,K_COM1,8,8,0.4,1);
    [score]=DMC(y_train1,K_COM2,K_COM1);
    F_1=score';
    Pre_value(II(i),JJ(i)) = F_1(II(i),JJ(i));
    y_train = lncRNA_disease_Y;

    if mod(i,500) == 0
        i      
    end
end
% y_train1=WKNKN( y_train,LF,DS, 5, 0.7 );
% [y_train1] = BR(y_train,LF,DS,5,5,0.1,0.9);
y_train1=y_train;
[DG,LG] = gaussiansimilarity(y_train1,nd,nm);

K_COM1=SNF({DS,DG},4,5,1);
K_COM2=SNF({LF,LG},4,5,1);

% y_train1=WKNKN( y_train,LF,DS, 5, 0.7 );
% [F_1] = BR(y_train1,K_COM2,K_COM1,8,8,0.4,1);
[score]=DMC(y_train1,K_COM2,K_COM1);
F_1=score';
for i =1:length(F_1(:))
    if lncRNA_disease_Y(i)==0
        Pre_value(i)=F_1(i);
    end
end
[X_1,Y_1,tpr,aupr_1] = perfcurve(lncRNA_disease_Y(:), Pre_value(:),1, 'xCrit', 'reca', 'yCrit', 'prec');
[X,Y,THRE,AUC_1,OPTROCPT,SUBY,SUBYNAMES] = perfcurve(lncRNA_disease_Y(:), Pre_value(:),1);

result = model_evaluate(lncRNA_disease_Y,F_1,y_train);

% SUM(l2,i)=AUC_1;
% SUM(i)=AUC_1;
%     end
%     a=sum(SUM,2)/10;
%     com(:,l1)=a;
R(t,1)=AUC_1;
end