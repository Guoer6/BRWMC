%tju cs for bioinformatics 
clear;
RESULT=zeros(10,7);
% com(10,10)=0;
% Kr=zeros(10,10);
% for K=1:10
%     K
%     for R=1:10
%         R
%         r=R*0.1;
prediction=zeros(828,314);
for i=1:1
            i
%% 数据集
% load dataSet1;
% load dataSet2;
% load dataSet3;
load interaction;
load disSemSim;
% load interactionA;

%%
% lncRNA_disease_Y = LD_adjmat ;
lncRNA_disease_Y = interaction;
y=lncRNA_disease_Y;
[nm,nd]=size(y);
DS = disSemSim;
LF = lncRNAfunsim(DS,y);

a=0.1;
nfolds = 5 ;
nruns=1;
crossval_idx = crossvalind('Kfold',y(:),nfolds);
    for fold = 1:nfolds
    fold
     y_train = lncRNA_disease_Y;
    test_idx  = find(crossval_idx==fold);
    y_train(test_idx) = 0;
% 	y_train1=WKNKN( y_train,LF,DS, 10, 0.3 );
    [y_train1] = BR(y_train,LF,DS,5,5,0.1,0.9);   %随机游走做预处理
% y_train1= y_train;

    [DG,LG] = gaussiansimilarity(y_train1,nd,nm);
    
%     neighbor_num = floor(a*nd);  
%     DK = KSNS_opt(y_train1',neighbor_num,DS,0);  
%     neighbor_num = floor(a*nm);  
%     LK = KSNS_opt(y_train1,neighbor_num,LF,0);
%%%%%%%%%%%%%%%-q-kernl-%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      [sdistd,sdistm] = SDIST(y_train1,nm,nd);    %归一化
%     %q-kernel function
%     [kd,km] = PQD3(sdistd,sdistm,nm,nd,0.1,0.6);    %Q核函数
%     sd =(1- mapminmax(kd,0,1));
%     qm =(1- mapminmax(km,0,1));
%     QL=(sd+sd')/2;
%     QD=(qm+qm')/2;
    
%%%%%%%%%%%%%%%%%%-相似性网络的融合-%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     beta=0.5;
%     K_COM1=confusion(DS,DG,beta);
%     K_COM2=confusion(LF,LG,beta);
    
%     [K_COM1,K_COM2]=networkfusion(DS,DG,LF,LG);
    
    K_COM1=SNF({DS,DG},4,5,1);
    K_COM2=SNF({LF,LG},4,5,1);
    
%     K_COM1=SKF({DS,DG},20,10,0.1);		
%     K_COM2=SKF({LF,LG},20,10,0.1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     y_train1=WKNKN( y_train,K_COM2,K_COM1, 5, 0.7 );
%     [K_COM1,K_COM2]=CFR(K_COM1,K_COM2);
%     [A]=IMC(y_train1,K_COM2,K_COM1,100); 

%     [A] = BR(y_train1,K_COM2,K_COM1,5,5,0.1,0.9);
%     [y_train1,Rl,Rd] = BR(y_train,K_COM2,K_COM1,5,5,0.1,0.9);   %随机游走做预处理
    [score]=DMC(y_train1,K_COM2,K_COM1);    %双矩阵补全进行预测
%     [score]=DMC(y_train1,Rl,Rd);    %双矩阵补全进行预测
    
    F_1=score';
%     [F_1] = LapRLS(K_COM2,K_COM1,y_train1,2^(-5),20,1);
    
%     [F_1] = BR(y_train1,LF,DS,5,5,0.5,0.9);
%     [F_1]=LP(A,K_COM1,K_COM2,nd,nm,0.6);
%     [F_1]=NCPHLDA(A,A,K_COM1,K_COM2);

    prediction=prediction+F_1;
    y(test_idx)= F_1(test_idx);
%     result = model_evaluate(lncRNA_disease_Y,F_1,y_train);
    end

[X_1,Y_1,tpr,aupr_F_1] = perfcurve(lncRNA_disease_Y(:),y(:),1, 'xCrit', 'reca', 'yCrit', 'prec');
[X,Y,THRE,AUC_F_1,OPTROCPT,SUBY,SUBYNAMES] = perfcurve(lncRNA_disease_Y(:),y(:),1);
result = model_evaluate(lncRNA_disease_Y,F_1,y_train);
% SUM(l2,i)=AUC_F_1;
% SUM(i)=AUC_F_1;
% SUM=SUM';

RESULT(i,:)=result;
MEAN = mean(RESULT());
% AUC=result(1,2);

% if (AUC>0.869 && AUC<0.87)
% 	result = model_evaluate1(lncRNA_disease_Y,F_1,y_train);
%     break;
% end

end

% ave=mean(SUM);

% Kr(K,R)=ave;
%     end
% end
% result = model_evaluate(lncRNA_disease_Y,F_1,y_train);

