% %%����ͼ
% x=0:0.1:1.0;%x���ϵ����ݣ���һ��ֵ�������ݿ�ʼ���ڶ���ֵ��������������ֵ������ֹ
%  a=[0.9574,0.9639,0.9639,0.9638,0.9639,0.9639,0.9639,0.9639,0.9639,0.9638,0.9638]; %a����yֵ
% %  b=[334.4,143.2,297.4,487.2,596.2]; %b����yֵ
%  plot(x,a,'-*b'); %���ԣ���ɫ�����
% axis([0.1,0.9,0.957,0.966])  %ȷ��x����y���ͼ��С
% set(gca,'XTick',[0.1:0.1:0.9]) %x�᷶Χ1-6�����1
% set(gca,'YTick',[0.957:0.001:0.966]) %y�᷶Χ0-700�����100
% legend('AUC');   %���ϽǱ�ע
% xlabel('��')  %x����������
% ylabel('AUC') %y����������


% %%����ͼ
% x=result_LP_KNNC_MKL(:,1);
% y=result_LP_KNNC_MKL(:,2);
% z=result_LP_KNNC_MKL(:,3);
% [xx,yy]=meshgrid(1:1:10,0:0.1:1);
% zz=griddata(result_LP_KNNC_MKL(:,1),result_LP_KNNC_MKL(:,2),result_LP_KNNC_MKL(:,3),xx,yy);
% surf(xx,yy,zz)
% shading interp
% label1 = {'0', '0.1', '0.2', '0.3','0.4','0.5','0.6','0.7','0.8','0.9','1'};%��ǩ
% label2 = {'1', '2', '3', '4','5','6','7','8','9','10'};%��ǩ
% 
% 
% set(gca, 'xticklabel', label2);
% set(gca, 'yticklabel', label1);


%��״ͼ
% ACC��ʾ��ά����
% SNF=zeros(10,10);
figure;h=bar3(SNF);  
for n=1:numel(h)  
    cdata=get(h(n),'zdata');  
    set(h(n),'cdata',cdata,'facecolor','interp')  
end
