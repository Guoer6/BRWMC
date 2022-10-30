% %%折线图
% x=0:0.1:1.0;%x轴上的数据，第一个值代表数据开始，第二个值代表间隔，第三个值代表终止
%  a=[0.9574,0.9639,0.9639,0.9638,0.9639,0.9639,0.9639,0.9639,0.9639,0.9638,0.9638]; %a数据y值
% %  b=[334.4,143.2,297.4,487.2,596.2]; %b数据y值
%  plot(x,a,'-*b'); %线性，颜色，标记
% axis([0.1,0.9,0.957,0.966])  %确定x轴与y轴框图大小
% set(gca,'XTick',[0.1:0.1:0.9]) %x轴范围1-6，间隔1
% set(gca,'YTick',[0.957:0.001:0.966]) %y轴范围0-700，间隔100
% legend('AUC');   %右上角标注
% xlabel('θ')  %x轴坐标描述
% ylabel('AUC') %y轴坐标描述


% %%曲面图
% x=result_LP_KNNC_MKL(:,1);
% y=result_LP_KNNC_MKL(:,2);
% z=result_LP_KNNC_MKL(:,3);
% [xx,yy]=meshgrid(1:1:10,0:0.1:1);
% zz=griddata(result_LP_KNNC_MKL(:,1),result_LP_KNNC_MKL(:,2),result_LP_KNNC_MKL(:,3),xx,yy);
% surf(xx,yy,zz)
% shading interp
% label1 = {'0', '0.1', '0.2', '0.3','0.4','0.5','0.6','0.7','0.8','0.9','1'};%标签
% label2 = {'1', '2', '3', '4','5','6','7','8','9','10'};%标签
% 
% 
% set(gca, 'xticklabel', label2);
% set(gca, 'yticklabel', label1);


%柱状图
% ACC表示二维数组
% SNF=zeros(10,10);
figure;h=bar3(SNF);  
for n=1:numel(h)  
    cdata=get(h(n),'zdata');  
    set(h(n),'cdata',cdata,'facecolor','interp')  
end
