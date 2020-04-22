%研究不同定位算法的定位性能差异
%200*10m的道路环境空旷区域
%移动车辆匀速在道路中间行驶(1m/s)

%%
%1m/s速度行驶，间隔10m，两排
lengh=length(1:1:200)
trace_1=[1:1:200;repmat(5,lengh,1)']'    %真实轨迹,2.5m中间
leng=length(0:15:200)
AP_1=[0:15:200;repmat(0.5,leng,1)']'                                    %AP位置(下面这排)
AP_2=[0:15:200;repmat(9.5,leng,1)']'                                   %AP位置(上面这排)
AP=[AP_1;AP_2]                                                      %全部的AP位置，AP间隔10m部署

%%
%计算轨迹点到AP的距离
for j=1:size(AP(:,1))
    for i=1:size(trace_1(:,1))
        d(i,j)=sqrt((trace_1(i,1)-AP(j,1)).^2+(trace_1(i,2)-AP(j,2)).^2)   %所有轨迹点到第j个AP的真实距离
    end
end    %得到200*42的结构体，每一行代表第i个轨迹点到所有AP的真实距离

%%
%根据输入噪声生成RSSI值
A=2  %衰减因子
rssi=-37.5721-10*1.2*log10(d)    %1m处RSSI为-.37.5721，衰减因子为1.26
rssi_noise=rssi+normrnd(0,2,lengh,leng*2)      %生成方差为2dBm的矩阵,调整这里的噪音

%%
%对数阴影衰减模型，估计距离
intial_test=abs(-37.5721)
distance= 10.^((abs(rssi_noise)-intial_test)/(10 * A))  %所有轨迹点到第j个AP的估计距离

save('AP','AP')
save('onlinedata','rssi','distance')
save('trace_1','trace_1')

%%
clc;
clear;
load('AP.mat')
load('onlinedata.mat')
%最小二乘法应用前先选取4个最近的点
distance_sort=sort(distance,2)  %按行排序
for k=1:length(distance(:,1))
    index(k,:)=find(distance(k,:)<=distance_sort(k,5)) %计算出最小索引的四个数，设置这里的4
end
distance_dim=distance_sort(:,1:5)    %4个最近的点的估计距离，设置这里的4

%4个最近的点的AP坐标位置
AP_X=AP(:,1)  %AP的x轴
AP_Y=AP(:,2)  %AP的y轴
for n=1:length(distance(:,1))    %n为200,n为第n个定位点
    for m=1:length(index(1,:))    %m为4，总共有4列
       mm=index(n,m)   %取出第n行第m列的索引
       AP_X_dim(n,m)=AP_X(mm)   %对于第n个定位点，取出AP的x轴
       AP_Y_dim(n,m)=AP_Y(mm)   %对于第n个定位点，取出AP的y轴
    end
end
AP_Z_dim=zeros(length(AP_Y_dim),size(AP_X_dim,2))  %AP的z轴

%最小二乘法
for p=1:length(distance_dim(:,1))  %d代表第d个轨迹点
    X=AP_X_dim(p,:)
    Y=AP_Y_dim(p,:)
    Z=AP_Z_dim(p,:)
    d_LSE=distance_dim(p,:)'  %将距离从行向量变成列向量
    distance_LSE(:,p)=Algo_LSE(X,Y,Z,d_LSE)
end

%最小二乘+传统泰勒展开
len=size(AP_X_dim,2)    
for i=1:length(distance_dim(:,1)) 
     xtyt_LSE(:,i)=distance_LSE(1:2,i)   %设置泰勒级数展开的初始点，这里是最小二乘法
    
     delta=[0.001;0.001];%初始假设
%      while (abs(delta(1,1))+abs(delta(2,1)))>=0.001%定位误差阈值
k=0
     while k<50 
        for j=1:len
         AP_dim(:,j)=[AP_X_dim(i,j);AP_Y_dim(i,j)]
         d(j)=sqrt(sum((xtyt_LSE(:,i)- AP_dim(:,j)).^2))
         A1(j,:)=[(xtyt_LSE(1,i)-AP_X_dim(i,j))/d(j),(xtyt_LSE(2,i)-AP_Y_dim(i,j))/d(j)];
         B1(j,1)=distance_dim(i,j)-d(j);
        end
        delta=(A1'*A1)\A1'*B1;
        xtyt_LSE(:,i)=xtyt_LSE(:,i)+delta;
        
        k=k+1
     end
end
%save('distance_MLE','distance_MLE')

%%
% clc;
% clear;
% load('AP.mat')
% load('onlinedata.mat')
% %加权最小二乘法应用前先选取4个最近的点
% distance_sort=sort(distance,2)  %按行排序
% for k=1:length(distance(:,1))
%     index(k,:)=find(distance(k,:)<=distance_sort(k,4)) %计算出最小索引的四个数，设置这里的4
% end
% distance_dim=distance_sort(:,1:4)    %4个最近的点的估计距离，设置这里的4
% 
% %4个最近的点的AP坐标位置
% AP_X=AP(:,1)  %AP的x轴
% AP_Y=AP(:,2)  %AP的y轴
% for n=1:length(distance(:,1))    %n为200,n为第n个定位点
%     for m=1:length(index(1,:))    %m为4，总共有4列
%        mm=index(n,m)   %取出第n行第m列的索引
%        AP_X_dim(n,m)=AP_X(mm)   %对于第n个定位点，取出AP的x轴
%        AP_Y_dim(n,m)=AP_Y(mm)   %对于第n个定位点，取出AP的y轴
%     end
% end
% AP_Z_dim=zeros(length(AP_Y_dim),size(AP_X_dim,2))  %AP的z轴

%加权最小二乘法
for p=1:length(distance_dim(:,1))  %d代表第d个轨迹点
    X=AP_X_dim(p,:)
    Y=AP_Y_dim(p,:)
    Z=AP_Z_dim(p,:)
    d_LSE=distance_dim(p,:)'  %将距离从行向量变成列向量
    distance_WLSE(:,p)=Algo_WLSE(X,Y,Z,d_LSE)
end
%save('distance_WLSE','distance_WLSE')

%%
% clc;
% clear;
% load('AP.mat')
% load('onlinedata.mat')
% 
% %极大似然定位算法应用之前先选取5个最近的点
% distance_sort=sort(distance,2)  %按行排序
% for k=1:length(distance(:,1))
%     index(k,:)=find(distance(k,:)<=distance_sort(k,4)) %计算出最小索引的四个数，设置这里的5
% end
% distance_dim=distance_sort(:,1:4)    %5个最近的点的估计距离，设置这里的5

%3个最近的点的AP坐标位置
% AP_X=AP(:,1)  %AP的x轴
% AP_Y=AP(:,2)  %AP的y轴
% for n=1:length(distance(:,1))    %n为200,n为第n个定位点
%     for m=1:length(index(1,:))    %m为3，总共有4列
%        mm=index(n,m)   %取出第n行第m列的索引
%        AP_X_dim(n,m)=AP_X(mm)   %对于第n个定位点，取出AP的x轴
%        AP_Y_dim(n,m)=AP_Y(mm)   %对于第n个定位点，取出AP的y轴
%     end
% end


%极大似然法
% len=size(AP_X_dim,2)    
% for i=1:length(distance_dim(:,1))   
%     for j=1:len-1
%         A(j,:)=[2*(AP_X_dim(i,j)-AP_X_dim(i,len)),2*(AP_Y_dim(i,j)-AP_Y_dim(i,len))]
%         B(j,1)=[distance_dim(i,len)^2-distance_dim(i,j)^2+AP_X_dim(i,j)^2-AP_X_dim(i,len)^2+AP_Y_dim(i,j)^2-AP_Y_dim(i,len)^2]  
%     end
%     distance_MLE(:,i)=(A'*A)\A'*B    
% end     
%  

%改进极大似然法
len=size(AP_X_dim,2)    
for i=1:length(distance_dim(:,1)) 
     xtyt_WMLE(:,i)=distance_MLE(1:2,i)   %设置泰勒级数展开的初始点，这里是最小二乘法
    
     delta=[0.001;0.001];%初始假设
%      while (abs(delta(1,1))+abs(delta(2,1)))>=0.001%定位误差阈值
     k=0
     T=sum(distance_dim(i,:))
     Q=(1/16)*diag(T-distance_dim(i,:))
     while k<100 
        for j=1:len
         AP_dim(:,j)=[AP_X_dim(i,j);AP_Y_dim(i,j)]
         d(j)=sqrt(sum((xtyt_WMLE(:,i)- AP_dim(:,j)).^2))
         A1(j,:)=[(xtyt_WMLE(1,i)-AP_X_dim(i,j))/d(j),(xtyt_WMLE(2,i)-AP_Y_dim(i,j))/d(j)];
         B1(j,1)=distance_dim(i,j)-d(j);
        end
        delta=(A1'*A1)\A1'*Q*B1;
        xtyt_WMLE(:,i)=xtyt_WMLE(:,i)+delta;
        
        k=k+1
     end
end
%%
% clc;
% clear;
% load('AP.mat')
% load('onlinedata.mat')
% 
% %加权质心算法应用之前先选取3个最近的点
% distance_sort=sort(distance,2)  %按行排序
% for k=1:length(distance(:,1))
%     index(k,:)=find(distance(k,:)<=distance_sort(k,4)) %计算出最小索引的四个数，设置这里的3
% end
% distance_dim=distance_sort(:,1:4)    %4个最近的点的估计距离，设置这里的3

%3个最近的点的AP坐标位置
% AP_X=AP(:,1)  %AP的x轴
% AP_Y=AP(:,2)  %AP的y轴
% for n=1:length(distance(:,1))    %n为200,n为第n个定位点
%     for m=1:length(index(1,:))    %m为3，总共有4列
%        mm=index(n,m)   %取出第n行第m列的索引
%        AP_X_dim(n,m)=AP_X(mm)   %对于第n个定位点，取出AP的x轴
%        AP_Y_dim(n,m)=AP_Y(mm)   %对于第n个定位点，取出AP的y轴
%     end
% end

%加权质心算法
len=size(AP_X_dim,2)    
for i=1:length(distance_dim(:,1))  
    x1=0.1;
    y1=0.1;
    dq=0.1;
    for j=1:len
        x1=x1+AP_X_dim(i,j)/distance_dim(i,j);
        y1=y1+AP_Y_dim(i,j)/distance_dim(i,j);
        dq=dq+1/distance_dim(i,j);
    end
    distance_WZX(:,i)= [x1/dq;y1/dq];
end

%加权质心算法+改进泰勒
len=size(AP_X_dim,2)    
for i=1:length(distance_dim(:,1)) 
     xtyt_WWZX(:,i)=distance_WZX(1:2,i)   %设置泰勒级数展开的初始点，这里是最小二乘法
    
     delta=[0.001;0.001];%初始假设
%      while (abs(delta(1,1))+abs(delta(2,1)))>=0.001%定位误差阈值
     k=0
     T=sum(distance_dim(i,:))
     Q=(1/16)*diag(T-distance_dim(i,:))
     while k<100 
        for j=1:len
         AP_dim(:,j)=[AP_X_dim(i,j);AP_Y_dim(i,j)]
         d(j)=sqrt(sum((xtyt_WWZX(:,i)- AP_dim(:,j)).^2))
         A1(j,:)=[(xtyt_WWZX(1,i)-AP_X_dim(i,j))/d(j),(xtyt_WWZX(2,i)-AP_Y_dim(i,j))/d(j)];
         B1(j,1)=distance_dim(i,j)-d(j);
        end
        delta=(A1'*A1)\A1'*Q*B1;
        xtyt_WWZX(:,i)=xtyt_WWZX(:,i)+delta;
        
        k=k+1
     end
end

%%
%画图
load('trace_1.mat')
error_WLSE=sqrt(sum((distance_WLSE(1:2,:)-trace_1').^2))./2   %误差
mean_error_WLSE=mean(error_WLSE)           %平均定位误差
rmse_error_WLSE=(sqrt(mean((distance_WLSE(1,:)-trace_1(:,1)').^2))+sqrt(mean((distance_WLSE(2,:)-trace_1(:,2)').^2)))/2

error_WZX=sqrt(sum((distance_WZX-trace_1').^2))./2   %误差
mean_error_WZX=mean(error_WZX)           %平均定位误差
rmse_error_WZX=(sqrt(mean((distance_WZX(1,:)-trace_1(:,1)').^2))+sqrt(mean((distance_WZX(2,:)-trace_1(:,2)').^2)))/2

error_MLE=sqrt(sum((distance_MLE-trace_1').^2))./2   %误差
mean_error_MLE=mean(error_MLE)           %平均定位误差
rmse_error_MLE=(sqrt(mean((distance_MLE(1,:)-trace_1(:,1)').^2))+sqrt(mean((distance_MLE(2,:)-trace_1(:,2)').^2)))/2

error_xtyt_LSE=sqrt(sum((xtyt_LSE-trace_1').^2))./2   %误差
mean_error_xtyt_LSE=mean(error_xtyt_LSE)           %平均定位误差
rmse_error_xtyt_LSE=(sqrt(mean((xtyt_LSE(1,:)-trace_1(:,1)').^2))+sqrt(mean((xtyt_LSE(2,:)-trace_1(:,2)').^2)))/2

error_xtyt_WWZX=sqrt(sum((xtyt_WWZX-trace_1').^2))./2   %误差
mean_error_xtyt_WWZX=mean(error_xtyt_WWZX)           %平均定位误差
rmse_error_xtyt_WWZX=(sqrt(mean((xtyt_WWZX(1,:)-trace_1(:,1)').^2))+sqrt(mean((xtyt_WWZX(2,:)-trace_1(:,2)').^2)))/2


% figure(1)
% h1=cdfplot(error_MLE)
% hold on
% h2=cdfplot(error_WZX)
% hold on
% h3=cdfplot(error_WLSE)
% hold on
% h4=cdfplot(error_xtyt_LSE)
% hold on
% h5=cdfplot(error_xtyt_WWZX)
% legend('MLE','WCL','WLS','TSLS','WTS-WCL')
% set(gca,'ytick',0:0.1:1)
% set(gca, 'linewidth', 1.0, 'fontsize', 12, 'fontname', 'Times New Roman')
% xlabel('Error (m)')
% ylabel('CDF')
% set(gca,'XLim',[0 10]);%X轴的数据显示范围
% 
% 
% 
