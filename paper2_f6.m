clear;
% run codes in "Model simulation programs" to get simulation data
% run codes in "model analysis" to obtain the data for figures
% including optical_flow_model.m, mousePatternDetectionmodel.m

n=10;
orderD=zeros(5999,6,n);
orderS=zeros(5999,6,n);
orderScv=zeros(5999,6,n);
speed_mm=zeros(6,n);
for p=1:6
    for n=1:10  
    cd( 'N:\Junhao\freq_0p5to12');
    load (['result_of_GP_trial' num2str(n) '_50_50_6000_p' num2str(6-p) '_0.5_12.mat'],'vfs','vx','vy');
    for t=1:size(vfs,3)      
          %speed
          speed=vx(:,:,t) + 1i*vy(:,:,t);
          speed_sum=sum(sum(abs(speed)));
          speed_mean=mean(mean(abs(speed)));
          orderD(t,p,n)=abs(sum(sum(speed)))/speed_sum;
          orderS(t,p,n)=mean((abs(reshape(speed,[size(speed,1)*size(speed,2),1]))-speed_mean).^2);          
          orderScv(t,p,n)=sqrt(mean((abs(reshape(speed,[size(speed,1)*size(speed,2),1]))-speed_mean).^2))/speed_mean;          
    end  
    speed_mm(p,n)=mean(mean(mean(abs(vx + 1i*vy))));
    end
end


%% singularity
pn=zeros(6,6,10);% pattern number (time*number)
pnt=zeros(6,6,10);% pattern number (time)
pnn=zeros(6,6,10);% pattern number (number)

ppn=zeros(6,5,10);% proportion of pattern number

for M=1:10%M=6 is 2017 Apr18  
    for exp=1:6
    cd( 'N:\Junhao\freq_0p5to12');
    load (['result_of_GP_trial' num2str(M) '_50_50_6000_p' num2str(6-exp) '_0.5_12r3d2.mat'],'patterns','activeArray');
        for i=1:5
        pnt(exp,i,M)=pnt(exp,i,M)+mean(patterns(find(patterns(:,1)==i),4));
        pnn(exp,i,M)=pnn(exp,i,M)+length(patterns(find(patterns(:,1)==i),4));
        end
    
    pn(exp,1:5,M)=sum(activeArray');
    pn(exp,6,M)=length(find(sum(activeArray)==0));
    pnx=tabulate(cumsum(single((sum(activeArray))>0)));
    b=find(pnx(:,2)>1);
    pnt(exp,6,M)=mean(pnx(b,2));
    pnn(exp,6,M)=length((b));
        for i=1:5
            if i<=4
            ppn(exp,i,M)=length(find(sum(activeArray(3:5,:))==i-1));
            else
            ppn(exp,i,M)=length(find(sum(activeArray(3:5,:))>3));
            end
        end
    end
end

%% plot vfs in certain space 
interval=1;
%raw
rawmatrix=zeros(50,50);
for i=1:interval+1:50
    rawmatrix(i,:)=1;
end
%figure;imagesc(rawmatrix)

%column
colmatrix=zeros(50,50);
for i=1:interval+1:50
    colmatrix(:,i)=1;
end
%figure;imagesc(colmatrix)
intermatrix=rawmatrix & colmatrix;
figure;imagesc(intermatrix)

% load intermedia data for reproduce the results in a quicker way (Deposited data links provided in the paper)
load('fig6a.mat')
load('fig6.mat')

%% Fig.6a,plot figures 
figure
M=3;
p=[5 0];

for i=1:2
    for j=1:3
    load (['./model_data/result_of_p' num2str(p(i)) '_trial' num2str(M) '_50_50_3000.mat'],'vx','vy','Ve');
    subplot(2,3,(i-1)*3+j) 
    itime=200+8*j;
    imagesc(Ve(1:30,1:30,itime))
    caxis([-65 -55])

    % Plot velocity field
    hold on
    ha2 = quiver(vx(1:50,1:50,itime).*intermatrix, vy(1:50,1:50,itime).*intermatrix);
    set(ha2,'color',[0 0 0],'ShowArrowHead','on','MaxHeadSize',0.8,'AutoScale','on', 'AutoScaleFactor', 2)
    grid off
    axis off
    hold off

    end
end

%% mean and sd for 10 trials
xx=zeros(6,6);
for i=1:6
    for j=1:6
        xx(i,j)=(i-0.38)+(j-1)*0.145;
    end
end
xx=[0.66,0.78,0.93,1.055,1.200,1.335;1.66,1.78,1.93,2.055,2.200,2.335;2.66,2.78,2.93,3.055,3.200,3.335;3.66,3.78,3.93,4.055,4.200,4.335;4.66,4.78,4.93,5.055,5.200,5.335;5.66,5.78,5.93,6.055,6.200,6.335];
figure
% subplot(3,2,1)
errorbar(mean(10*speed_mm(:,:),2),std(10*speed_mm(:,:),0,2)./sqrt(10),'-o','LineWidth',1.5)
xlim([0.5,6.5]);
ylabel('Average speed (mm/s)')
xlabel('Brain states from anesthetized to wakefulness')
box off
% set(gca,'xdir','reverse') 
% title('Mouse 1')

% subplot(3,2,3)
errorbar(squeeze(mean(mean(orderScv(:,:,:),1),3)),squeeze(std(mean(orderScv(:,:,:),1),0,3))./sqrt(10),'-o','LineWidth',1.5)%Heterogeneity
xlim([0.5,6.5]);
ylabel('Heterogeneity')
xlabel('Brain states from anesthetized to wakefulness')
box off
% set(gca,'xdir','reverse')

% subplot(3,2,5)
errorbar(squeeze(mean(mean(orderD(:,:,:),1),3)),squeeze(std(mean(orderD(:,:,:),1),0,3))./sqrt(10),'-o','LineWidth',1.5)%Heterogeneity
xlim([0.5,6.5]);
ylabel('Homogeneity')
xlabel('Brain states from anesthetized to wakefulness')
box off
% set(gca,'xdir','reverse') 
% mean duration time of patterns


subplot(3,2,2)
pnt(isnan(pnt))=0;
pnt=pnt/Fs;
seg=[-0.5,-0.3,-0.1,0.1,0.3,0.5];
cols=[0,114,189;
    217,83,23;
    237,177,32;
    126,47,142;
    119,172,48;
    77,190,238]/256;
for mu=1:6
    for type=1:6
        boxplot(squeeze(pnt(mu,type,:)), 'Positions', mu*2+seg(type),'Colors',cols(type,:),'Symbol','+')
        hold on
    end
end
xlim([0.25,13.75])
ylim([-0.1,1])
xticks([2:2:13])
xticklabels({'1','2','3','4','5','6'})

ylabel('Average duration(s)')
box off
% set(gca,'xdir','reverse') 


%times of pattern number/sec
subplot(3,2,4)
pnn=pnn/30;
for mu=1:6
    for type=1:6
        boxplot(squeeze(pnn(mu,type,:)), 'Positions', mu*2+seg(type),'Colors',cols(type,:),'Symbol','+')
        hold on
    end
end
xlim([0.25,13.75])
ylim([0,50])
xticks([2:2:13])
xticklabels({'1','2','3','4','5','6'})
set(gca, 'YScale', 'Linear')
ylabel('Pattern number/sec')
box off


%% fit occupation time, wave number
clear pnnM pnnMr
% subplot(3,2,6)
figure
pnnM=squeeze(mean(ppn(:,:,2:9),3))/(60*100);
for i=1:5
   pnnMr(:,i)=pnnM(:,5+1-i); 
end
bar(pnnMr,'stacked')
legend('>3','3','2','1','0')
% hold on
% y=pnnM(:,5);
% x=1:1:size(pnnM,1); 
% p=fittype('poly1');  
% f=fit(x',y,p); 
% plot(f);
xlim([0.5,6.5]);
ylim([0,1]);
xlabel('Brain states from anesthetized to wakefulness')
ylabel('Proportion of pattern number/frame')
% set(gca,'xdir','reverse') 
box off

%% fig6b insert, need codes in 'model simulation programs'
p = 0.5:-0.1:0;
for id=1:length(p)
    [k,alpha,f]=stability_anesthesia_fieldmodel(p(id),'bottom');
    [TheoryWaveVol(id),TheoryWaveVol_largeWave(id)]=speed_by_LR(k,alpha,f);
end


close all;
fig=figure();
set(gcf,'Position',[100,50,300,300]);
fontsize=11;

plot(p,TheoryWaveVol*10,'o-.','linewidth',1.5);hold on;
plot(p,TheoryWaveVol_largeWave*10,'o-.','linewidth',1.5);hold on;
h=legend('Overall','Large waves');
set(h,'fontsize',fontsize)
xlim([-0.1,0.6]);
ylim([4,32]);
xlabel('Anesthetic degree p','fontsize',fontsize);
ylabel('Wave speed (mm/s)','fontsize',fontsize);
title('Theoretical wave speed','fontsize',fontsize);
set(gca, 'XDir','reverse');
box off;

[~,p1]=Mann_Kendall(TheoryWaveVol,0.05)
[~,p2]=Mann_Kendall(TheoryWaveVol_largeWave,0.05)
