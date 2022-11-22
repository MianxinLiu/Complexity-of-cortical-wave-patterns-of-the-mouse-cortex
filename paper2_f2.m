clear;


%% Fig2a,anes,histogram of wave direction
cd('/media/user/Elements/NC Data/2017 APR 19 publish');
load('Exp001_Fluo_002_001_sequenceDataFiltered_bandpass0.5_12_box_gp_done.mat','rsDataT','cortexMask','vfsT');
%left hemisphere
thetaLeft=[];
for i=1:size(rsDataT,1)
    for j=1:27
        if(cortexMask(i,j)==1)
           thetaLeft=[thetaLeft angle(vfsT(i,j,:))];
        end
    end
end
thetaLeft(find(thetaLeft==0))=[];
figure;
polarhistogram(2*pi-thetaLeft, 24,'Normalization','probability');

set(gca,'ThetaTick',     [0 90 180 270], ...
        'ThetaTickLabel',{0 '\pi/2' '\pi' '3/2\pi'})
title('Wave propagation directions');
hold on

%right hemisphere
thetaRight=[];
for i=1:size(rsDataT,1)
    for j=28:size(rsDataT,2)
        if(cortexMask(i,j)==1)
           thetaRight=[thetaRight angle(vfsT(i,j,:))];
        end
    end
end
thetaRight(find(thetaRight==0))=[];
polarhistogram(2*pi-thetaRight, 24,'Normalization','probability');
legend('left','right')
set(gca,'ThetaTick',     [0 90 180 270], ...
        'ThetaTickLabel',{0 '\pi/2' '\pi' '3/2\pi'})
    
cd('/media/user/Elements/NC Data/2017 APR 19 publish');
load('Exp001_Fluo_018_001_sequenceDataFiltered_bandpass0.5_12_box_gp_done.mat','rsDataT','cortexMask','vfsT');
%left hemisphere
thetaLeft=[];
for i=1:size(rsDataT,1)
    for j=1:27
        if(cortexMask(i,j)==1)
           thetaLeft=[thetaLeft angle(vfsT(i,j,:))];
        end
    end
end
thetaLeft(find(thetaLeft==0))=[];
figure;
polarhistogram(2*pi-thetaLeft, 24,'Normalization','probability');

set(gca,'ThetaTick',     [0 90 180 270], ...
        'ThetaTickLabel',{0 '\pi/2' '\pi' '3/2\pi'})
title('Wave propagation directions');
hold on

%right hemisphere
thetaRight=[];
for i=1:size(rsDataT,1)
    for j=28:size(rsDataT,2)
        if(cortexMask(i,j)==1)
           thetaRight=[thetaRight angle(vfsT(i,j,:))];
        end
    end
end
thetaRight(find(thetaRight==0))=[];
polarhistogram(2*pi-thetaRight, 24,'Normalization','probability');
legend('left','right')
set(gca,'ThetaTick',     [0 90 180 270], ...
        'ThetaTickLabel',{0 '\pi/2' '\pi' '3/2\pi'})

% fully awake
cd('/media/user/Elements/NC Data/2017 Apr 18 publish');
load('Exp001_Fluo_001_001_sequenceDataFiltered_bandpass0.5_12_box_gp_done.mat','rsDataT','cortexMask','vfsT');
%left hemisphere
thetaLeft=[];
for i=1:size(rsDataT,1)
    for j=1:27
        if(cortexMask(i,j)==1)
           thetaLeft=[thetaLeft angle(vfsT(i,j,:))];
        end
    end
end
thetaLeft(find(thetaLeft==0))=[];
figure;
polarhistogram(2*pi-thetaLeft, 24,'Normalization','probability');

set(gca,'ThetaTick',     [0 90 180 270], ...
        'ThetaTickLabel',{0 '\pi/2' '\pi' '3/2\pi'})
title('Wave propagation directions');
hold on

%right hemisphere
thetaRight=[];
for i=1:size(rsDataT,1)
    for j=28:size(rsDataT,2)
        if(cortexMask(i,j)==1)
           thetaRight=[thetaRight angle(vfsT(i,j,:))];
        end
    end
end
thetaRight(find(thetaRight==0))=[];
polarhistogram(2*pi-thetaRight, 24,'Normalization','probability');
legend('left','right')
set(gca,'ThetaTick',     [0 90 180 270], ...
        'ThetaTickLabel',{0 '\pi/2' '\pi' '3/2\pi'})
    
%% speed from anesthesia to awake state
explabel_a2w={[2,18],[2,7],[1,3],[1,3],[1,6]};
path_a2w={'2017 APR 19 publish','2018Feb 04 M2553M publish','2018Feb 04 M2555M publish','2018Feb 04 M2556M publish','2018Feb 04 M2560F publish'};

path_fa={'2017 Apr 18 publish','2017Aug 05 M2242M publish','2017July24 M2259F publish','2017 July24 M2242M publish'};
explabel_fa={[1:10],[1],[5],[1:7]};


speed_a2w=zeros(5,2);
speed_a2w_std=zeros(5,2);
for M=1:5
    for state=1:2
        exp=explabel_a2w{M}(state);

        cd(['/media/user/Elements/NC Data/' path_a2w{M}]);
        if exp<10
            load (['Exp001_Fluo_00' num2str(exp) '_001_sequenceDataFiltered_bandpass0.5_12_box_gp_done.mat'],'vfsT_a','cortexMask');
        else
             load (['Exp001_Fluo_0' num2str(exp) '_001_sequenceDataFiltered_bandpass0.5_12_box_gp_done.mat'],'vfsT_a','cortexMask');
        end
        speed=reshape(abs(vfsT_a),size(vfsT_a,1)*size(vfsT_a,2),size(vfsT_a,3));
        cortexMaskin=reshape(cortexMask,size(vfsT_a,1)*size(vfsT_a,2),1);
        speed(cortexMaskin==0,:)=[];
        save(['Exp' num2str(exp) '_speed.mat'], 'speed')
        speed_a2w(M,state)=mean(mean(speed));
        speed_a2w_std(M,state)=std(17.4*reshape(speed,[],1));
    end
end

% for fully awake 
speed_fa=zeros(4,1);
speed_fa_std=zeros(4,1);

for M=1:4
    n=length(explabel_fa{M});
    speed_tmp=[];
    for state=1:n
        exp=explabel_fa{M}(state);

        cd(['/media/user/Elements/NC Data/' path_fa{M}]);
        if exp<10
            load (['Exp001_Fluo_00' num2str(exp) '_001_sequenceDataFiltered_bandpass0.5_12_box_gp_done.mat'],'vfsT_a','cortexMask');
        else
             load (['Exp001_Fluo_0' num2str(exp) '_001_sequenceDataFiltered_bandpass0.5_12_box_gp_done.mat'],'vfsT_a','cortexMask');
        end
        speed=reshape(abs(vfsT_a),size(vfsT_a,1)*size(vfsT_a,2),size(vfsT_a,3));
        cortexMaskin=reshape(cortexMask,size(vfsT_a,1)*size(vfsT_a,2),1);
        speed(cortexMaskin==0,:)=[];
        save(['Exp' num2str(exp) '_speed.mat'], 'speed')
        speed_tmp=[speed_tmp, speed];
    end
    speed_fa(M,1)=mean(mean(speed_tmp));
    speed_fa_std(M,1)=std(17.4*reshape(speed_tmp,[],1))/sqrt(n);
end

%% fig b, histogram of wave speed
cd('/media/user/Elements/NC Data/2017 APR 19 publish');
load('Exp2_speed.mat');
speed=reshape(speed,1,[]);
[N,X]=hist(speed*17.4,1000);
plot(X,N);

sptest=[];
sptest=[sptest;speed'];


hold on

cd('/media/user/Elements/NC Data/2017 APR 19 publish');
load('Exp18_speed.mat');
speed=reshape(speed,1,[]);
[N,X]=hist(speed*17.4,1000);
plot(X,N);

sptest=[];
sptest=[sptest;speed'];


set(gca,'Yscale','log')
xlim([0,40])
xlabel('Speed (mm/s)')
ylabel('Histogram (log)')
legend({'Anesthetized','post woken'})

figure;
cd('/media/user/Elements/NC Data/2017 Apr 18 publish');
load('Exp1_speed.mat');
speed=reshape(speed,1,[]);
sptest=[];
sptest=[sptest;speed'];

[N,X]=hist(speed*17.4,1000);
plot(X,N);
set(gca,'Yscale','log')
xlim([0,40])
xlabel('Speed (mm/s)')
ylabel('Histogram (log)')
legend({'Fully awake'})

% friedman test
speed_fri=[];
cd('/media/user/Elements/NC Data/2017 APR 19 publish');
load ('Exp001_Fluo_002_001_sequenceDataFiltered_bandpass0.5_12_box_gp_done.mat','vfsT_a','cortexMask');
speed=reshape(abs(vfsT_a),size(vfsT_a,1)*size(vfsT_a,2),size(vfsT_a,3));
cortexMaskin=reshape(cortexMask,size(vfsT_a,1)*size(vfsT_a,2),1);
speed(cortexMaskin==0,:)=[];
speed=speed*17.4;
speed=reshape(speed,1,[]);
speed_fri=[speed_fri,speed'];

load ('Exp001_Fluo_018_001_sequenceDataFiltered_bandpass0.5_12_box_gp_done.mat','vfsT_a');
speed=reshape(abs(vfsT_a),size(vfsT_a,1)*size(vfsT_a,2),size(vfsT_a,3));
cortexMaskin=reshape(cortexMask,size(vfsT_a,1)*size(vfsT_a,2),1);
speed(cortexMaskin==0,:)=[];
speed=speed*17.4;
speed=reshape(speed,1,[]);
speed_fri=[speed_fri,speed'];

cd('/media/user/Elements/NC Data/2017 Apr 18 publish');
load ('Exp001_Fluo_001_001_sequenceDataFiltered_bandpass0.5_12_box_gp_done.mat','vfsT_a');
speed=reshape(abs(vfsT_a),size(vfsT_a,1)*size(vfsT_a,2),size(vfsT_a,3));
cortexMaskin=reshape(cortexMask,size(vfsT_a,1)*size(vfsT_a,2),1);
speed(cortexMaskin==0,:)=[];
speed=speed*17.4;
speed=reshape(speed,1,[]);
speed_fri=[speed_fri(1:length(speed),:),speed'];

p1 = friedman(speed_fri(:,1:2),1); % anes vs post woken
p2 = friedman(speed_fri(:,[1,3]),1); % anes vs fully awake

%% fig. c
displace=[-0.1,-0.05,0,0.05,0.1];
cols=[0,114,189;
    217,83,23;
    237,177,32;
    126,47,142;
    119,172,48;
    77,190,238;
    162,20,47
    256, 0, 256
    0, 0, 0;
    ]/256;

figure
for M=1:5
   plot(displace(M)+(1:2),17.4*speed_a2w(M,1:2),'-o', 'Color', cols(M,:));
   hold on
   errorbar(displace(M)+(1:2),17.4*speed_a2w(M,1:2),speed_a2w_std(M,1:2), 'Color', cols(M,:))
end
for M=1:4
   plot(displace(M)+2.5, 17.4*speed_fa(M,1),'-o','Color', cols(M+5,:));
   hold on
   errorbar(2.5+displace(M),17.4*speed_fa(M,1),speed_fa_std(M,1), 'Color', cols(M+5,:))
end
% legend('Mouse 1','Mouse 2','Mouse 3','Mouse 4','Mouse 5')
xlabel('Cortical state');
ylabel('Average speed (mm/s)');
xticks([1,2,2.5])
xticklabels({'Anesthetized','Post woken','Fully awake'})

[p,h]= signrank(speed_a2w(:,1),speed_a2w(:,2),'tail','right')

%% order parameter on direction and speed
for M=1:5
    orderD=zeros(26921,2);
    orderS=zeros(26921,2);
    orderScv=zeros(26921,2);

    for state=1:2 %1 anes 2 awake
        exp=explabel_a2w{M}(state);
        cd(['/media/user/Elements/NC Data/' path_a2w{M}]);
        if exp<10
            load (['Exp001_Fluo_00' num2str(exp) '_001_sequenceDataFiltered_bandpass0.5_12_box_gp_done.mat'],'vfsT_a','cortexMask');
        else
             load (['Exp001_Fluo_0' num2str(exp) '_001_sequenceDataFiltered_bandpass0.5_12_box_gp_done.mat'],'vfsT_a','cortexMask');
        end
        for t=1:size(vfsT_a,3)      
              %speed
              speed_incortex2=vfsT_a(:,:,t).*cortexMask;
              speed_incortex2(find(cortexMask~=1))=[];
              speed_sum=sum(abs(speed_incortex2));
              orderD(t,state)=abs(sum(speed_incortex2))/speed_sum;
              orderS(t,state)=mean((abs(speed_incortex2)-(speed_sum/size(speed_incortex2,2))).^2)/(speed_sum/size(speed_incortex2,2));          
              orderScv(t,state)=sqrt(mean((abs(speed_incortex2)-(speed_sum/size(speed_incortex2,2))).^2))/(speed_sum/size(speed_incortex2,2));          
        end  
    end
    save order.mat orderD orderS orderScv
end 

% for fully awake 
for M=3:3
    n=length(explabel_fa{M});
    orderD=zeros(26921,n);
    orderS=zeros(26921,n);
    orderScv=zeros(26921,n);

    for state=1:n
        exp=explabel_fa{M}(state);
        cd(['/media/user/Elements/NC Data/' path_fa{M}]);
        if exp<10
            load (['Exp001_Fluo_00' num2str(exp) '_001_sequenceDataFiltered_bandpass0.5_12_box_gp_done.mat'],'vfsT_a','cortexMask');
        else
             load (['Exp001_Fluo_0' num2str(exp) '_001_sequenceDataFiltered_bandpass0.5_12_box_gp_done.mat'],'vfsT_a','cortexMask');
        end
        for t=1:size(vfsT_a,3)      
              %speed
              speed_incortex2=vfsT_a(:,:,t).*cortexMask;
              speed_incortex2(find(cortexMask~=1))=[];
              speed_sum=sum(abs(speed_incortex2));
              orderD(t,state)=abs(sum(speed_incortex2))/speed_sum;
              orderS(t,state)=mean((abs(speed_incortex2)-(speed_sum/size(speed_incortex2,2))).^2)/(speed_sum/size(speed_incortex2,2));          
              orderScv(t,state)=sqrt(mean((abs(speed_incortex2)-(speed_sum/size(speed_incortex2,2))).^2))/(speed_sum/size(speed_incortex2,2));          
        end  
    end
    save order.mat orderD orderS orderScv
end 

%% re-load order parameter results 

all_order_a2w=zeros(5,6);
all_order_a2w_std=zeros(5,6);

for M=1:5
    cd(['/media/user/Elements/NC Data/' path_a2w{M}]);
    load ('order.mat');
    
    all_order_a2w(M,1)=mean(orderD(:,1));
    all_order_a2w(M,2)=mean(orderD(:,2));
    all_order_a2w(M,3)=mean(orderScv(:,1));
    all_order_a2w(M,4)=mean(orderScv(:,2));
    all_order_a2w(M,5)=mean(orderS(:,1));
    all_order_a2w(M,6)=mean(orderS(:,2));
    
    all_order_a2w_std(M,1)=std(orderD(:,1));
    all_order_a2w_std(M,2)=std(orderD(:,2));
    all_order_a2w_std(M,3)=std(orderScv(:,1));
    all_order_a2w_std(M,4)=std(orderScv(:,2));
    all_order_a2w_std(M,5)=std(orderS(:,1));
    all_order_a2w_std(M,6)=std(orderS(:,2));
end

all_order_fa=zeros(4,3);
all_order_fa_std=zeros(4,3);

for M=1:4
    n=length(explabel_fa{M});
    cd(['/media/user/Elements/NC Data/' path_fa{M}]);
    load ('order.mat');

    all_order_fa(M,1)=mean(mean(orderD));
    all_order_fa(M,2)=mean(mean(orderScv));
    all_order_fa(M,3)=mean(mean(orderS));
    
    all_order_fa_std(M,1)=std(reshape(orderD,[],1))/sqrt(n);
    all_order_fa_std(M,2)=std(reshape(orderScv,[],1))/sqrt(n);
    all_order_fa_std(M,3)=std(reshape(orderS,[],1))/sqrt(n);
end


%% fig 2d
figure
for M=1:5
   plot(displace(M)+(1:2),all_order_a2w(M,1:2),'-o', 'Color', cols(M,:));
   hold on
   errorbar(displace(M)+(1:2),all_order_a2w(M,1:2),all_order_a2w_std(M,1:2), 'Color', cols(M,:))
end
for M=1:4
   plot(2.5+displace(M),all_order_fa(M,1),'-o', 'Color', cols(M+5,:));
   hold on
   errorbar(2.5+displace(M),all_order_fa(M,1),all_order_fa_std(M,1), 'Color', cols(M+5,:))
   all_order_fa_std(M,1)
end
% legend('Mouse 1','Mouse 2','Mouse 3','Mouse 4','Mouse 5')
xlabel('cortical state');
ylabel('Homogeneity');
xticks([1,2,2.5])
xticklabels({'Anesthetized','Post woken','Fully awake'})

% statistical test
[p,h]= signrank(all_order_a2w(:,1),all_order_a2w(:,2),'tail','right')

%% fig 2e
figure
for M=1:5
   plot(displace(M)+(1:2),all_order_a2w(M,3:4),'-o', 'Color', cols(M,:));
   hold on
   errorbar(displace(M)+(1:2),all_order_a2w(M,3:4),all_order_a2w_std(M,3:4), 'Color', cols(M,:))
end
for M=1:4
   plot(2.5+displace(M),all_order_fa(M,2),'-o', 'Color', cols(M+5,:));
   hold on
   errorbar(2.5+displace(M),all_order_fa(M,2),all_order_fa_std(M,2), 'Color', cols(M+5,:))
end
% legend('Mouse 1','Mouse 2','Mouse 3','Mouse 4','Mouse 5')
xlabel('cortical state');
ylabel('Heterogeneity');
xticks([1,2,2.5])
xticklabels({'Anesthetized','Post woken','Fully awake'})
ylim([0.3,1])

[p,h]= signrank(all_order_a2w(:,3),all_order_a2w(:,4),'tail','left')


