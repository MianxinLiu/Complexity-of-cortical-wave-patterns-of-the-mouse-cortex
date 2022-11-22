clear
%% Fig. 1A
%% power specturm on every point
% rsData is the spatial resized voltage signals
% cortexMask is also spatial resized, cortexMask==1 means inside the cortex
explabel={[2,18],[2,7],[1,3],[1,3],[1,6]};
path={'2017 APR 19 publish','2018Feb 04 M2553M publish','2018Feb 04 M2555M publish','2018Feb 04 M2556M publish','2018Feb 04 M2560F publish'};

path_fa={'2017 Apr 18 publish','2017Aug 05 M2242M publish','2017July24 M2259F publish','2017 July24 M2242M publish'};
explabel_fa={[1:10],[1],[5],[1:7]};

P1m_aw=zeros(2,5,13463);
for M=1:5
    for state=1:2
        exp=explabel{M}(state);
        cd(['/media/user/Elements/NC Data/' path{M}]);
        if exp<10
            load (['Exp001_Fluo_00' num2str(exp) '_001_sequenceDataFiltered_bandpass0.5_12_box_gp_done.mat'],'rsDataT','cortexMask');
        else
             load (['Exp001_Fluo_0' num2str(exp) '_001_sequenceDataFiltered_bandpass0.5_12_box_gp_done.mat'],'rsDataT','cortexMask');
        end
        L=size(rsDataT,3);
        Fs=150;
        for i=1:size(rsDataT,1)
            for j=1:size(rsDataT,2)
                if cortexMask(i,j)==1
                Y = fft(squeeze(rsDataT(i,j,:)));   
                P2 = abs(Y);
                P1(i,j,:) = P2(1:L/2+1).^2/L;
                P1(i,j,2:end-1) = 2*P1(i,j,2:end-1);
                f = Fs*(0:(L/2))/L;
                end
            end
        end
        %% moving average of power specturm
        P1m=length(f);
        a=squeeze(mean(mean(P1,1),2));% mean power specturm
        n=20;
        for i=1:length(a)
            if(i-n<=0)
                P1m(i,1)=mean(a(1:i+n,1));
            end
            if(i+n>=length(a))
                P1m(i,1)=mean(a(i-n:length(a),1));
            end

            if(i-n>0&&i+n<length(a))
                P1m(i,1)= mean(a(i-n:i+n,1));
            end
        end
        P1m_aw(state,M,:)=P1m;
    end
end

P1m_fw=zeros(19,13463);
count=1;
for M=1:4
    n=length(explabel_fa{M});
    for state=1:n
        exp=explabel_fa{M}(state);
        cd(['/media/user/Elements/NC Data/' path_fa{M}]);
        if exp<10
            load (['Exp001_Fluo_00' num2str(exp) '_001_sequenceDataFiltered_bandpass0.5_12_box_gp_done.mat'],'rsDataT','cortexMask');
        else
             load (['Exp001_Fluo_0' num2str(exp) '_001_sequenceDataFiltered_bandpass0.5_12_box_gp_done.mat'],'rsDataT','cortexMask');
        end
        L=size(rsDataT,3);
        Fs=150;
        for i=1:size(rsDataT,1)
            for j=1:size(rsDataT,2)
                if cortexMask(i,j)==1
                Y = fft(squeeze(rsDataT(i,j,:)));   
                P2 = abs(Y);
                P1(i,j,:) = P2(1:L/2+1).^2/L;
                P1(i,j,2:end-1) = 2*P1(i,j,2:end-1);
                f = Fs*(0:(L/2))/L;
                end
            end
        end
        %% moving average of power specturm
        P1m=length(f);
        a=squeeze(mean(mean(P1,1),2));% mean power specturm
        n=20;
        for i=1:length(a)
            if(i-n<=0)
                P1m(i,1)=mean(a(1:i+n,1));
            end
            if(i+n>=length(a))
                P1m(i,1)=mean(a(i-n:length(a),1));
            end

            if(i-n>0&&i+n<length(a))
                P1m(i,1)= mean(a(i-n:i+n,1));
            end
        end
        P1m_fw(count,:)=P1m;
        count=count+1;
    end
end


figure
plot(f,mean(squeeze(P1m_aw(1,:,:)),1),'LineWidth',1.5)
hold on
plot(f,mean(squeeze(P1m_aw(2,:,:)),1),'LineWidth',1.5)
plot(f,mean(P1m_fw,1),'-g','LineWidth',1.5)

legend('Anesthetized','Post woken','Fully awake', 'AutoUpdate','off');

shadedErrorBar(f,mean(squeeze(P1m_aw(1,:,:)),1),std(squeeze(P1m_aw(1,:,:)),1)./sqrt(5),'lineprops','-b')
shadedErrorBar(f,mean(squeeze(P1m_aw(2,:,:)),1),std(squeeze(P1m_aw(2,:,:)),1)./sqrt(5),'lineprops','-r')
shadedErrorBar(f,mean(P1m_fw,1),std(P1m_fw,1)./sqrt(19),'lineprops','-g')

set(gca, 'Xlim', [0 20])
xlabel('Frequency(Hz)')
ylabel('Power')


%% Fig. 1B  Time-frequency spectrum of voltage signals, wavelet
%%anes
cd('/media/user/Elements/NC Data/2017 APR 19 publish');
load('Exp001_Fluo_002_001_sequenceDataFiltered_bandpass0.5_12_box_gp_done.mat','rsDataT','cortexMask','vfsT');
rsData=rsDataT(:,:,9001:18000);
vfs=vfsT(:,:,9001:18000);

%% fig1c
figure;
set (gcf,'Position',[200,100,420*0.3,340*0.3], 'color','w')
itime=3250+130;
useAmplitude=false;

signalGrid = rsData(:,:,itime);
sigLims = [-0.005 0.005];

% Apply mask to data
cortexMask(find(cortexMask<0.7))=NaN;
cortexMask(find(cortexMask>=0.7))=1;
signalGrid = signalGrid .* cortexMask;
% Plot signal grid
imagesc(signalGrid, sigLims)
colormap(darkb2r(-0.005,0.005))
h=imagesc(signalGrid, sigLims);
set(h,'alphadata',~isnan(signalGrid))
hold on
line([19,19],[2,40],'color','k')
hold on
plot(19,6,'ro')
hold on
plot(19,22,'go')
hold on
plot(19,36,'bo')
grid off
axis off
hold off
%% Fig1b 
figure 
[cfs,frequencies] = cwt(rsData(36,19,3250:4000),150./[0.1:0.01:10],'cmor1-1',1/150);
t = linspace(0,2,300);
imagesc(t,frequencies,abs(cfs(:,1:300))); 
grid off;
set(gca,'YDir','reverse');
xlabel('Time(s)'); 
ylabel('Frequency(Hz)');
shading interp
colormap(parula)
%colorbar('south')

%% surf of the amplitude,fig.1 d
cortexMask_in=zeros(size(cortexMask,1),size(cortexMask,2));
index=0;
for r=1:size(cortexMask,1)
    for c=1:size(cortexMask,2)
        index=index+1;
        cortexMask_in(r,c)=index;                         
    end
end

cortex_index=reshape(cortexMask_in,[size(rsData,1)*size(rsData,2) 1]);

array2d=zeros(size(cortexMask,1)*size(cortexMask,2),size(rsData,3));
for k=1:size(rsData,3)
    array=reshape(rsData(:,:,k),[size(rsData,1)*size(rsData,2) 1]);
    array2d(:,k)=array;
end
figure
set(groot,'DefaultFigureColormap',jet)
imagesc(array2d(794:832,3250:3550))% 794-832 col: 19
colormap(darkb2r(-0.005,0.005))
shading interp
%% plot oscillation of three different pixels
%  global signal
ys=zeros(1,size(rsData,3)-1);
for i=1:size(rsData,3)
    ys(i)=sum(nansum(rsData(:,:,i).*cortexMask))/sum(nansum(cortexMask));
end

figure
plot(ys(3250:3550));
legend('average voltage amplitude')
hold on
plot(squeeze(rsData(6,19,3250:3550)));
hold on
plot(squeeze(rsData(22,19,3250:3550)));
hold on
plot(squeeze(rsData(36,19,3250:3550)));
hold on
line([0,300],[0,0])
hold on
line([149,149],[-0.005,0.005],'color','k')
% hold on
% line([162,162],[-0.005,0.005],'color','k')
hold on
line([175,175],[-0.005,0.005],'color','k')
% hold on
% line([130,130],[-0.005,0.005],'color','k')

%% plot vfs in certain space 
interval=1;
%raw
rawmatrix=zeros(size(cortexMask,1),size(cortexMask,2));
for i=1:interval+1:size(cortexMask,1)
    rawmatrix(i,:)=1;
end
%figure;imagesc(rawmatrix)

%column
colmatrix=zeros(size(cortexMask,1),size(cortexMask,2));
for i=1:interval+1:size(cortexMask,2)
    colmatrix(:,i)=1;
end
%figure;imagesc(colmatrix)
intermatrix=rawmatrix & colmatrix;
%figure;imagesc(intermatrix)
%% plot phase velocity field, amplitude and mask
%%149
figure;
set (gcf,'Position',[200,100,420*0.65,340*0.65], 'color','w')
itime=3250+149;
useAmplitude=false;

signalGrid = rsData(:,:,itime);
sigLims = [-0.005 0.005];

% Apply mask to data
cortexMask(find(cortexMask<0.7))=NaN;
cortexMask(find(cortexMask>=0.7))=1;
signalGrid = signalGrid .* cortexMask;
% Plot signal grid
imagesc(signalGrid, sigLims)
colormap(darkb2r(-0.005,0.005))
h=imagesc(signalGrid, sigLims);
set(h,'alphadata',~isnan(signalGrid))

% Plot velocity field
hold on
vf = vfs(:,:,itime) * vfScale.*cortexMask;
ha2 = quiver(real(vf).*intermatrix, imag(vf).*intermatrix);
set(ha2,'color',[0 0 0],'ShowArrowHead','on','MaxHeadSize',0.8,'AutoScale','on', 'AutoScaleFactor', 2)
set(gca,'YDir','reverse');
%title(num2str(itime));
grid off
axis off
hold off

%% 175
figure;
set (gcf,'Position',[200,100,420*0.65,340*0.65], 'color','w')
itime=3250+175;
useAmplitude=false;

signalGrid = rsData(:,:,itime);
sigLims = [-0.005 0.005];

% Apply mask to data
cortexMask(find(cortexMask<0.7))=NaN;
cortexMask(find(cortexMask>=0.7))=1;
signalGrid = signalGrid .* cortexMask;
% Plot signal grid
imagesc(signalGrid, sigLims)
colormap(darkb2r(-0.005,0.005))
h=imagesc(signalGrid, sigLims);
set(h,'alphadata',~isnan(signalGrid))

% Plot velocity field
hold on
vf = vfs(:,:,itime) * vfScale.*cortexMask;
ha2 = quiver(real(vf).*intermatrix, imag(vf).*intermatrix);
set(ha2,'color',[0 0 0],'ShowArrowHead','on','MaxHeadSize',0.8,'AutoScale','on', 'AutoScaleFactor', 2)
set(gca,'YDir','reverse');
%title(num2str(itime));
grid off
axis off
hold off

%% post woken
cd('/media/user/Elements/NC Data/2017 APR 19 publish');
load('Exp001_Fluo_018_001_sequenceDataFiltered_bandpass0.5_12_box_gp_done.mat','rsDataT','cortexMask','vfsT');

rsData=rsDataT(:,:,9001:18000);
vfs=vfsT(:,:,9001:18000);
%%Fig1b 
figure 
[cfs,frequencies] = cwt(rsData(36,19,3250:4000),150./[0.1:0.01:10],'cmor1-1',1/150);
t = linspace(0,2,300);
imagesc(t,frequencies,abs(cfs(:,1:300))); 
grid off;
set(gca,'YDir','reverse');
xlabel('Time(s)'); 
ylabel('Frequency(Hz)');
shading interp
colormap(parula)
%colorbar('south')

%% surf of the amplitude,fig.1 d
cortexMask_in=zeros(size(cortexMask,1),size(cortexMask,2));
index=0;
for r=1:size(cortexMask,1)
    for c=1:size(cortexMask,2)
        index=index+1;
        cortexMask_in(r,c)=index;                         
    end
end

cortex_index=reshape(cortexMask_in,[size(rsData,1)*size(rsData,2) 1]);

array2d=zeros(size(cortexMask,1)*size(cortexMask,2),size(rsData,3));
for k=1:size(rsData,3)
    array=reshape(rsData(:,:,k),[size(rsData,1)*size(rsData,2) 1]);
    array2d(:,k)=array;
end
figure
set(groot,'DefaultFigureColormap',jet)
imagesc(array2d(794:832,3250:3550))% 794-832 col: 19
colormap(darkb2r(-0.005,0.005))
shading interp
%% plot oscillation of three different pixels
%  global signal
ys=zeros(1,size(rsData,3)-1);
for i=1:size(rsData,3)
    ys(i)=sum(nansum(rsData(:,:,i).*cortexMask))/sum(nansum(cortexMask));
end

figure
plot(ys(3250:3550));
legend('average voltage amplitude')
hold on
plot(squeeze(rsData(6,19,3250:3550)));
hold on
plot(squeeze(rsData(22,19,3250:3550)));
hold on
plot(squeeze(rsData(36,19,3250:3550)));
hold on
line([0,300],[0,0])
hold on
line([130,130],[-0.005,0.005],'color','k')
% hold on
% line([162,162],[-0.005,0.005],'color','k')
hold on
line([175,175],[-0.005,0.005],'color','k')
% hold on
% line([130,130],[-0.005,0.005],'color','k')

%% plot vfs in certain space 
interval=1;
%raw
rawmatrix=zeros(size(cortexMask,1),size(cortexMask,2));
for i=1:interval+1:size(cortexMask,1)
    rawmatrix(i,:)=1;
end
%figure;imagesc(rawmatrix)

%column
colmatrix=zeros(size(cortexMask,1),size(cortexMask,2));
for i=1:interval+1:size(cortexMask,2)
    colmatrix(:,i)=1;
end
%figure;imagesc(colmatrix)
intermatrix=rawmatrix & colmatrix;
%figure;imagesc(intermatrix)
%% plot phase velocity field, amplitude and mask
%%149
figure;
set (gcf,'Position',[200,100,420*0.65,340*0.65], 'color','w')
itime=3250+130;
useAmplitude=false;

signalGrid = rsData(:,:,itime);
sigLims = [-0.005 0.005];

% Apply mask to data
cortexMask(find(cortexMask<0.7))=NaN;
cortexMask(find(cortexMask>=0.7))=1;
signalGrid = signalGrid .* cortexMask;
% Plot signal grid
imagesc(signalGrid, sigLims)
colormap(darkb2r(-0.005,0.005))
h=imagesc(signalGrid, sigLims);
set(h,'alphadata',~isnan(signalGrid))

% Plot velocity field
hold on
vf = vfs(:,:,itime) * vfScale.*cortexMask;
ha2 = quiver(real(vf).*intermatrix, imag(vf).*intermatrix);
set(ha2,'color',[0 0 0],'ShowArrowHead','on','MaxHeadSize',0.8,'AutoScale','on', 'AutoScaleFactor', 2)
set(gca,'YDir','reverse');
%title(num2str(itime));
grid off
axis off
hold off

%% 175
figure;
set (gcf,'Position',[200,100,420*0.65,340*0.65], 'color','w')
itime=3250+175;
useAmplitude=false;

signalGrid = rsData(:,:,itime);
sigLims = [-0.005 0.005];

% Apply mask to data
cortexMask(find(cortexMask<0.7))=NaN;
cortexMask(find(cortexMask>=0.7))=1;
signalGrid = signalGrid .* cortexMask;
% Plot signal grid
imagesc(signalGrid, sigLims)
colormap(darkb2r(-0.005,0.005))
h=imagesc(signalGrid, sigLims);
set(h,'alphadata',~isnan(signalGrid))

% Plot velocity field
hold on
vf = vfs(:,:,itime) * vfScale.*cortexMask;
ha2 = quiver(real(vf).*intermatrix, imag(vf).*intermatrix);
set(ha2,'color',[0 0 0],'ShowArrowHead','on','MaxHeadSize',0.8,'AutoScale','on', 'AutoScaleFactor', 2)
set(gca,'YDir','reverse');
%title(num2str(itime));
grid off
axis off
hold off
%% fully awake

cd('/media/user/Elements/NC Data/2017 Apr 18 publish');
load('Exp001_Fluo_002_001_sequenceDataFiltered_bandpass0.5_12_box_gp_done.mat','rsDataT','cortexMask','vfsT');

rsData=rsDataT(:,:,9001:18000);
vfs=vfsT(:,:,9001:18000);
%%Fig1b 
figure 
[cfs,frequencies] = cwt(rsData(36,19,3250:4000),150./[0.1:0.01:10],'cmor1-1',1/150);
t = linspace(0,2,300);
imagesc(t,frequencies,abs(cfs(:,1:300))); 
grid off;
set(gca,'YDir','reverse');
xlabel('Time(s)'); 
ylabel('Frequency(Hz)');
shading interp
colormap(parula)
%colorbar('south')

%% surf of the amplitude,fig.1 d, fully awake
cortexMask_in=zeros(size(cortexMask,1),size(cortexMask,2));
index=0;
for r=1:size(cortexMask,1)
    for c=1:size(cortexMask,2)
        index=index+1;
        cortexMask_in(r,c)=index;                         
    end
end

cortex_index=reshape(cortexMask_in,[size(rsData,1)*size(rsData,2) 1]);

array2d=zeros(size(cortexMask,1)*size(cortexMask,2),size(rsData,3));
for k=1:size(rsData,3)
    array=reshape(rsData(:,:,k),[size(rsData,1)*size(rsData,2) 1]);
    array2d(:,k)=array;
end
figure
set(groot,'DefaultFigureColormap',jet)
imagesc(array2d(794:832,3250:3550))% 794-832 col: 19
colormap(darkb2r(-0.005,0.005))
shading interp
%% plot oscillation of three different pixels
%  global signal
ys=zeros(1,size(rsData,3)-1);
for i=1:size(rsData,3)
    ys(i)=sum(sum(rsData(:,:,i).*cortexMask))/length(find(cortexMask~=0));
end

figure
plot(ys(3250:3550));
legend('average voltage amplitude')
hold on
plot(squeeze(rsData(6,19,3250:3550)));
hold on
plot(squeeze(rsData(22,19,3250:3550)));
hold on
plot(squeeze(rsData(36,19,3250:3550)));
hold on
line([0,300],[0,0])
% hold on
% line([149,149],[-0.005,0.005],'color','k')
hold on
line([162,162],[-0.005,0.005],'color','k')
hold on
line([175,175],[-0.005,0.005],'color','k')
% hold on
% line([130,130],[-0.005,0.005],'color','k')

%% plot vfs in certain space 
interval=1;
%raw
rawmatrix=zeros(size(cortexMask,1),size(cortexMask,2));
for i=1:interval+1:size(cortexMask,1)
    rawmatrix(i,:)=1;
end
%figure;imagesc(rawmatrix)

%column
colmatrix=zeros(size(cortexMask,1),size(cortexMask,2));
for i=1:interval+1:size(cortexMask,2)
    colmatrix(:,i)=1;
end
%figure;imagesc(colmatrix)
intermatrix=rawmatrix & colmatrix;
%figure;imagesc(intermatrix)
%% plot phase velocity field, amplitude and mask
%%162
figure;
set (gcf,'Position',[200,100,420*0.65,340*0.65], 'color','w')
itime=3250+162;
useAmplitude=false;

signalGrid = rsData(:,:,itime);
sigLims = [-0.005 0.005];

% Apply mask to data
cortexMask(find(cortexMask<0.7))=NaN;
cortexMask(find(cortexMask>=0.7))=1;
signalGrid = signalGrid .* cortexMask;
% Plot signal grid
imagesc(signalGrid, sigLims)
colormap(darkb2r(-0.005,0.005))
h=imagesc(signalGrid, sigLims);
set(h,'alphadata',~isnan(signalGrid))

% Plot velocity field
hold on
vf = vfs(:,:,itime) * vfScale.*cortexMask;
ha2 = quiver(real(vf).*intermatrix, imag(vf).*intermatrix);
set(ha2,'color',[0 0 0],'ShowArrowHead','on','MaxHeadSize',0.8,'AutoScale','on', 'AutoScaleFactor', 2)
set(gca,'YDir','reverse');
%title(num2str(itime));
grid off
axis off
hold off

%% 175
figure;
set (gcf,'Position',[200,100,420*0.65,340*0.65], 'color','w')
itime=3250+175;
useAmplitude=false;

signalGrid = rsData(:,:,itime);
sigLims = [-0.005 0.005];

% Apply mask to data
cortexMask(find(cortexMask<0.7))=NaN;
cortexMask(find(cortexMask>=0.7))=1;
signalGrid = signalGrid .* cortexMask;
% Plot signal grid
imagesc(signalGrid, sigLims)
colormap(darkb2r(-0.005,0.005))
h=imagesc(signalGrid, sigLims);
set(h,'alphadata',~isnan(signalGrid))

% Plot velocity field
hold on
vf = vfs(:,:,itime) * vfScale.*cortexMask;
ha2 = quiver(real(vf).*intermatrix, imag(vf).*intermatrix);
set(ha2,'color',[0 0 0],'ShowArrowHead','on','MaxHeadSize',0.8,'AutoScale','on', 'AutoScaleFactor', 2)
set(gca,'YDir','reverse');
%title(num2str(itime));
grid off
axis off
hold off



