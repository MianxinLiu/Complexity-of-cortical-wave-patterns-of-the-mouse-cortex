% Script to find all patterns in a sequence of velocity vector fields and
% display some of their properties
%% Fig 4a

% plot vfs in certain space 
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

cd('/media/user/Elements/NC Data/2017 APR 19 publish');
load('Exp001_Fluo_002_001_sequenceDataFiltered_bandpass0.5_12_box_gp_done.mat','vfsT','cortexMask','rsDataT');

% plane wave
figure;
set (gcf,'Position',[200,100,420*0.8,340*0.8], 'color','w')
itime=9000+7030;
useAmplitude=false;

signalGrid = rsDataT(:,:,itime);
colorMapSpec = jet;
sigLims = [-0.005 0.005];
vfColor = [0 0 0];
vfScale =1;

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
vf = vfsT(:,:,itime) * vfScale.*cortexMask;
ha2 = quiver(real(vf).*intermatrix, imag(vf).*intermatrix);
set(ha2,'color',[0 0 0],'ShowArrowHead','on','MaxHeadSize',0.8,'AutoScale','on', 'AutoScaleFactor', 2)
set(gca,'YDir','reverse');
%title(num2str(itime));
grid off
axis off
title('Plane wave','FontWeight','normal')

% standing wave

figure;
set (gcf,'Position',[200,100,420*0.8,340*0.8], 'color','w')
itime=15328;
useAmplitude=false;

signalGrid = rsDataT(:,:,itime);
colorMapSpec = jet;
sigLims = [-0.005 0.005];
vfColor = [0 0 0];

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
vf = vfsT(:,:,itime) * vfScale.*cortexMask;
ha2 = quiver(real(vf).*intermatrix, imag(vf).*intermatrix);
set(ha2,'color',[0 0 0],'ShowArrowHead','on','MaxHeadSize',0.8,'AutoScale','on', 'AutoScaleFactor', 2)
set(gca,'YDir','reverse');
%title(num2str(itime));
grid off
axis off
title('Standing wave','FontWeight','normal')

% Source, sink and saddle

figure;
set (gcf,'Position',[200,100,420*0.8,340*0.8], 'color','w')
itime=22780;
useAmplitude=false;

signalGrid = rsDataT(:,:,itime);
colorMapSpec = jet;
sigLims = [-0.005 0.005];
vfColor = [0 0 0];

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
vf = vfsT(:,:,itime) * vfScale.*cortexMask;
ha2 = quiver(real(vf).*intermatrix, imag(vf).*intermatrix);
set(ha2,'color',[0 0 0],'ShowArrowHead','on','MaxHeadSize',0.8,'AutoScale','on', 'AutoScaleFactor', 2)
set(gca,'YDir','reverse');
%title(num2str(itime));
grid off
axis off
title('Source, sink and saddle','FontWeight','normal')
scatter(11,30,100,'k+','LineWidth',2.5) %source
scatter(37,41,100,'ko','LineWidth',2.5) %sink
scatter(37,31,100,'kx','LineWidth',2.5) %saddle

% Unclassified

figure;
set (gcf,'Position',[200,100,420*0.8,340*0.8], 'color','w')
itime=25188;
useAmplitude=false;

signalGrid = rsDataT(:,:,itime);
colorMapSpec = jet;
sigLims = [-0.005 0.005];
vfColor = [0 0 0];

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
vf = vfsT(:,:,itime) * vfScale.*cortexMask;
ha2 = quiver(real(vf).*intermatrix, imag(vf).*intermatrix);
set(ha2,'color',[0 0 0],'ShowArrowHead','on','MaxHeadSize',0.8,'AutoScale','on', 'AutoScaleFactor', 2)
set(gca,'YDir','reverse');
%title(num2str(itime));
grid off
axis off
title('Unclassified','FontWeight','normal')


%% Fig 4 b anethestized
cd('/media/user/Elements/NC Data/2017 APR 19 publish');
load('Exp001_Fluo_002_001_sequenceDataFiltered_bandpass0.5_12_box_gp_done.mat','vfsT','cortexMask','rsDataT');

Fs=150;
params = setPatternParams(Fs);
mid_ver=26;
plane_num=zeros(1,4);
mid_hor=22;%26;
activeArray_lo=zeros(4,size(vfsT,3));
vfs_relative = vfsT;
vfs_relative = vfs_relative.*cortexMask;

% %% Find all patterns
%% left
[patterns, pattTypes, colNames, pattLocs] = ...
    findAllPatterns(double(real(vfs_relative(:,1:mid_ver,:))), double(imag(vfs_relative(:,1:mid_ver,:))), params);

% Also make a binary array showing the patterns active in every time step
activeArray_l = makeActivePatternsArray(patterns, length(pattTypes), ...
    size(vfs_relative,3));

%% right
[patterns, pattTypes, colNames, pattLocs] = ...
    findAllPatterns(double(real(vfs_relative(:,mid_ver:end,:))), double(imag(vfs_relative(:,mid_ver:end,:))), params);

% Also make a binary array showing the patterns active in every time step
activeArray_r = makeActivePatternsArray(patterns, length(pattTypes), ...
    size(vfs_relative,3));

%% 
activeArray=activeArray_l+activeArray_r;


%% plot the timecorse of the wave types
activewave=zeros(6,size(activeArray,2));
activewave(1:5,:)=activeArray;
activewave(6,find(sum(activeArray)==0))=1;
% Visualise when patterns are active over time
figure
wavetype={'Plane wave','Standing wave','Source','Sink','Saddle','Unclassified'};
realTime = (1:150)/Fs;
imagesc(realTime, [1 2 3 4 5 6], activewave(:,3700+2950:3850+2950), [0,8])
colormap(1-bone)
set(gca, 'YTick', 1:length(wavetype), 'YTickLabel', wavetype)
title('Pattern activity')
xlabel('Time (s)')
title('Anesthestized')

%% Fig 4 b post woken
cd('/media/user/Elements/NC Data/2017 APR 19 publish');
load('Exp001_Fluo_018_001_sequenceDataFiltered_bandpass0.5_12_box_gp_done.mat','vfsT','cortexMask','rsDataT');

Fs=150;
params = setPatternParams(Fs);
mid_ver=26;
plane_num=zeros(1,4);
mid_hor=22;%26;
activeArray_lo=zeros(4,size(vfsT,3));
vfs_relative = vfsT;
vfs_relative = vfs_relative.*cortexMask;

% %% Find all patterns
%% left
[patterns, pattTypes, colNames, pattLocs] = ...
    findAllPatterns(double(real(vfs_relative(:,1:mid_ver,:))), double(imag(vfs_relative(:,1:mid_ver,:))), params);

% Also make a binary array showing the patterns active in every time step
activeArray_l = makeActivePatternsArray(patterns, length(pattTypes), ...
    size(vfs_relative,3));

%% right
[patterns, pattTypes, colNames, pattLocs] = ...
    findAllPatterns(double(real(vfs_relative(:,mid_ver:end,:))), double(imag(vfs_relative(:,mid_ver:end,:))), params);

% Also make a binary array showing the patterns active in every time step
activeArray_r = makeActivePatternsArray(patterns, length(pattTypes), ...
    size(vfs_relative,3));

%% 
activeArray=activeArray_l+activeArray_r;


%% plot the timecorse of the wave types
activewave=zeros(6,size(activeArray,2));
activewave(1:5,:)=activeArray;
activewave(6,find(sum(activeArray)==0))=1;
% Visualise when patterns are active over time
figure
wavetype={'Plane wave','Standing wave','Source','Sink','Saddle','Unclassified'};
realTime = (1:150)/Fs;
imagesc(realTime, [1 2 3 4 5 6], activewave(:,3700+2950:3850+2950), [0,8])
colormap(1-bone)
set(gca, 'YTick', 1:length(wavetype), 'YTickLabel', wavetype)
title('Pattern activity')
xlabel('Time (s)')
title('Post woken')

%% Fig 4 b fully awake
cd('/media/user/Elements/NC Data/2017 Apr 18 publish');
load('Exp001_Fluo_001_001_sequenceDataFiltered_bandpass0.5_12_box_gp_done.mat','rsDataT','cortexMask','vfsT');

Fs=150;
params = setPatternParams(Fs);
mid_ver=26;
plane_num=zeros(1,4);
mid_hor=22;%26;
activeArray_lo=zeros(4,size(vfsT,3));
vfs_relative = vfsT;
vfs_relative = vfs_relative.*cortexMask;

% %% Find all patterns
%% left
[patterns, pattTypes, colNames, pattLocs] = ...
    findAllPatterns(double(real(vfs_relative(:,1:mid_ver,:))), double(imag(vfs_relative(:,1:mid_ver,:))), params);

% Also make a binary array showing the patterns active in every time step
activeArray_l = makeActivePatternsArray(patterns, length(pattTypes), ...
    size(vfs_relative,3));

%% right
[patterns, pattTypes, colNames, pattLocs] = ...
    findAllPatterns(double(real(vfs_relative(:,mid_ver:end,:))), double(imag(vfs_relative(:,mid_ver:end,:))), params);

% Also make a binary array showing the patterns active in every time step
activeArray_r = makeActivePatternsArray(patterns, length(pattTypes), ...
    size(vfs_relative,3));

%% 
activeArray=activeArray_l+activeArray_r;


%% plot the timecorse of the wave types
activewave=zeros(6,size(activeArray,2));
activewave(1:5,:)=activeArray;
activewave(6,find(sum(activeArray)==0))=1;
% Visualise when patterns are active over time
figure
wavetype={'Plane wave','Standing wave','Source','Sink','Saddle','Unclassified'};
realTime = (1:150)/Fs;
imagesc(realTime, [1 2 3 4 5 6], activewave(:,3700+2950:3850+2950), [0,8])
colormap(1-bone)
set(gca, 'YTick', 1:length(wavetype), 'YTickLabel', wavetype)
title('Pattern activity')
xlabel('Time (s)')
title('Fully awake')


