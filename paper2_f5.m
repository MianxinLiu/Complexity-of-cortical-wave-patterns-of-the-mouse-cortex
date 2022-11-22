clear;

explabel={[2,18],[2,7],[1,3],[1,3],[1,6]};
path={'2017 APR 19 publish','2018Feb 04 M2553M publish','2018Feb 04 M2555M publish','2018Feb 04 M2556M publish','2018Feb 04 M2560F publish'};

%% Set pattern detection parameters
Fs=150;
params = setPatternParams(Fs);
mid_ver=26;
mid_hor=22;%26;

pn=zeros(2,6,5);% pattern number (time*number)
pnt=zeros(2,6,5);% pattern number (time)
pnn=zeros(2,6,5);% pattern number (number)

ppn=zeros(2,5,5);% proportion of pattern number
params.minCritRadius=3;
params.minDuration=2;


for M=1:5
    for state=1:2
        exp=explabel{M}(state);
        cd(['/media/user/Elements/NC Data/' path{M}]);
        if exp<10
            load (['Exp001_Fluo_00' num2str(exp) '_001_sequenceDataFiltered_bandpass0.5_12_box_gp_done.mat'],'vfsT','cortexMask');
        else
             load (['Exp001_Fluo_0' num2str(exp) '_001_sequenceDataFiltered_bandpass0.5_12_box_gp_done.mat'],'vfsT','cortexMask');
        end
        vfs_relative = vfsT;
        vfs_relative = vfs_relative.*cortexMask;

        % %% Find all patterns
        %% left
        [patternsl, pattTypes, colNames, pattLocs] = ...
            findAllPatterns(double(real(vfs_relative(:,1:mid_ver,:))), double(imag(vfs_relative(:,1:mid_ver,:))), params);

        % Also make a binary array showing the patterns active in every time step
        activeArray_l = makeActivePatternsArray(patternsl, length(pattTypes), ...
            size(vfs_relative,3));

        %% right
        [patternsr, pattTypes, colNames, pattLocs] = ...
            findAllPatterns(double(real(vfs_relative(:,mid_ver:end,:))), double(imag(vfs_relative(:,mid_ver:end,:))), params);

        % Also make a binary array showing the patterns active in every time step
        activeArray_r = makeActivePatternsArray(patternsr, length(pattTypes), ...
            size(vfs_relative,3));

        %% 
        activeArray=activeArray_l+activeArray_r;
        patterns=[patternsl;patternsr];

        for i=1:5
            pnt(state,i,M)=pnt(1,i,M)+mean(patterns(find(patterns(:,1)==i),4));
            pnn(state,i,M)=pnn(1,i,M)+length(patterns(find(patterns(:,1)==i),4));
        end

        pn(state,1:5,M)=sum(activeArray');
        pn(state,6,M)=length(find(sum(activeArray)==0));
        pnx=tabulate(cumsum(single((sum(activeArray))>0)));
        b=find(pnx(:,2)>1);
        pnt(state,6,M)=mean(pnx(b,2));
        pnn(state,6,M)=length((b));
        for i=1:5
            if i<=4
            ppn(state,i,M)=length(find(sum(activeArray(3:5,:))==i-1));
            else
            ppn(state,i,M)=length(find(sum(activeArray(3:5,:))>3));
            end
        end 
    end
end
save fig5.mat pn pnt pnn ppn
 
%fully awake
path_fa={'2017 Apr 18 publish','2017Aug 05 M2242M publish','2017July24 M2259F publish','2017 July24 M2242M publish'};
explabel_fa={[1:10],[1],[5],[1:7]};

Fs=150;
params = setPatternParams(Fs);
mid_ver=26;
mid_hor=22;%26;

pn_fw=zeros(6,19);% pattern number (time*number)
pnt_fw=zeros(6,19);% pattern number (time)
pnn_fw=zeros(6,19);% pattern number (number)

ppn_fw=zeros(5,19);% proportion of pattern number
params.minCritRadius=3;
params.minDuration=2;

count=1;
for M=1:4
    n=length(explabel_fa{M});
    for state=1:n
        exp=explabel_fa{M}(state);
        cd(['/media/user/Elements/NC Data/' path_fa{M}]);
        if exp<10
            load (['Exp001_Fluo_00' num2str(exp) '_001_sequenceDataFiltered_bandpass0.5_12_box_gp_done.mat'],'vfsT','cortexMask');
        else
             load (['Exp001_Fluo_0' num2str(exp) '_001_sequenceDataFiltered_bandpass0.5_12_box_gp_done.mat'],'vfsT','cortexMask');
        end
        vfs_relative = vfsT;
        vfs_relative = vfs_relative.*cortexMask;

        % %% Find all patterns
        %% left
        [patternsl, pattTypes, colNames, pattLocs] = ...
            findAllPatterns(double(real(vfs_relative(:,1:mid_ver,:))), double(imag(vfs_relative(:,1:mid_ver,:))), params);

        % Also make a binary array showing the patterns active in every time step
        activeArray_l = makeActivePatternsArray(patternsl, length(pattTypes), ...
            size(vfs_relative,3));

        %% right
        [patternsr, pattTypes, colNames, pattLocs] = ...
            findAllPatterns(double(real(vfs_relative(:,mid_ver:end,:))), double(imag(vfs_relative(:,mid_ver:end,:))), params);

        % Also make a binary array showing the patterns active in every time step
        activeArray_r = makeActivePatternsArray(patternsr, length(pattTypes), ...
            size(vfs_relative,3));

        %% 
        activeArray=activeArray_l+activeArray_r;
        patterns=[patternsl;patternsr];

        for i=1:5
            pnt_fw(i,count)=pnt_fw(i,count)+mean(patterns(find(patterns(:,1)==i),4));
            pnn_fw(i,count)=pnn_fw(i,count)+length(patterns(find(patterns(:,1)==i),4));
        end

        pn_fw(1:5,count)=sum(activeArray');
        pn_fw(6,count)=length(find(sum(activeArray)==0));
        pnx_fw=tabulate(cumsum(single((sum(activeArray))>0)));
        b=find(pnx_fw(:,2)>1);
        pnt_fw(6,count)=mean(pnx_fw(b,2));
        pnn_fw(6,count)=length((b));
        for i=1:5
            if i<=4
            ppn_fw(i,count)=length(find(sum(activeArray(3:5,:))==i-1));
            else
            ppn_fw(i,count)=length(find(sum(activeArray(3:5,:))>3));
            end
        end 
        count=count+1;
    end
end

save fig5_2.mat pn_fw pnt_fw pnn_fw ppn_fw

% load intermedia data for reproduce the results in a quicker way (Deposited data links provided in the paper)
load('fig5.mat')

%% plot fig 5a
figure

pnt(isnan(pnt))=0;
pnM=pnt/150*1000;
pnt_fw(isnan(pnt_fw))=0;
pnM_fw=pnt_fw/150*1000;
for i=1:6
    boxplot(squeeze(pnM(1,i, :)), 'Positions', i-0.2,'Colors','b')
    hold on
    boxplot(squeeze(pnM(2,i, :)), 'Positions',i,'Colors','r')
    boxplot(squeeze(pnM_fw(i, :)), 'Positions',i+0.2,'Colors',[1,0.54,0])
end

for m=1:5
    scatter((1:6)-0.2+0.1,squeeze(pnM(1,:,m)),64,'b.')
end
for m=1:5
    scatter((1:6)+0.1,squeeze(pnM(2,:,m)),64,'r.')
end
for m=1:19
    scatter((1:6)+0.2+0.1,squeeze(pnM_fw(:,m)),64,[1,0.54,0],'.')
end

xlim([0,6.5])
ylim([-0.5,120])
xticks([1:6])
xticklabels({'plane','standing','sink','source','saddle','unclassified'})
ylabel('Average duration(ms)')


%% Fig5b,times of pattern number/sec
figure
pnM=pnn(:,:,:)/180;
pnM_fw=pnn_fw(:,:,:)/180;

for i=1:6
    boxplot(squeeze(pnM(1,i, :)), 'Positions', i-0.2,'Colors','b')
    hold on
    boxplot(squeeze(pnM(2,i, :)), 'Positions',i,'Colors','r')
    boxplot(squeeze(pnM_fw(i, :)), 'Positions',i+0.2,'Colors',[1,0.54,0])
end

for m=1:5
    scatter((1:6)-0.2+0.1,squeeze(pnM(1,:,m)),64,'b.')
end
for m=1:5
    scatter((1:6)+0.1,squeeze(pnM(2,:,m)),64,'r.')
end
for m=1:19
    scatter((1:6)+0.2+0.1,squeeze(pnM_fw(:,m)),64,[1,0.54,0],'.')
end

xlim([0,6.5])
ylim([-3,105])
xticks([1:6])
xticklabels({'plane','standing','sink','source','saddle','unclassified'})
ylabel('Pattern number/sec')

%% Fig5c,fit occupation time, wave number
clear pnnM pnnMr pnnM_fw pnnMr_fw
figure
pnnM=squeeze(nanmean(ppn(:,:,:),3))/(180*150);
pnnM(sum(pnnM')==0,:)=[];

pnnM_fw=squeeze(nanmean(ppn_fw(:,:),2))/(180*150);
pnnM_fw(sum(pnnM_fw')==0,:)=[];

for i=1:5
   pnnMr(:,i)=pnnM(:,5+1-i); 
end

for i=1:5
   pnnMr_fw(i)=pnnM_fw(5+1-i); 
end

bar([1, 1.5, 2],[pnnMr; pnnMr_fw],'stacked')
xticks([1,1.5,2])
xticklabels({'anesthetized','post woken','fully awake'})

legend('p>3','p3','p2','p1','p0')
ylim([0,1]);
ylabel('Proportion of pattern number/frame')