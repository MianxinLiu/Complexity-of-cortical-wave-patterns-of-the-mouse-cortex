%% load allen atlas2015, count ROI number of the whole brain

% load atlas matrix
load('/media/user/Elements/NC codes/023_2015_ROIsABM.mat');
myROI=ROI;
RC = [2 17 90 91 151 152 172 173 174 176 177 195 196 197 ...
    198 199 200 209 224 297 298 358 359 379 380 381 ...
    383 384 402 403 404 405 406 407];

% Removing regions with no data.
n=0;
for i=[9,24,27,28,29,30,33,46,53,54]
    myROI(i-n,:,:)=[];
    Acromym(i-n)=[];
    Labels(i-n)=[];
    n=n+1;
end
% combine VISa and VISrl into PTLp
for i=[27,5]
    myROI(i,:,:)=myROI(i,:,:)+myROI(i+1,:,:);
    myROI(i+1,:,:)=[];
    Acromym(i)={'PTLp'};
    Acromym(i+1)=[];
end
boundary=squeeze(sum(myROI,1));
boundary(boundary>1)=NaN;

%% adjust ROI
mouseResize = 0.5;
keepRows = 20:107;
keepCols = 6:108;
% Crop out data outside selected area
if ~isempty(keepRows) && ~isscalar(keepRows)
myROI = myROI(:,keepRows,:);
end
if ~isempty(keepCols) && ~isscalar(keepCols)
myROI = myROI(:,:,keepCols);
end
firstFrame = imresize(squeeze(myROI(1,:,:)), mouseResize);
rsmyROI = zeros([size(myROI, 1),size(firstFrame)]);
for itime = 1:size(myROI,1)
rsmyROI(itime,:,:) = imresize(squeeze(myROI(itime,:,:)), mouseResize);
end


%% coherernce between anes and awake
%% regional voltage data 

load('fig7model.mat','rsmyROI')
explabel={[2,18],[2,7],[1,3],[1,3],[1,6]};
path={'2017 APR 19 publish','2018Feb 04 M2553M publish','2018Feb 04 M2555M publish','2018Feb 04 M2556M publish','2018Feb 04 M2560F publish'};


%% statistics of coherence
for M=1:5
    for state=1:2
        exp=explabel{M}(state);
        cd(['/media/user/Elements/NC Data/' path{M}]);
        if exp<10
            load (['Exp001_Fluo_00' num2str(exp) '_001_sequenceDataFiltered_bandpass0.5_12_box_gp_done.mat'],'rsDataT','cortexMask');
        else
            load (['Exp001_Fluo_0' num2str(exp) '_001_sequenceDataFiltered_bandpass0.5_12_box_gp_done.mat'],'rsDataT','cortexMask');
        end
    %         cd( 'N:\Kiki');
    %         load (['cortexMask_five_mice' '.mat']);

        rsDataT_reg=zeros(round(size(rsmyROI,1)/2),size(rsDataT,3));
        for i=1:round(size(rsmyROI,1)/2)
            for j=1:size(rsDataT,3)
            rsDataT_reg(i,j)=sum(sum(rsDataT(:,:,j).*squeeze(rsmyROI(i,:,:)),1),2);
            end
        end

        cxf_all=zeros(round(size(rsmyROI,1)/2),round(size(rsmyROI,1)/2),129);
        for i=1:round(size(rsmyROI,1)/2)
           for j=1:round(size(rsmyROI,1)/2)
           [cxy,f] = mscohere(rsDataT_reg(i,9001:18000),rsDataT_reg(j,9001:18000),60,40,[],150); 
           cxf_all(i,j,:)=cxy;
           end
        end
        filenm = ['exp_' num2str(exp) 'cxf_all6040'  '.mat' ];
        save(filenm, 'cxf_all',  '-mat');
        clear rsDataT_reg cxf_all
    end
end

reg(:,:,1)=[8,3;4,12;7,4];
reg(:,:,2)=[8,6;4,11;7,8];
for M=1:5%M=6 is M1 awake states
    for state=1:2
        exp=explabel{M}(state);
        cd(['/media/user/Elements/NC Data/' path{M}]);
        load (['exp_' num2str(exp) 'cxf_all6040'  '.mat']);
        for j=1:3
            for k=1:2
            reg_cxf(j,:,state,k,M)=squeeze(cxf_all(reg(j,1,k),reg(j,2,k),:));
            end
        end
    end
end

load('fig7bd.mat')
%% coherence on difference frequency band
% delta
cxf_delta_a=zeros(2,15);
cxf_delta_w=zeros(2,15);
n=0;
for M=1:5
    for j=1:3
        n=n+1;
        for k=1:2                  
            cxf_delta_a(k,n)=nanmean(squeeze(reg_cxf(j,2:8,1,k,M)));                      
            cxf_delta_w(k,n)=nanmean(squeeze(reg_cxf(j,2:8,2,k,M)));                         
        end  
    end
end
% meandelta=[nanmean(cxf_delta_w');nanmean(cxf_delta_a')];
% stddelta=[nanstd(cxf_delta_w',0,1);nanstd(cxf_delta_a',0,1)];
  
%% fig7b
% mean+std
M=1;
figure
for j=1
    for k=1:2
    subplot(1,2,3-k)%weak,strong
         
    plot(f,nanmean(squeeze(reg_cxf(j,:,1,k,:))'));
    hold on
    plot(f,nanmean(squeeze(reg_cxf(j,:,2,k,:))'));
    hold on

    legend('anes','awake', 'AutoUpdate','off');
    hold on
    shadedErrorBar(f,nanmean(squeeze(reg_cxf(j,:,1,k,:))'),nanstd(squeeze(reg_cxf(j,:,1,k,:))'./sqrt(5)),'lineprops','-b')
    hold on
    shadedErrorBar(f,nanmean(squeeze(reg_cxf(j,:,2,k,:))'),nanstd(squeeze(reg_cxf(j,:,2,k,:))'./sqrt(5)),'lineprops','-r')

    xlim([0 12])
    ylim([0 1]) 
    box off
    end  
end
%% fig7d,bar
trial=15;
% delta
cols=[0,114,189;
    217,83,23;
    237,177,32;
    126,47,142;
    119,172,48;
    77,190,238]/256;
% subplot(1,3,1)
cxf_delta_a=zeros(2,15);
cxf_delta_w=zeros(2,15);
n=0;
for M=1:5
    for j=1:3
        n=n+1;
        for k=1:2                  
            cxf_delta_a(k,n)=mean(squeeze(reg_cxf(j,2:8,1,k,M)));  % k=1 strong                    
            cxf_delta_w(k,n)=mean(squeeze(reg_cxf(j,2:8,2,k,M)));                      
          
        end  
    end
end

for i=1:2
    boxplot(squeeze(cxf_delta_a(3-i,:)), 'Positions', i-0.2,'Colors',cols(1,:),'Symbol','+')
    hold on
    boxplot(squeeze(cxf_delta_w(3-i,:)), 'Positions', i+0.2,'Colors',cols(2,:),'Symbol','+')
end
sigstar({[1,2]-0.2 [1,2]+0.2},[signrank(cxf_delta_a(1,:),cxf_delta_a(2,:)),signrank(cxf_delta_w(1,:),cxf_delta_w(2,:))])
xlim([0.5,2.5])
ylim([0,1])
xticks([1,2])
xticklabels({'weak','strong'});
ylabel('coherence');
xlabel('Strength of long-range connection');
title('0.5-4Hz') 
box off

% 4-8hz
subplot(1,3,2)
figure
cxf_delta_a=zeros(2,15);
cxf_delta_w=zeros(2,15);
n=0;
for M=1:5
    for j=1:3
        n=n+1;
        for k=1:2                  
            cxf_delta_a(k,n)=mean(squeeze(reg_cxf(j,8:16,1,k,M)));                      
            cxf_delta_w(k,n)=mean(squeeze(reg_cxf(j,8:16,2,k,M)));                      
   
        end  
    end
end
for i=1:2
    boxplot(squeeze(cxf_delta_a(3-i,:)), 'Positions', i-0.2,'Colors',cols(1,:),'Symbol','+')
    hold on
    boxplot(squeeze(cxf_delta_w(3-i,:)), 'Positions', i+0.2,'Colors',cols(2,:),'Symbol','+')
end
sigstar({[1,2]-0.2 [1,2]+0.2},[signrank(cxf_delta_a(1,:),cxf_delta_a(2,:)),signrank(cxf_delta_w(1,:),cxf_delta_w(2,:))])
xlim([0.5,2.5])
ylim([0,1])
xticks([1,2])
xticklabels({'weak','strong'});
ylabel('coherence');
xlabel('Strength of long-range connection');
title('4-8Hz')
 
box off
% 8-12Hz
subplot(1,3,3)
cxf_delta_a=zeros(2,15);
cxf_delta_w=zeros(2,15);
n=0;
for M=1:5
    for j=1:3
        n=n+1;
        for k=1:2                  
            cxf_delta_a(k,n)=mean(squeeze(reg_cxf(j,16:21,1,k,M)));                      
            cxf_delta_w(k,n)=mean(squeeze(reg_cxf(j,16:21,2,k,M)));                      
            xlim([0 12])
            ylim([0 0.9])    
        end  
    end
end
for i=1:2
    boxplot(squeeze(cxf_delta_a(3-i,:)), 'Positions', i-0.2,'Colors',cols(1,:),'Symbol','+')
    hold on
    boxplot(squeeze(cxf_delta_w(3-i,:)), 'Positions', i+0.2,'Colors',cols(2,:),'Symbol','+')
end
sigstar({[1,2]-0.2 [1,2]+0.2},[signrank(cxf_delta_a(1,:),cxf_delta_a(2,:)),signrank(cxf_delta_w(1,:),cxf_delta_w(2,:))])
xlim([0.5,2.5])
ylim([0,0.8])
xticks([1,2])
xticklabels({'weak','strong'});
ylabel('coherence');
xlabel('Strength of long-range connection');
title('8-12Hz')  
box off


