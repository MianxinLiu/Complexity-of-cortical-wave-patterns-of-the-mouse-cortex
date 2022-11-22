
%% large wave speed
explabel={[2,18],[2,7],[1,3],[1,3],[1,6]};
path={'2017 APR 19 publish','2018Feb 04 M2553M publish','2018Feb 04 M2555M publish','2018Feb 04 M2556M publish','2018Feb 04 M2560F publish'};

for M=1:5
    for state=1:2
        rsDatam=[];
        exp=explabel{M}(state);
        cd(['/media/user/Elements/NC Data/' path{M}]);
        if exp<10
            load (['Exp001_Fluo_00' num2str(exp) '_001_sequenceDataFiltered_bandpass0.5_12_box_gp_done.mat'],'rsDataT','cortexMask');
        else
            load (['Exp001_Fluo_0' num2str(exp) '_001_sequenceDataFiltered_bandpass0.5_12_box_gp_done.mat'],'rsDataT','cortexMask');
        end
       
        %  global signal
        ys=zeros(1,size(rsDataT,3)-1);
        for i=1:size(rsDataT,3)
            ys(i)=sum(sum(rsDataT(:,:,i).*cortexMask))/length(find(cortexMask~=0));
        end

        [pks,locs] = findpeaks(ys,'MinPeakHeight',0.001) ;
        lw=zeros(size(locs,2),2);
        for i=1:size(locs,2)
            for j=1:100
                if locs(i)-j>0
                    if ys(locs(i)-j)<0
                        lw(i,1)=locs(i)-j+1;
                        break
                    end
                end      
            end
            for j=1:100
                if locs(i)+j<size(ys,2)
                    if ys(locs(i)+j)<0
                        lw(i,2)=locs(i)+j-1;
                        break
                    end
                end        
            end
        end
        for ic=size(lw,1):-1:1
            if lw(ic,1)*lw(ic,2)==0
                lw(ic,:)=[];
            end
        end
        rsDataf={};
        rsDatag={};
        rsDatam1=[];
        for ii=1:size(lw,1)
            rsData=[];
            for jj=2:lw(ii,2)-lw(ii,1)+1
                rsData(:,:,jj)=imgaussfilt(squeeze(rsDataT(:,:,jj+lw(ii,1)-1)),2);
            end    
            max_idx=zeros(size(rsData,1),size(rsData,2));
            for iii=1:size(rsData,1)
                for jjj=1:size(rsData,2)
                    if cortexMask(iii,jjj)>0
                    [max_value(iii,jjj),max_idx(iii,jjj)]=max(squeeze(rsData(iii,jjj,:)));     
                    end
                end
            end

            rsDataf{ii}=max_idx/Fs;
            [fx,fy]=gradient(max_idx/Fs);
            rsDatag{ii}=abs(fx+i*fy);
            gradient_in=reshape(abs(fx+i*fy),[1 size(cortexMask,1)*size(cortexMask,2)]);
            rsDatam1(ii)=median(gradient_in(find(gradient_in>0)));
        end
        rsDatam=[rsDatam,rsDatam1];
        save(['exp' num2str(exp) 'speedlarge.mat'], 'rsDataf', 'rsDatag', 'rsDatam');
    end 
end

path_fa={'2017 Apr 18 publish','2017Aug 05 M2242M publish','2017July24 M2259F publish','2017 July24 M2242M publish'};
explabel_fa={[1:10],[1],[5],[1:7]};

for M=1:4
    n=length(explabel_fa{M});
    rsDatam=[];
    for state=1:n
        exp=explabel_fa{M}(state);
        cd(['/media/user/Elements/NC Data/' path_fa{M}]);
        if exp<10
            load (['Exp001_Fluo_00' num2str(exp) '_001_sequenceDataFiltered_bandpass0.5_12_box_gp_done.mat'],'rsDataT','cortexMask');
        else
            load (['Exp001_Fluo_0' num2str(exp) '_001_sequenceDataFiltered_bandpass0.5_12_box_gp_done.mat'],'rsDataT','cortexMask');
        end
       
        %  global signal
        ys=zeros(1,size(rsDataT,3)-1);
        for i=1:size(rsDataT,3)
            ys(i)=sum(sum(rsDataT(:,:,i).*cortexMask))/length(find(cortexMask~=0));
        end

        [pks,locs] = findpeaks(ys,'MinPeakHeight',0.001) ;
        lw=zeros(size(locs,2),2);
        for i=1:size(locs,2)
            for j=1:100
                if locs(i)-j>0
                    if ys(locs(i)-j)<0
                        lw(i,1)=locs(i)-j+1;
                        break
                    end
                end      
            end
            for j=1:100
                if locs(i)+j<size(ys,2)
                    if ys(locs(i)+j)<0
                        lw(i,2)=locs(i)+j-1;
                        break
                    end
                end        
            end
        end
        for ic=size(lw,1):-1:1
            if lw(ic,1)*lw(ic,2)==0
                lw(ic,:)=[];
            end
        end
        rsDataf={};
        rsDatag={};
        rsDatam1=[];
        for ii=1:size(lw,1)
            rsData=[];
            for jj=2:lw(ii,2)-lw(ii,1)+1
                rsData(:,:,jj)=imgaussfilt(squeeze(rsDataT(:,:,jj+lw(ii,1)-1)),2);
            end    
            max_idx=zeros(size(rsData,1),size(rsData,2));
            for iii=1:size(rsData,1)
                for jjj=1:size(rsData,2)
                    if cortexMask(iii,jjj)>0
                    [max_value(iii,jjj),max_idx(iii,jjj)]=max(squeeze(rsData(iii,jjj,:)));     
                    end
                end
            end

            rsDataf{ii}=max_idx/Fs;
            [fx,fy]=gradient(max_idx/Fs);
            rsDatag{ii}=abs(fx+i*fy);
            gradient_in=reshape(abs(fx+i*fy),[1 size(cortexMask,1)*size(cortexMask,2)]);
            rsDatam1(ii)=median(gradient_in(find(gradient_in>0)));
        end
        rsDatam=[rsDatam,rsDatam1];
    end 
    save speedlarge.mat rsDataf rsDatag rsDatam 
end

all_large_speed=zeros(5,2);
all_large_speed_std=zeros(5,2);

for M=1:5
    for state=1:2
        exp=explabel{M}(state);
        cd(['/media/user/Elements/NC Data/' path{M}]);
        load (['exp' num2str(exp) 'speedlarge.mat']);

        all_large_speed(M,state)=mean(rsDatam);
        all_large_speed_std(M,state)=std(17.4./rsDatam);
        
    end
end
    
% all_large_speed_fw=17.4./[all_large_speed_0, all_large_speed_2,all_large_speed_3, all_large_speed_1];
% all_large_speed_std_fw=[all_large_speed_std_0/sqrt(10), all_large_speed_std_2,all_large_speed_std_3, all_large_speed_std_1/sqrt(7)];
all_large_speed_fw=zeros(1,4);
all_large_speed_std_fw=zeros(1,4);

for M=1:4%M=6 is M1 awake states
    n=length(explabel_fa{M});
    for state=1:n
        exp=explabel_fa{M}(state);
        cd(['/media/user/Elements/NC Data/' path_fa{M}]);
        load (['speedlarge.mat']);
        
        all_large_speed_fw(M)=17.4/mean(rsDatam);
        all_large_speed_std_fw(M)=std(rsDatam)/sqrt(n);
    end
end

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
%    plot(displace(M)+(1:2),17.4./all_large_speed(M,1:2),'-o', 'Color', cols(M,:));
   hold on
   errorbar(displace(M)+(1:2),17.4./all_large_speed(M,1:2),all_large_speed_std(M,1:2), '-o', 'Color', cols(M,:))
end

for M=1:4
%    plot(displace(M)+2.5, all_large_speed_fw(M,1),'-o','Color', cols(M+5,:));
   hold on
   errorbar(2.5+displace(M),all_large_speed_fw(M),all_large_speed_std_fw(M), '-o', 'Color', cols(M+5,:))
end
% legend('Mouse 1','Mouse 2','Mouse 3','Mouse 4','Mouse 5')
xlabel('Cortical state');
ylabel('Average speed (mm/s)');
xticks([1,2,2.5])
xticklabels({'anesthetized','post woken','fully awake'})

% statistical tests
[p,h]= signrank((17.4./all_large_speed(:,1)),(17.4./all_large_speed(:,2)),'tail','left');
