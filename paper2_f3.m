clear;
%% svd from all mice
%% using standard cortex mask

explabel={[2,18],[2,7],[1,3],[1,3],[1,6]};
path={'2017 APR 19 publish','2018Feb 04 M2553M publish','2018Feb 04 M2555M publish','2018Feb 04 M2556M publish','2018Feb 04 M2560F publish'};

vfs_anes=[];
for M=1:5%M=6 is M1 awake states
    exp=explabel{M}(1);
    cd(['/media/user/Elements/NC Data/' path{M}]);
    if exp<10
        load (['Exp001_Fluo_00' num2str(exp) '_001_sequenceDataFiltered_bandpass0.5_12_box_gp_done.mat'],'vfsT');
    else
        load (['Exp001_Fluo_0' num2str(exp) '_001_sequenceDataFiltered_bandpass0.5_12_box_gp_done.mat'],'vfsT');
    end
    vfs_anes=cat(3,vfs_anes,vfsT(3:40,1:48,5000:7000));
end


vfs_awak=[];
for M=1:5%M=6 is M1 awake states
    exp=explabel{M}(2);
    cd(['/media/user/Elements/NC Data/' path{M}]);
    if exp<10
        load (['Exp001_Fluo_00' num2str(exp) '_001_sequenceDataFiltered_bandpass0.5_12_box_gp_done.mat'],'vfsT');
    else
        load (['Exp001_Fluo_0' num2str(exp) '_001_sequenceDataFiltered_bandpass0.5_12_box_gp_done.mat'],'vfsT');
    end
    vfs_awak=cat(3,vfs_awak,vfsT(3:40,1:48,5000:7000));
end

path_fa={'2017 Apr 18 publish','2017Aug 05 M2242M publish','2017July24 M2259F publish','2017 July24 M2242M publish'};
explabel_fa={[1:10],[1],[5],[1:7]};

vfs_wake=[];
for M=1:4%M=6 is M1 awake states
    n=length(explabel_fa{M});
    for state=1:n
        exp=explabel_fa{M}(state);
        cd(['/media/user/Elements/NC Data/' path_fa{M}]);
        if exp<10
            load (['Exp001_Fluo_00' num2str(exp) '_001_sequenceDataFiltered_bandpass0.5_12_box_gp_done.mat'],'vfsT');
        else
            load (['Exp001_Fluo_0' num2str(exp) '_001_sequenceDataFiltered_bandpass0.5_12_box_gp_done.mat'],'vfsT');
        end
        if M==1
            vfs_wake=cat(3,vfs_wake,vfsT(3:40,1:48,5000:7000));
        else
            vfs_wake=cat(3,vfs_wake,vfsT(1:38,1:48,5000:7000));
        end
    end
end

vfsT=[];
vfsT=cat(3,vfs_anes,vfs_awak);
vfsT=cat(3,vfsT,vfs_wake);
T=size(vfsT,3);
cortexMaskS=zeros(38,48);
% cd( '\\158.182.15.58\test2\Kiki\2017 Apr 19\0.5-12Hz\a0.5b10');
cd('/media/user/Elements/NC Data/2017 Apr 18 publish' )
load (['Exp001_Fluo_00' num2str(2) '_001_sequenceDataFiltered_bandpass0.5_12_box_gp_done' '.mat'],'cortexMask');
cortexMaskS=cortexMaskS+cortexMask(3:40,1:48);
% cd( 'N:\Kiki\2017 July 24\M2242M');
cd('/media/user/Elements/NC Data/2017July24 M2259F publish' )
load (['Exp001_Fluo_00' num2str(5) '_001_sequenceDataFiltered_bandpass0.5_12_box_gp_done' '.mat'],'cortexMask');
cortexMaskS=cortexMaskS+cortexMask(1:38,1:48);
% cd( 'N:\Kiki\2017 Aug 5');

cd('/media/user/Elements/NC Data/2017Aug 05 M2242M publish' )
load (['Exp001_Fluo_00' num2str(1) '_001_sequenceDataFiltered_bandpass0.5_12_box_gp_done' '.mat'],'cortexMask');
cortexMaskS=cortexMaskS+cortexMask(1:38,3:50);

load ('/media/user/Elements/NC codes/cortexMask_five_mice.mat');
cortexMaskS=cortexMaskS+cortexMask(1:38,3:50);
cortexMaskS(find(cortexMaskS<4))=0;
cortexMaskS(find(cortexMaskS>=4))=1;
% cortexMaskS(find(cortexMaskS<3))=0;
% cortexMaskS(find(cortexMaskS>=3))=1;

interval=1;
%raw
rawmatrix=zeros(size(cortexMaskS,1),size(cortexMaskS,2));
for i=1:interval+1:size(cortexMaskS,1)
    rawmatrix(i,:)=1;
end
%figure;imagesc(rawmatrix)

%column
colmatrix=zeros(size(cortexMaskS,1),size(cortexMaskS,2));
for i=1:interval+1:size(cortexMaskS,2)
    colmatrix(:,i)=1;
end
%figure;imagesc(colmatrix)
intermatrix=rawmatrix & colmatrix;

%% fig3a
[U0, S0, V0, U, S, V, reUav] = plotcsvd(vfsT.*cortexMaskS, 5, 1:T, false, 1, intermatrix);
save (['exp' num2str(exp) '_svdMode_all_temp.mat'], 'U0','S0','V0','U','S','V','reUav');

%% fig3b,projection

%% svd from all mice
W_all_anes=[];
period=1000:2000;

for M=1:5%M=6 is M1 awake states
    exp=explabel{M}(1);
    cd(['/media/user/Elements/NC Data/' path{M}]);
    if exp<10
        load (['Exp001_Fluo_00' num2str(exp) '_001_sequenceDataFiltered_bandpass0.5_12_box_gp_done.mat'],'vfsT');
    else
        load (['Exp001_Fluo_0' num2str(exp) '_001_sequenceDataFiltered_bandpass0.5_12_box_gp_done.mat'],'vfsT');
    end
    [nr, nc, nt, ~] = size(vfsT(3:40,1:48,period).*cortexMaskS);
    svdmat = reshape(vfsT(3:40,1:48,period).*cortexMaskS, nr*nc, [])';

    % Separate x and y components into separate variables for SVD
    svdmat = cat(2, real(svdmat), imag(svdmat));

    W0=svdmat/V0()';
    W=W0.^2./repmat(sum(W0.^2,2),1,size(V0,1));
    W=mean(W,1);
    W_all_anes=[W_all_anes;W];
end
        

W_all_awak=[];
for M=1:5%M=6 is M1 awake states
    exp=explabel{M}(2);
    cd(['/media/user/Elements/NC Data/' path{M}]);
    if exp<10
        load (['Exp001_Fluo_00' num2str(exp) '_001_sequenceDataFiltered_bandpass0.5_12_box_gp_done.mat'],'vfsT');
    else
        load (['Exp001_Fluo_0' num2str(exp) '_001_sequenceDataFiltered_bandpass0.5_12_box_gp_done.mat'],'vfsT');
    end
    
    [nr, nc, nt, ~] = size(vfsT(3:40,1:48,period).*cortexMaskS);
    svdmat = reshape(vfsT(3:40,1:48,period).*cortexMaskS, nr*nc, [])';

    % Separate x and y components into separate variables for SVD
    svdmat = cat(2, real(svdmat), imag(svdmat));

    W0=svdmat/V0()';
    W=W0.^2./repmat(sum(W0.^2,2),1,size(V0,1));
    W=mean(W,1);
    W_all_awak=[W_all_awak;W];
end
        
W_all_wake=[];
for M=1:4
    n=length(explabel_fa{M});
    for state=1:n
        exp=explabel_fa{M}(state);
        cd(['/media/user/Elements/NC Data/' path_fa{M}]);
        if exp<10
            load (['Exp001_Fluo_00' num2str(exp) '_001_sequenceDataFiltered_bandpass0.5_12_box_gp_done.mat'],'vfsT');
        else
            load (['Exp001_Fluo_0' num2str(exp) '_001_sequenceDataFiltered_bandpass0.5_12_box_gp_done.mat'],'vfsT');
        end
        
        if M==1
            [nr, nc, nt, ~] = size(vfsT(3:40,1:48,period).*cortexMaskS);
            svdmat = reshape(vfsT(3:40,1:48,period).*cortexMaskS, nr*nc, [])';
        else
            [nr, nc, nt, ~] = size(vfsT(:,1:48,period).*cortexMaskS);
            svdmat = reshape(vfsT(:,1:48,period).*cortexMaskS, nr*nc, [])';
        end

        % Separate x and y components into separate variables for SVD
        svdmat = cat(2, real(svdmat), imag(svdmat));

        W0=svdmat/V0()';
        W=W0.^2./repmat(sum(W0.^2,2),1,size(V0,1));
        W=mean(W,1);
        W_all_wake=[W_all_wake;W];        
    end
end
        
figure
errorbar(1:15,mean(W_all_anes(:,1:15)),std(W_all_anes(:,1:15))/sqrt(5),'-o','LineWidth',1.5);
hold on
errorbar(1:15,mean(W_all_awak(:,1:15)),std(W_all_awak(:,1:15))/sqrt(5),'-o','LineWidth',1.5);
hold on
errorbar(1:15,mean(W_all_wake(:,1:15)),std(W_all_wake(:,1:15))/sqrt(19),'-o','LineWidth',1.5);
legend('Anesthetized','post woken','fully awake')
xlabel('SVD modes')
ylabel('Variance (%)')

x1=sum(W_all_anes(:,1:5),2);
x2=sum(W_all_awak(:,1:5),2);
x3=sum(W_all_wake(:,1:5),2);
p1 = ranksum(x1,x2);
p2 = ranksum(x2,x3);
p3 = ranksum(x1,x3);