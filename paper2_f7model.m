% run codes in "Model simulation programs" to get simulation data
% run codes in "model analysis" to obtain the data for figures

% load intermedia data for reproduce the results in a quicker way (Deposited data links provided in the paper)
load('fig7ac.mat')

%% Fig7a,plot example
% plot example mean+std
% variable f is given by 'mscohere' function
figure
subplot 121
plot(f,mean(cell2mat(cxy_anes(2,:))'),'linewidth',2);hold on;
hold on
plot(f,mean(cell2mat(cxy_awake(2,:))'),'linewidth',2);hold on;
legend('anes','awake', 'AutoUpdate','off');
hold on
shadedErrorBar(f,mean(cell2mat(cxy_anes(2,:))'),std(cell2mat(cxy_anes(2,:))')./sqrt(200),'lineprops','-b')
hold on
shadedErrorBar(f,mean(cell2mat(cxy_awake(2,:))'),std(cell2mat(cxy_awake(2,:))')./sqrt(200),'lineprops','-r')
xlim([0.5,12]);ylim([0,0.3]);
box off;
title('Weak long-range connection');

subplot 122
plot(f,mean(cell2mat(cxy_anes(6,:))'),'linewidth',2);hold on;
hold on
plot(f,mean(cell2mat(cxy_awake(6,:))'),'linewidth',2);hold on;
legend('anes','awake', 'AutoUpdate','off');
hold on
shadedErrorBar(f,mean(cell2mat(cxy_anes(6,:))'),std(cell2mat(cxy_anes(6,:))')./sqrt(200),'lineprops','-b')
hold on
shadedErrorBar(f,mean(cell2mat(cxy_awake(6,:))'),std(cell2mat(cxy_awake(6,:))')./sqrt(200),'lineprops','-r')
xlim([0.5,12]);ylim([0,0.3]);
box off;
title('Strong long-range connection');

%% Fig7c,[0,5, 10],[20,25,30]

trial=100;
cols=[0,114,189;
    217,83,23;
    237,177,32;
    126,47,142;
    119,172,48;
    77,190,238]/256;
for m=2:4
    anes=coh_anes{m};
    anes_c1=anes(1:3,:);
    anes_c1=reshape(anes_c1,[1,3*size(anes,2)]);
    anes_c2=anes(5:7,:);
    anes_c2=reshape(anes_c1,[1,3*size(anes,2)]);
    anes_c=[anes_c1;anes_c2];
    awak=coh_awake{m};
    awak_c1=awak(1:3,:);
    awak_c1=reshape(awak_c1,[1,3*size(awak,2)]);
    awak_c2=awak(5:7,:);
    awak_c2=reshape(awak_c2,[1,3*size(awak,2)]);
    awak_c=[awak_c1;awak_c2];
    figure;
    for i=1:2
        boxplot(squeeze(anes_c(i,:)), 'Positions', i-0.2,'Colors',cols(1,:),'Symbol','+')
        hold on
        boxplot(squeeze(awak_c(i,:)), 'Positions', i+0.2,'Colors',cols(2,:),'Symbol','+')
    end
    sigstar({[1,2]-0.2 [1,2]+0.2},[signrank(anes_c(1,:),anes_c(2,:)),signrank(awak_c(1,:),awak_c(2,:))])
    xlim([0.5,2.5])
    if m==2
        ylim([0,0.3])
    else
        ylim([0,0.065])
    end
    xticks([1,2])
    xticklabels({'weak','strong'});
    ylabel('coherence');
    xlabel('Strength of long-range connection');
    switch m
        case 2
            title('0.5-4Hz') 
        case 3
            title('4-8Hz') 
        case 4
            title('8-12Hz') 
    end
    box off
end
