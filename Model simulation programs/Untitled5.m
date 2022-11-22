

p = 0.5:-0.1:0;
for id=1:length(p)
    [k,alpha,f]=stability_anesthesia_fieldmodel(p(id),'bottom');
    [TheoryWaveVol(id),TheoryWaveVol_largeWave(id)]=speed_by_LR(k,alpha,f);
end

% figure;
% plot(p,TheoryWaveVol,'o-');hold on;
% plot(p,TheoryWaveVol_largeWave,'o-');hold on;
% set(gca,'xdir','reverse');
% ylabel('10mm/s');

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