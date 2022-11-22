function y=view_spatialtemporal_evolve(Ve)
% look at the spatial evolution of the simulated data 
% Ve is a N x N x L array to record the simulated data in L steps

figure;

set(gcf,'Position',[100,50,450,400]);
L=size(Ve,3);
% wrap = vertcat(jet,flipud(jet));
% Ve=Ve(:,2:L-1);
for step=1:L
    pcolor(Ve(:,:,step));  shading interp; %
    colorbar; 
    colormap jet;
%     colormap(wrap);
%     caxis([0,2*pi]);
    caxis([-65,-52]);
    set(gca,'YDir','reverse'); 
%     caxis([0,2]);
%     xlim([1,40]);ylim([1,40]);
%     hold on;
%     ha2 =quiver(Vy(:,:,step),Vx(:,:,step));
%     set(ha2,'color',[0 0 0],'ShowArrowHead','on','MaxHeadSize',2,'AutoScale','on', 'AutoScaleFactor', 4);    
    drawnow;  
    frame=getframe(gcf);  
    im=frame2im(frame); % make gif
    [I,map]=rgb2ind(im,20);          
    if step==1
        imwrite(I,map,'movie_awak.gif','gif', 'Loopcount',inf,'DelayTime',0.1);% create
    else
        imwrite(I,map,'movie_awak.gif','gif','WriteMode','append','DelayTime',0.1);
    end
%     
%     pause(0.1);
end