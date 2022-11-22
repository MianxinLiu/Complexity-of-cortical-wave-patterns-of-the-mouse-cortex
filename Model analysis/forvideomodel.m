fileName = 'p0Ve';
saveVideo = true;
% Video file location
vidLoc = '.';
% Video file name
vidName = fileName(1:end-4);
% Video frame rate
vidFrameRate = 20;
% Choose scale of velocity vector fields
vfScale = 1;

if saveVideo
disp('Writing video to file...'); tic
fig = figure('Name', 'fileName');
% Create video file
vidTitle = strcat(vidName, datestr(now,'ddmm'), '_', ...
datestr(now,'HHMM'), '.mp4');
vidObj = VideoWriter(vidTitle,'MPEG-4');
vidObj.FrameRate = vidFrameRate;
open(vidObj);
% Advance time and save video frames
start=150;
for itime = start:start+450
% Set signal and colors based on whether you are plotting amplitude
% or phase
signalGrid = Ve(:,:,itime);

% colorMapSpec = pmkmp_new;
colorMapSpec=jet;
sigLims = [-66,-50];
vfColor = [0 0 0];

% Plot signal grid
imagesc(signalGrid, sigLims)
% colormap(gca, colorMapSpec)

h=imagesc(signalGrid, sigLims);
set(h,'alphadata',~isnan(signalGrid))
% colormap(gca,darkb2r(-0.005,0.005))
colorbar
% Plot velocity field
% hold on
% vf = vfs(:,:,itime) * vfScale;
% quiver((1:size(signalGrid,2)), (1:size(signalGrid,1)), ...
% real(vf).*intermatrix, imag(vf).*intermatrix, 0, 'Color', vfColor)
% set(gca,'YDir','reverse');
% axis off
% hold off
% Update title with current time
title(sprintf('%0.2f s', (itime-start)/Fs))
% Write to video file
writeVideo(vidObj, im2frame(print(fig,'-RGBImage')));
end
% Close video file
close(vidObj);
toc
end

fileName = 'p0Qe';
saveVideo = true;
% Video file location
vidLoc = '.';
% Video file name
vidName = fileName(1:end-4);
% Video frame rate
vidFrameRate = 20;
% Choose scale of velocity vector fields
vfScale = 1;

if saveVideo
disp('Writing video to file...'); tic
fig = figure('Name', 'fileName');
% Create video file
vidTitle = strcat(vidName, datestr(now,'ddmm'), '_', ...
datestr(now,'HHMM'), '.mp4');
vidObj = VideoWriter(vidTitle,'MPEG-4');
vidObj.FrameRate = vidFrameRate;
open(vidObj);
% Advance time and save video frames
start=150;
for itime = start:start+450
% Set signal and colors based on whether you are plotting amplitude
% or phase
signalGrid = Qe(:,:,itime);

% colorMapSpec = pmkmp_new;
colorMapSpec=jet;
sigLims = [0.0002,0.03];
vfColor = [0 0 0];

% Plot signal grid
imagesc(signalGrid, sigLims)
% colormap(gca, colorMapSpec)

h=imagesc(signalGrid, sigLims);
set(h,'alphadata',~isnan(signalGrid))
% colormap(gca,darkb2r(-0.005,0.005))
colorbar
% Plot velocity field
% hold on
% vf = vfs(:,:,itime) * vfScale;
% quiver((1:size(signalGrid,2)), (1:size(signalGrid,1)), ...
% real(vf).*intermatrix, imag(vf).*intermatrix, 0, 'Color', vfColor)
% set(gca,'YDir','reverse');
% axis off
% hold off
% Update title with current time
title(sprintf('%0.2f s', (itime-start)/Fs))
% Write to video file
writeVideo(vidObj, im2frame(print(fig,'-RGBImage')));
end
% Close video file
close(vidObj);
toc
end