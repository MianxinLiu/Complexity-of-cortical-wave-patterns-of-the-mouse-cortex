function [vx,vy,vfs,wvcfs] = optical_flow_model_gp_surrogate(data,fmin,fmax)


% Sampling frequency (in Hz)
Fs = 100;
% Scale by which to downsample spatially
mouseResize = 1;
% Flag to use mask located in meta.traceMask
useMask = 0;
% Spatial selection of recording to keep (set to empty or zero to use full
% recording)
%used in v4
keepRows = 1:50;
keepCols = 1:50;
%% Filtering parameters
% Centre frequency of filter (in Hz)
cfreq = 3;
% Morlet wavelet scale parameter: higher values have better frequency
% resolution but lower temporal resolution (around 5-8 is typical)
morletParam = 5;
%% Optical flow parameters
% Smoothing parameter: higher values will give a smoother velocity field
% (typically 0<OPALPHA<2).
opAlpha = 0.5;%0.5;
% Non-linearity penalty parameter: close to zero will be highly non-linear,
% large values will be approximately linear, resulting in faster
% computations but possibly less accurate flow fields
opBeta = 10;
% Use flag to calculate amplitude velocity fields rather than phase
useAmplitude = false;
%% Video output parameters
% Flag to save video of signal and velocity field
saveVideo = true;
% % Video file location
% vidLoc = '.';
% % Video file name
% vidName = fileName(1:end-4);
% % Video frame rate
% vidFrameRate = 20;
% Choose scale of velocity vector fields
vfScale = 1;
% %% Load mouse data file
%
% fprintf('Loading file %s...\n', fileName); tic
% load(fullfile(fileLoc, fileName))
% Choose mask to use
if useMask
cortexMask = mask; % Use provided mask for data
else
cortexMask = true(size(data,1), size(data,2)); % Use no mask
end
% Crop out data outside selected area
if ~isempty(keepRows) && ~isscalar(keepRows)
data = data(keepRows,:,:);
cortexMask = cortexMask(keepRows,:,:);
end
if ~isempty(keepCols) && ~isscalar(keepCols)
data = data(:,keepCols,:);
cortexMask = cortexMask(:,keepCols,:);
end

general.fmin = fmin; %Hz 0.2;%
general.fmax = fmax;%5; %Hz
general.attenuation = 20;%10;%dB
general.order = 3; %4 did not work
general.srate = 100; %sampling rate [Hz]
[b, a] = cheby2(general.order,general.attenuation,[general.fmin*2/general.srate general.fmax*2/general.srate]);%Riedner 2007

for j=1:size(data,2)
    for i=1:size(data,1)
        cache=filtfilt(b,a,squeeze(squeeze(double(data(i,j,:)))));
        mcache = mean(cache);
        rsData(i,j,:) = cache - mcache;
    end
end
clearvars data
% Find any recording channels that contain NaNs or have no signal: these
% can be smoothed over when calculating optical flow
nanChans = any(isnan(rsData),3);
zeroChans = all(rsData==0, 3);
badChannels = find(nanChans | zeroChans);
disp('Computing generalized phase'); tic
wvcfs=zeros(size(rsData,1),size(rsData,2),size(rsData,3));
lp=0.5;
for i=1:size(rsData,1)
    for j=1:size(rsData,2)
        if cortexMask(i,j)==1
        xf= squeeze(rsData(i,j,1:size(rsData,3)));
        [xgp,wt] = generalized_phase_vector(xf, Fs, lp );
        wvcfs(i,j,1:size(rsData,3))=xgp;
        end
    end
end

toc
disp('Computing optical flow fields...'); tic
% Calculate velocity fields
[vx, vy, csteps] = opticalFlowO(wvcfs, badChannels, opAlpha, opBeta, ...
~useAmplitude);
vfs = vx + 1i*vy;
vfs = vfs / mean(abs(vfs(:)));
% Apply mask to remove vectors outside the cortex
vfs = vfs .* repmat(cortexMask, 1, 1, size(vfs,3));
toc
% %% Make video of raw data and optical flow fields
% if saveVideo
%     disp('Writing video to file...'); tic
%     fig = figure('Name', 'fileName');
%     % Create video file
%     vidTitle = strcat(vidName, datestr(now,'ddmm'), '_', ...
%     datestr(now,'HHMM'), '.avi');
%     vidObj = VideoWriter(vidTitle);
%     vidObj.FrameRate = vidFrameRate;
%     open(vidObj);
%     % Advance time and save video frames
%     for itime = 1:500
%     % Set signal and colors based on whether you are plotting amplitude
%     % or phase
%     if useAmplitude
%     signalGrid = abs(wvcfs(:,:,itime));
%     colorMapSpec = parula;
%     sigLims = [min(abs(wvcfs(:))), max(abs(wvcfs(:)))];
%     vfColor = [1 1 1];
%     else
%     signalGrid = angle(wvcfs(:,:,itime));
%     colorMapSpec = pmkmp_new;
%     sigLims = [-pi pi];
%     vfColor = [0 0 0];
%     end
%     % Apply mask to data
%     signalGrid = signalGrid .* cortexMask;
%     % Plot signal grid
%     imagesc(signalGrid, sigLims)
%     colormap(gca, colorMapSpec)
%     % Plot velocity field
%     hold on
%     vf = vfs(:,:,itime) * vfScale;
%     quiver((1:size(signalGrid,2)), (1:size(signalGrid,1)), ...
%     real(vf), imag(vf), 0, 'Color', vfColor)
%     set(gca,'YDir','reverse');
%     hold off
%     % Update title with current time
%     title(sprintf('%0.2f s', itime/Fs))
%     % Write to video file
%     writeVideo(vidObj, im2frame(print(fig,'-RGBImage')));
%     end
%     % Close video file
%     close(vidObj);
%     toc
% end

