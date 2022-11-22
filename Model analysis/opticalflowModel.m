tr=4;
for p=0:5
fileLoc = '\\158.182.15.58\test2\Junhao';
% % File name
fileName = ['result_of_p' num2str(p) '_trial' num2str(tr) ];
cd(fileLoc)
load(fileName)
% Sampling frequency (in Hz)
Fs = 100;
% Spatial selection of recording to keep (set to empty or zero to use full
% recording)
keepRows = 1:50;
keepCols = 1:50;
time=1:3000;
cortexMask=ones(length(keepRows),length(keepCols),length(time));
%% Optical flow parameters
% Smoothing parameter: higher values will give a smoother velocity field
% (typically 0<OPALPHA<2).
opAlpha = 0.5;%0.5;
% Non-linearity penalty parameter: close to zero will be highly non-linear,
% large values will be approximately linear, resulting in faster
% computations but possibly less accurate flow fields
opBeta = 1;
% Use flag to calculate amplitude velocity fields rather than phase
useAmplitude = false;
badChannels=[];
% %% Load mouse data file
%
% fprintf('Loading file %s...\n', fileName); tic
% load(fullfile(fileLoc, fileName))

% Crop out data outside selected area
if ~isempty(keepRows) && ~isscalar(keepRows)
Ve = VolE(keepRows,:,time);
end
if ~isempty(keepCols) && ~isscalar(keepCols)
Ve = Ve(:,keepCols,time);
end

wvcfs=filterphase(Ve,Fs);
disp('Computing optical flow fields...'); tic
% Calculate velocity fields
[vx, vy, csteps] = opticalFlow(wvcfs, badChannels, opAlpha, opBeta, ...
~useAmplitude);
vfs = vx + 1i*vy;
vfs = vfs / mean(abs(vfs(:)));
toc
% Script to find all patterns in a sequence of velocity vector fields and
% display some of their properties

%% Set pattern detection parameters
params = setPatternParams(Fs);

% %% Find all patterns
% [patterns, pattTypes, colNames, pattLocs] = ...
%     findAllPatterns(real(vfs), imag(vfs), params);
% change from single to double
[patterns, pattTypes, colNames, pattLocs] = ...
    findAllPatterns(double(real(vfs)), double(imag(vfs)), params);

%% Also make a binary array showing the patterns active in every time step
activeArray = makeActivePatternsArray(patterns, length(pattTypes), ...
    size(vfs,3));
save (['result_of_p' num2str(p) '_trial' num2str(tr) '_50_50_3000.mat'])
end