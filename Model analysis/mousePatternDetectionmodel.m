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
% Visualise when patterns are active over time
figure
realTime = (1:size(vfs,3))/Fs;
imagesc(realTime, 1:length(pattTypes), activeArray)
colormap(1-bone)
set(gca, 'YTick', 1:length(pattTypes), 'YTickLabel', pattTypes)
title('Pattern activity')
xlabel('Time (s)')

%% Display the spatial distributions of patterns
figure
plotPatternLocs(pattLocs, pattTypes, realTime, 1, size(vfs(:,:,1)), cortexMask)

%% Also plot the durations and mean displacement of each pattern type
figure
unqPatts = unique(patterns(:,1));
maxDur = max(patterns(:, strcmp(colNames, 'duration')))/Fs;
maxDisp = max(patterns(:, strcmp(colNames, 'meanDisplacement')))*Fs;
for ipatt = 1:length(unqPatts)
    thisPatt = patterns(patterns(:,1)==unqPatts(ipatt), :);
    thisType = pattTypes{unqPatts(ipatt)};
    
    % Histogram of durations
    subplot(2, length(unqPatts), ipatt)
    thisDurs = thisPatt(:, strcmp(colNames, 'duration'));
    histogram(thisDurs/Fs, linspace(0, maxDur, 10))
    xlabel('Duration (s)')
    ylabel('Counts')
    title(sprintf('%s, mean %0.3g', thisType, mean(thisDurs/Fs)))
    
    % Histogram of displacements
    if ~strcmp(thisType, 'planeWave') && ~strcmp(thisType, 'synchrony')
        subplot(2, length(unqPatts), length(unqPatts)+ipatt)
        thisDisp = thisPatt(:, strcmp(colNames, 'meanDisplacement'));
        histogram(thisDisp*Fs, linspace(0, maxDisp, 10))
        xlabel('Displacement (grid/s)')
        ylabel('Counts')
        title(sprintf('Mean %0.3g', mean(thisDisp*Fs)))
    end
end
