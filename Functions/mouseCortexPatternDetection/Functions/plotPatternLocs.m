function plotPatternLocs(allLocs, pattTypes, realTime, ntrials, nspacebins, cortexMask, smoothScale)
% Takes the patterns cell array ALLPATTERNLOCS from findAllPatterns.m and
% produces a number of plots

plotTime = true;
if isscalar(realTime)
    realTime = 1:realTime;
    plotTime = false;
end
if ~exist('ntrials', 'var')
    ntrials = 1;
end
if ~exist('smoothScale', 'var')
    smoothScale = 1;
end

% Combine patterns across all trials
npat = length(pattTypes);
combLocs = cell(1, npat);
noPatts = false(1, npat);
if ~iscell(allLocs{1})
    allLocs = {allLocs};
end

for ipatt = 1:npat
    % Reshape cell array to combine across all trials
    thisPatt = cellfun(@(x) x{ipatt}, allLocs, 'UniformOutput', 0);
    thisPatt = cat(1, thisPatt{:});
    if isempty(thisPatt)
        noPatts(ipatt) = true;
    else
        combLocs{ipatt} = thisPatt;
    end
end

combLocs = combLocs(~noPatts);
pattTypes = pattTypes(~noPatts);
npat = length(pattTypes);

% Plot distributions over time and space
for ipatt = 1:npat
    thisPatt = combLocs{ipatt};
    
    % Take histogram of pattern times
    nActive = histcounts(thisPatt(:,3), 1:length(realTime));
    nActive = nActive / ntrials;
    percActive = 100 * sum(nActive) / length(realTime);
    if plotTime
        subplot(2, npat, ipatt)
        plot(realTime(1:end-1), smooth(nActive, smoothScale))
        axis tight
        tightLims = ylim;
        ylim([tightLims(1)-0.1, tightLims(2)+0.1])
        xlabel('Time (s)')
        ylabel('# patterns active')
        title(sprintf('%s, %0.1f%% active', pattTypes{ipatt}, percActive))
        subplot(2, npat, npat+ipatt)
    else
        subplot(2, ceil(npat/2), ipatt)
    end
    
    % Take 2D histogram of pattern locations
    thisComplexLoc = thisPatt(:,1) + 1i*thisPatt(:,2);
    if all(abs(abs(thisComplexLoc(:))-1) < 10*eps)
        % Check if pattern is plane wave with direction instead of location
        % Kiki here I reverse Y-axis
        polarhistogram(2*pi-angle(thisComplexLoc), 24)
        % polarhistogram(angle(thisComplexLoc), 24)
        title('Plane wave propagation directions')
    else
        a=histogram2(thisPatt(:,1), thisPatt(:,2), nspacebins, 'FaceColor', 'flat', ...
            'ShowEmptyBins', 'on', 'DisplayStyle', 'tile');
        counts = a.Values.*cortexMask;
        imagesc(counts);
        axis equal
        set(gca, 'XTick', [], 'YTick', [])
        xlabel('X location')
        ylabel('Y location')
        title('Spatial centre locations')
        colorbar
    end
    
    if ~plotTime
        title(sprintf('%s, %0.1f%% active', pattTypes{ipatt}, percActive))
    end
    

    
end

