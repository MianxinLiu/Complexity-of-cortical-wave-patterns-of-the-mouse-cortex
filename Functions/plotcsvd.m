function [U0, S0, V0, U, S, V, reUav] = plotcsvd(vf, nmodeplot, realTime, useComplexSVD, vectorScale, intermx)
% PLOTCSVD performs a singular vector decomposition of a vector field.
%   Plots most dominant spatial modes and their time courses, and the
%   trajectory of the top 3 most dominant modes.
%
% INPUTS:
% VF, a complex velocity field., Must be a 3D or 4D vector containing
%   [row,column,time,trial].
% NMODESPLOT, the number of spatial modes to plot (default 5).
% REALTIME, vector containing the real time in seconds of all time points
%   in VF. If REALTIME is empty or zero, only spatial modes will be 
%   plotted (default 1:NT, i.e. sampling frequency 1).
% USECOMPLEXSVD, flag to decompose vectors as complex values instead of
%   separating x- and y-coordinates (default false).
% VECTORSCALE, number to increase or decrease vector lengths in plots
%   (default 1, corresponding to MATLAB default size).
%
% OUTPUTS:
% U0,S0,VO original dimention  
% U,S,V dimensionality reduction results
% Rory Townsend, Oct 2017
% rory.townsend@sydney.edu.au

%% Set default values
if ~exist('nmodeplot', 'var')
    nmodeplot = 5;
end
if ~exist('realTime', 'var')
    realTime = 1:size(vf,3);
end
if ~exist('useComplexSVD', 'var')
    useComplexSVD = false;
end
if ~exist('vectorScale', 'var')
    vectorScale = 1;
end

% Check whether to plot time courses or not
if isempty(realTime) || (~isvector(realTime) && realTime == 0)
    onlySpatial = true;
else
    onlySpatial = false;
end

[nr, nc, nt, ~] = size(vf);


%% Perform SVD with complex velocity field
svdmat = reshape(vf, nr*nc, [])';
    
if ~useComplexSVD
    % Separate x and y components into separate variables for SVD
    svdmat = cat(2, real(svdmat), imag(svdmat));
end

[U0,S0,V0] = svd(svdmat, 0);

V=V0;
if ~useComplexSVD
    V = V(1:nr*nc,:) - 1i*V(nr*nc+(1:nr*nc),:);
end

% Calculate total energy then crop to specified number of modes
eigVals = diag(S0);
totVar = sum(eigVals.^2);
prctVar = 100 * eigVals.^2 / totVar;
U = U0(:,1:nmodeplot);
S = S0(1:nmodeplot, 1:nmodeplot);
V = V(:,1:nmodeplot);

% Make time modes as positive as possible
isNegative = mean(real(U), 1) < 0;
U(:, isNegative) = -U(:, isNegative);
V(:, isNegative) = -V(:, isNegative);

% Find average trajectory across all trials
reUav = zeros(nt, nmodeplot);
imUav = reUav;
absUav = imUav;
for it = 1:nt
    Uav = U(it + (0:nt:size(U,1)-nt), :);
    reUav(it, :) = mean(real(Uav),1);
    imUav(it, :) = mean(imag(Uav),1);
    absUav(it, :) = mean(abs(Uav),1);
end
if useComplexSVD
    UavLims = [min([reUav(:); imUav(:); absUav(:)]), ...
        max([reUav(:); imUav(:); absUav(:)])];
else
    UavLims = [min(reUav(:)), max(reUav(:))];
end

%% Plot spatial modes containing most energy
set (gcf,'Position',[400,100,900,150], 'color','w')
% [ha, pos] = tight_subplot(1,nmodeplot,[.000001 .000001],[.001 .001],[.001 .001]);

for imode=1:nmodeplot
    % Plot spatial mode
%     if onlySpatial
%         subplot(1, ceil(nmodeplot/2), imode)
%     else
%         subplot(1,nmodeplot,imode)
%     end
    figure
    set (gcf,'Position',[200,100,420*0.65,340*0.65], 'color','w')
    thisMode = reshape(V(:,imode), nr, nc);
    quiver(real(thisMode).*intermx, imag(thisMode).*intermx, vectorScale)
    ha2 = quiver(real(thisMode).*intermx, imag(thisMode).*intermx, vectorScale);
    set(ha2,'color',[0 0 0],'ShowArrowHead','on','MaxHeadSize',0.8,'AutoScale','on', 'AutoScaleFactor', 2)
    set(gca,'YDir','reverse');
    grid off
    
    axis(0.5+[0 nc 0 nr])
    axis off
    title(sprintf('Var = %0.1f%%', prctVar(imode)))
    %title(sprintf('Var = %0.1f%%, %0.2g%%', prctVar(imode), ...
    %    sum(prctVar(1:imode))))
    
end
% suptitle('Top SVD modes')

%% OPTIONAL: Also plot trajectory of top 3 modes
% figure
% plot3(smooth(reUav(:,1)), smooth(reUav(:,2)), smooth(reUav(:,3)))
% grid on
% xlabel('Mode 1')
% ylabel('Mode 2')
% zlabel('Mode 3')
