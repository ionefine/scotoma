% CompileKStats
%
% Compile the nine VoxelDemingSlopeMapStuff result files into allk.csv.
% k is the Deming-grid estimate; r and rExact are response-space effects.
% kResolved is NaN when the optimum lies on a search-grid boundary.

ROIList = 1:3;
sigmaList = 1:3;
blocks = cell(numel(ROIList)*numel(sigmaList),1);
b = 0;

for roi = ROIList
    for sigmaFactor = sigmaList
        fileName = sprintf('V%d_sigmaFac_%g.mat',roi,sigmaFactor);
        if ~isfile(fileName)
            error('CompileKStats:missingFile','Cannot find %s.',fileName);
        end
        D = load(fileName);
        required = {'radX','bestFill','rMean','rExactMean','nVoxUsed', ...
                    'atGridBoundary','subjectIDs'};
        missing = required(~isfield(D,required));
        if ~isempty(missing)
            error('CompileKStats:missingFields','%s lacks: %s.', ...
                  fileName,strjoin(missing,', '));
        end

        [nSub,nRad] = size(D.bestFill);
        expectedSize = [nSub,nRad];
        if nRad ~= numel(D.radX) || numel(D.subjectIDs) ~= nSub || ...
           ~isequal(size(D.rMean),expectedSize) || ...
           ~isequal(size(D.rExactMean),expectedSize) || ...
           ~isequal(size(D.nVoxUsed),expectedSize) || ...
           ~isequal(size(D.atGridBoundary),expectedSize)
            error('CompileKStats:badDimensions','Arrays in %s have inconsistent dimensions.',fileName);
        end

        subjectIDs = D.subjectIDs(:);
        n = nSub*nRad;
        subNum = repmat(subjectIDs,nRad,1);
        sub = arrayfun(@(v) sprintf('sub%02d',v),subNum,'UniformOutput',false);
        rad = repelem(D.radX(:),nSub);
        ROI = repmat({sprintf('V%d',roi)},n,1);
        sigmaFactorColumn = repmat(sigmaFactor,n,1);
        k = D.bestFill(:);
        atGridBoundary = logical(D.atGridBoundary(:));
        kResolved = k;
        kResolved(atGridBoundary) = NaN;
        r = D.rMean(:);
        rExact = D.rExactMean(:);
        nVox = D.nVoxUsed(:);

        b = b+1;
        blocks{b} = table(sub,subNum,rad,ROI,sigmaFactorColumn,k,kResolved, ...
                          r,rExact,nVox,atGridBoundary, ...
            'VariableNames',{'sub','subNum','rad','ROI','sigmaFactor','k', ...
                             'kResolved','r','rExact','nVox','atGridBoundary'});
    end
end
T = vertcat(blocks{:});
writetable(T,'allk.csv')

% Sanity-check plot. Boundary-limited estimates are shown in red but are
% excluded from the mean and SEM.
plotROI = 'V3';
plotSigma = 3;
plotMaxEcc = 2.4;
radList = unique(T.rad(strcmp(T.ROI,plotROI) & T.sigmaFactor == plotSigma));
radList = radList(radList <= plotMaxEcc);
m = nan(size(radList));
sem = nan(size(radList));
figure(1); clf; hold on
yl = [-0.1,1.1];
patch([2,plotMaxEcc,plotMaxEcc,2],[yl(1),yl(1),yl(2),yl(2)], ...
      [0.8,0.8,0.8],'LineStyle','none');
for i = 1:numel(radList)
    inBin = strcmp(T.ROI,plotROI) & T.sigmaFactor == plotSigma & T.rad == radList(i);
    y = T.kResolved(inBin);
    y = y(isfinite(y));
    if ~isempty(y)
        m(i) = mean(y);
        sem(i) = std(y)/sqrt(numel(y));
        plot(radList(i)*ones(size(y)),y,'k.');
    end
    onBoundary = inBin & T.atGridBoundary;
    plot(T.rad(onBoundary),T.k(onBoundary),'ro','MarkerSize',4);
end
errorbar(radList,m,sem,'k','LineStyle','none');
plot(radList,m,'ko-','MarkerFaceColor','k','MarkerSize',5);
yline(0,'k:'); xline(2,'k--');
xlabel('Eccentricity (deg)'); ylabel('Filling-in fraction (k)');
title(sprintf('%s, sigma factor %g',plotROI,plotSigma));
set(gca,'YLim',yl,'XLim',[0,plotMaxEcc],'XTick',0:0.2:plotMaxEcc)
