% VoxelDemingSlopeMapStuff
%
% Estimate Deming-grid k separately in each eccentricity bin and subject.
% Also save the response-space effects r = slope-(1-m) and rExact =
% slope-(the simulated k=0 slope). This script does not estimate the
% fixed-pRF regression k or k2; those are reported by prfShiftFigure.

compileRequest.ROI = 1;
compileRequest.radRange = [0,4];
compileRequest.minvexpl = 0.2;
compileRequest.minSigma = 0.1;
compileRequest.subList = 1:10;
compileRequest.edgeSigma = 0;

% Reuse compiled data only when it was made with this exact request.
needCompile = ~exist('subjectData','var') || ~exist('stimData','var') || ...
              ~exist('x','var') || ~exist('y','var') || ...
              ~exist('compiledOpts','var') || ~isequaln(compiledOpts,compileRequest) || ...
              ~iscell(stimData) || numel(stimData) ~= numel(subjectData);
if needCompile
    [subjectData,stimData,x,y] = compileStimAndSubData(compileRequest);
    compiledOpts = compileRequest;
end
compileOpts = compileRequest;

noiseSD = 0;
dR = 0.2;
radRangeList = 0:dR:4;
radX = radRangeList(1:end-1)+dR/2;
nRads = numel(radX);
nSub = numel(subjectData);
if ~all(cellfun(@(S) isfield(S,'subNum'),subjectData))
    error('VoxelDemingSlopeMapStuff:missingSubjectID','Each subject needs a subNum field.');
end
subjectIDs = cellfun(@(S) S.subNum,subjectData(:));
sigmaFac = 1;                 % robustness check for assumed pRF size
fillGrid = -0.1:0.05:1;

simOpts = struct('rngSeed',1,'verbose',false,'fillMeasure','k');
demingOpts = struct('lambda',1,'forceZeroIntercept',true, ...
                    'centerData',true,'combineRuns','concat');

bestFill = nan(nSub,nRads);
rMean = nan(nSub,nRads);
rExactMean = nan(nSub,nRads);
nVoxUsed = zeros(nSub,nRads);
atGridBoundary = false(nSub,nRads);

for s = 1:nSub
    figure(s); clf; hold on
    xline(2,'--k','LineWidth',1.5);
    xlabel('pRF eccentricity (deg)');
    ylabel('Voxel Deming slope');
    title(sprintf('Subject %d: measured and fitted slopes',s));
    grid on; set(gca,'YLim',[-0.25,1.25],'XLim',[0,4],'XTick',radRangeList)
end

for radNum = 1:nRads
    subsetOpts = compileOpts;
    subsetOpts.radRange = radRangeList(radNum:radNum+1);
    for s = 1:nSub
        fprintf('Bin %d/%d, subject %d/%d\n',radNum,nRads,s,nSub)
        subjectDataSub = subsetSubjectDataByVoxel(subjectData{s},stimData{s},subsetOpts);
        nVox = numel(subjectDataSub.sigma);
        if nVox < 3
            warning('VoxelDemingSlopeMapStuff:tooFewVoxels', ...
                    'Skipping bin %d, subject %d: only %d voxels.',radNum,s,nVox);
            continue
        end

        % Recalculate unit-area pRFs after applying the sigma robustness factor.
        G = zeros(size(subjectDataSub.Gprf));
        for v = 1:nVox
            pRF.center = subjectDataSub.prfXY(v,:);
            pRF.sig = subjectDataSub.sigma(v)*sigmaFac;
            pRF.ar = 1;
            G(:,v) = Gauss(pRF,x,y,1);
        end
        area = sum(G,1);
        if any(~isfinite(G(:))) || any(~isfinite(area) | area <= 0)
            error('VoxelDemingSlopeMapStuff:badGaussian', ...
                  'Invalid pRF mass in bin %d, subject %d.',radNum,s);
        end
        subjectDataSub.Gprf = G./area;

        fo = fitFillFracFromDemingMap(subjectDataSub,stimData{s},fillGrid, ...
                                      noiseSD,simOpts,demingOpts);
        bestFill(s,radNum) = fo.bestFillFrac;
        nVoxUsed(s,radNum) = nnz(fo.commonVoxels);
        atGridBoundary(s,radNum) = ismember(fo.bestIndex,[1,numel(fillGrid)]);
        q = isfinite(fo.r);
        if any(q), rMean(s,radNum) = mean(fo.r(q)); end
        q = isfinite(fo.rExact);
        if any(q), rExactMean(s,radNum) = mean(fo.rExact(q)); end

        % Diagnostic overlay using the same options and model as the grid fit.
        simData = simulateSubjectDataWithFilling({subjectDataSub},stimData(s), ...
                                                  noiseSD,fo.bestFillFrac,simOpts);
        measuredSlope = voxelDemingSlopeMap(subjectDataSub,demingOpts);
        fittedSlope = voxelDemingSlopeMap(simData{1},demingOpts);
        ecc = hypot(subjectDataSub.prfXY(:,1),subjectDataSub.prfXY(:,2));
        figure(s)
        plot(ecc,measuredSlope,'o','MarkerFaceColor','b','MarkerEdgeColor','none','MarkerSize',3);
        plot(ecc,fittedSlope,'o','MarkerFaceColor','r','MarkerEdgeColor','none','MarkerSize',3);
        drawnow
    end
end

meanFill = mean(bestFill,1,'omitnan');
semFill = std(bestFill,0,1,'omitnan')./sqrt(sum(isfinite(bestFill),1));
meanR = mean(rMean,1,'omitnan');
semR = std(rMean,0,1,'omitnan')./sqrt(sum(isfinite(rMean),1));
meanRExact = mean(rExactMean,1,'omitnan');
semRExact = std(rExactMean,0,1,'omitnan')./sqrt(sum(isfinite(rExactMean),1));
figure(nSub+1); clf; hold on
for radNum = 1:nRads
    plot(radX(radNum),bestFill(:,radNum),'k.');
end
errorbar(radX,meanFill,semFill,'k','LineStyle','none');
plot(radX,meanFill,'ko-','MarkerFaceColor','b');
xline(2,'--k','LineWidth',1.5);
xlabel('Eccentricity (deg)'); ylabel('Filling-in fraction (k)');
set(gca,'YLim',[-0.1,1.1],'XLim',[0,max(radX)+0.1]); grid on

fileName = sprintf('V%d_sigmaFac_%g.mat',compileOpts.ROI,sigmaFac);
save(fileName,'radX','radRangeList','bestFill','meanFill','semFill', ...
              'rMean','meanR','semR','rExactMean','meanRExact','semRExact', ...
              'nVoxUsed','atGridBoundary','fillGrid','sigmaFac', ...
              'subjectIDs','compileOpts','simOpts','demingOpts')
