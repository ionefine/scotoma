% VoxelDemingSlopeMapStuff
%
% Generates 'k' statistics across eccentricity bins for each subject,
% computes meand and SD's and saves the results in matlab file
% 'V%X_sigmaFac_Y'  where X is the visual area ('V1, V2 or V3)' and Y is
% the sigma scale factor (1, 2 or 3).  Should be run separately for each
% combination of visual area and sigma scale factor.
%
% After this, run 'CompileKStats.m' to load in each of the 9 files and save
% results as a csv file in long format to be used for plotting in
% 'MakeScotomaFigure.m'  in the ../text/papers/scotoma/matlab' directory.

sigmaFac = 1;  % control to test robusteness of pRF size estimates


if ~exist('subjectData','var')
    compileOpts.ROI = 3;
    compileOpts.radRange = [0,4];
    compileOpts.minvexpl = 0.2;  %
    compileOpts.minSigma = 0.1;
    compileOpts.subList = 1:10;
    compileOpts.edgeSigma = 0;
    [subjectData,stimData,x,y] = compileStimAndSubData(compileOpts);
end
noiseSD = 0;


dR = .2;
radRangeList = [0:dR:2.4];

radX = radRangeList(1:end-1)+dR/2;
nRads = length(radRangeList)-1;
imSize = [108,108];


for i=1:10
    figure(i)
    clf
end

nSub = numel(subjectData);
bestFill          = nan(nSub, nRads);
demVal            = nan(nSub, nRads);
tmpSlope_ff_only  = nan(nSub, nRads);
fitOut            = cell(nSub, 1);

opts = struct();
opts.edgeNSigma         = 2;
opts.combineRuns        = 'concat';
opts.lambda             = 1;
opts.forceZeroIntercept = true;
opts.centerData         = true;

simOpts = struct();
simOpts.useVoxelGain = false;
simOpts.rngSeed      = 1;
simOpts.verbose      = false;

demingOpts = struct();
demingOpts.lambda             = 1;
demingOpts.forceZeroIntercept = true;
demingOpts.centerData         = true;
demingOpts.combineRuns        = 'concat';
demingOpts.verbose            = false;

fillGrid = fliplr([-.1:.05:1]);

for radNum = 1:nRads

    for s = 1:nSub
        fprintf('Subject %d of %d\n',s,nSub)
        compileOpts.radRange = [radRangeList(radNum),radRangeList(radNum+1)];
        subjectDataSub = subsetSubjectDataByVoxel(subjectData{s}, stimData, compileOpts);

        % scale pRF sigmas by sigmaFac (single application)
        subjectDataSub.sigma = subjectDataSub.sigma * sigmaFac;

        % recompute G with the scaled sigmas
        pRF = [];
        G = zeros(size(subjectDataSub.Gprf));
        for i = 1:length(subjectDataSub.sigma)
            pRF.center = [subjectDataSub.prfXY(i,1), subjectDataSub.prfXY(i,2)];
            pRF.sig    = subjectDataSub.sigma(i);
            pRF.ar     = 1;
            G(:,i)     = Gauss(pRF, x, y, 1);
        end

        % Normalize the pRFs to have equal area (not height).
        G = G ./ repmat(sum(G), size(G,1), 1);
        subjectDataSub.Gprf = G;

        fitOut{s} = fitFillFracFromResidualDeming(subjectDataSub, stimData, fillGrid, noiseSD, simOpts, demingOpts);

        % deming slope predicted by feedforward only
        simData_ff_only = simulateSubjectDataWithFilling({subjectDataSub}, stimData, 0, 0, simOpts);
        tmpSlope_ff_only(s, radNum) = mean(voxelDemingSlopeMap(simData_ff_only{1}, demingOpts), 'omitnan');

        bestFill(s, radNum) = fitOut{s}.bestFillFrac;
        demVal(s,   radNum) = mean(fitOut{s}.realSlope, 'omitnan');

        fillFrac = bestFill(s, radNum);
        simData     = simulateSubjectDataWithFilling({subjectDataSub}, stimData, noiseSD, fillFrac, simOpts);
        simBetaVoxel = voxelDemingSlopeMap(simData{1}, demingOpts);

        figure(s)
        hold on
        ecc = sqrt(subjectDataSub.prfXY(:,1).^2 + subjectDataSub.prfXY(:,2).^2);

        plot(ecc, fitOut{s}.realSlope, 'wo', 'MarkerFaceColor','b', 'MarkerSize',3);
        plot(ecc, simBetaVoxel,        'wo', 'MarkerFaceColor','r', 'MarkerSize',3, 'MarkerEdgeColor','none');

        xline(2, '--k', 'LineWidth', 1.5);
        xlabel('pRF eccentricity');
        ylabel('Voxel Deming slope');
        title('Voxel-wise Deming slopes vs eccentricity');
        grid on;
        set(gca,'YLim',[-.25,1.25])
        set(gca,'XTick',radRangeList);
        drawnow
    end

    % per-radNum aggregate plots
    meanFill = mean(bestFill, 'omitnan');
    semFill  = std(bestFill, 0, 'omitnan') ./ sqrt(sum(isfinite(bestFill)));
    meanDem  = mean(demVal,   'omitnan');
    semDem   = std(demVal,   0, 'omitnan') ./ sqrt(sum(isfinite(demVal)));

    figure(11)
    clf
    hold on
    for i = 1:nRads
        plot(radX(i), bestFill(:,i), 'k.')
    end
    errorbar(radX, meanFill, semFill, 'k', 'LineStyle','none')
    plot(radX, meanFill, 'ko-', 'MarkerFaceColor','b')
    set(gca,'YLim',[-.1,1.1])
    xlabel('Eccentricity (deg)')
    ylabel('Filling In Factor (k)')
    set(gca,'XLim',[0, max(radX)+.1]);

    figure(12)
    clf
    hold on
    for i = 1:nRads
        plot(radX(i), demVal(:,i), 'k.')
    end
    errorbar(radX, meanDem, semDem, 'k', 'LineStyle','none')
    plot(radX, meanDem, 'ko-', 'MarkerFaceColor','b')
    set(gca,'YLim',[-.1,1.1])
    xlabel('Eccentricity (deg)')
    ylabel('Deming slope')
    set(gca,'XLim',[0, max(radX)+.1]);
    drawnow
end

fileName = sprintf('V%d_sigmaFac_%d',compileOpts.ROI,sigmaFac);
save(fileName,"radX","bestFill","compileOpts","opts","simOpts","demingOpts")

fileName = sprintf('V%d_deming_%d',compileOpts.ROI,sigmaFac);
save(fileName, 'meanDem','radX', 'semDem', 'meanFill', 'semFill','tmpSlope_ff_only');



