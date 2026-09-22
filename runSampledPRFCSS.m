% runSampledPRFCSS
% Use independent study-1 x/y locations. Fit linear sigma/beta and CSS
% sigma/n/beta from the three current full-field runs, then hold them fixed
% while k is estimated from the three scotoma runs.
%
% Retained original pRFs: requested ROI, centre eccentricity 0-5 deg,
% variance explained > 0.20, and sigma >= 0.10 deg. Both newly fitted pRF
% models must be usable for a voxel to enter the linear-versus-CSS k
% comparison. See fitSampledPRFCSS for the complete rules.
% The full-field fits and subsequent scotoma k fits use the same retained
% voxels, but the scotoma responses never influence sigma, n, or beta.

requiredFunctions = {'compileStimAndSubData','subData','loadScotomaStimuli', ...
                     'loadHRFParams','convHRF','convByHemisphere', ...
                     'hrf_twogamma','fitSampledPRFCSS'};
for i = 1:numel(requiredFunctions)
    if exist(requiredFunctions{i},'file') ~= 2
        error('runSampledPRFCSS:missingFunction', ...
              'Add the folder containing %s.m to the MATLAB path.',requiredFunctions{i});
    end
end

prfCodeDir = 'C:\Users\Ione Fine\Documents\code\scotoma\pRF-master';
if exist('fminsearchcon','file') ~= 2
    if ~isfolder(prfCodeDir)
        error('runSampledPRFCSS:missingPRFCode', ...
            'fminsearchcon is not on the path and pRF-master was not found at %s.',prfCodeDir);
    end
    optimizerFile = dir(fullfile(prfCodeDir,'**','fminsearchcon.m'));
    if isempty(optimizerFile)
        error('runSampledPRFCSS:missingOptimizer','Could not find fminsearchcon inside pRF-master.');
    end
    addpath(optimizerFile(1).folder);
end

fitOpts = struct();
fitOpts.eccEdges = 0:0.25:5;
fitOpts.nPerBin = Inf; % use every eligible voxel in each eccentricity bin
fitOpts.minVoxPerBin = 5;
fitOpts.rngSeed = 1;
fitOpts.sigmaBounds = [0.05 8];
fitOpts.nBounds = [0.05 1.35]; % keep n = 1 inside the CSS parameter space
fitOpts.cssStarts = [0.12 0.33 0.75 1.00 1.20];
fitOpts.minFullR2 = -Inf; % no second R2 cutoff beyond original vexpl > .20
fitOpts.kBounds = [0 1];
fitOpts.kGridStep = 0.05;
fitOpts.boundaryTol = 0.01;
fitOpts.useParallel = true;
fitOpts.verbose = true;

roiList = 1:3; % use 1 for a quick V1 test, then 1:3 for the full analysis
SampledPRFCSS = cell(max(roiList),1);
summaryBlocks = cell(max(roiList),1);
voxelBlocks = cell(max(roiList),1);
for roi = roiList
    compileOpts = struct();
    compileOpts.ROI = roi;
    compileOpts.radRange = [0 5];
    compileOpts.minvexpl = 0.2;
    compileOpts.minSigma = 0.1;
    compileOpts.subList = 1:10;
    compileOpts.edgeSigma = 0;
    [subjectData,stimData,x,y] = compileStimAndSubData(compileOpts);
    SampledPRFCSS{roi} = fitSampledPRFCSS(subjectData,stimData,x,y,roi,fitOpts);
    summaryBlocks{roi} = SampledPRFCSS{roi}.summaryTable;
    voxelBlocks{roi} = SampledPRFCSS{roi}.voxelTable;
    result = SampledPRFCSS{roi};
    save(sprintf('SampledPRFCSS_V%d.mat',roi),'result','fitOpts','compileOpts','-v7.3');
end
SampledPRFCSSSummary = vertcat(summaryBlocks{roiList});
hasVoxelTable = cellfun(@(T) istable(T) && width(T) > 0,voxelBlocks(roiList));
if any(hasVoxelTable)
    selectedBlocks = voxelBlocks(roiList);
    SampledPRFCSSVoxels = vertcat(selectedBlocks{hasVoxelTable});
else
    SampledPRFCSSVoxels = table();
end
save('SampledPRFCSS_All.mat','SampledPRFCSS','SampledPRFCSSSummary', ...
     'SampledPRFCSSVoxels','fitOpts','-v7.3');
disp(SampledPRFCSSSummary)
fprintf('\nCSS fits saved. Next run runSampledLinearShifts, then plotSampledPRFResults.\n');
