% runSampledLinearShifts
% Direct nonlinear fits of pRF shifts/size changes, with and without k.
% Run runSampledPRFCSS first so the saved files contain independent x/y and
% full-field-fitted linear sigma/beta values. No iterative displacement or
% first-order derivative approximation is used here.

requiredFunctions = {'compileStimAndSubData','subData','loadScotomaStimuli', ...
                     'loadHRFParams','convHRF','convByHemisphere', ...
                     'hrf_twogamma','fitSampledLinearShifts'};
for i = 1:numel(requiredFunctions)
    if exist(requiredFunctions{i},'file') ~= 2
        error('runSampledLinearShifts:missingFunction', ...
              'Add the folder containing %s.m to the MATLAB path.',requiredFunctions{i});
    end
end

prfCodeDir = 'C:\Users\Ione Fine\Documents\code\scotoma\pRF-master';
if exist('fminsearchcon','file') ~= 2
    optimizerFile = dir(fullfile(prfCodeDir,'**','fminsearchcon.m'));
    if isempty(optimizerFile)
        error('runSampledLinearShifts:missingOptimizer', ...
              'Could not find fminsearchcon inside %s.',prfCodeDir);
    end
    addpath(optimizerFile(1).folder);
end

shiftOpts = struct();
shiftOpts.eccEdges = 0:0.25:5;
shiftOpts.minVoxPerSubject = 5;
shiftOpts.shiftBounds = [-2 2;-2 2;-2 3];
shiftOpts.minSigma = 0.05;
shiftOpts.minShiftEcc = 1e-6; % radial direction is undefined only at exact zero
shiftOpts.kBounds = [0 1];
shiftOpts.boundaryTol = 0.01;
shiftOpts.useParallel = true;
shiftOpts.verbose = true;

roiList = 1:3; % use 1 for a quick V1 test
SampledLinearShifts = cell(max(roiList),1);
summaryBlocks = cell(max(roiList),1);
for roi = roiList
    sampleFile = sprintf('SampledPRFCSS_V%d.mat',roi);
    if ~isfile(sampleFile)
        error('runSampledLinearShifts:missingSampleFile', ...
              '%s was not found. Run runSampledPRFCSS first.',sampleFile);
    end
    saved = load(sampleFile,'result');
    if ~isfield(saved,'result') || ~isfield(saved.result,'voxelTable')
        error('runSampledLinearShifts:badSampleFile','%s lacks result.voxelTable.',sampleFile);
    end
    compileOpts = struct('ROI',roi,'radRange',[0 5],'minvexpl',0.2, ...
                         'minSigma',0.1,'subList',1:10,'edgeSigma',0);
    [subjectData,stimData,x,y] = compileStimAndSubData(compileOpts);
    result = fitSampledLinearShifts(subjectData,stimData,x,y, ...
                                    saved.result.voxelTable,shiftOpts);
    SampledLinearShifts{roi} = result;
    summaryBlocks{roi} = result.summaryTable;
    save(sprintf('SampledLinearShifts_V%d.mat',roi), ...
         'result','shiftOpts','compileOpts','-v7.3');
end

SampledLinearShiftsSummary = table();
for roi = roiList
    T = summaryBlocks{roi};
    T = addvars(T,repmat(roi,height(T),1),'Before',1,'NewVariableNames','ROI');
    SampledLinearShiftsSummary = [SampledLinearShiftsSummary;T]; %#ok<AGROW>
end
save('SampledLinearShifts_All.mat','SampledLinearShifts', ...
     'SampledLinearShiftsSummary','shiftOpts','-v7.3');
disp(SampledLinearShiftsSummary)
fprintf('\nShift fits saved. Run plotSampledPRFResults to create the figures.\n');
