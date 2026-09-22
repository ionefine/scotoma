% runKResponseFigure  Fast fixed-pRF k/k2 and response decomposition.
%
% Unlike runPRFShiftFigure, this does not fit pRF shifts or size changes.
% It is therefore the appropriate runner for the bottom-up/predictive plot.

compileOpts.ROI = 3;                    % 1=V1, 2=V2, 3=V3
compileOpts.radRange = [0 4];
compileOpts.minvexpl = 0.2;
compileOpts.minSigma = 0.1;
compileOpts.subList = 1:10;
compileOpts.edgeSigma = 0;
[subjectData,stimData] = compileStimAndSubData(compileOpts);

eccEdges = 1.25:0.25:3;
nBoot = 2000;
T = fitFixedKByEcc(subjectData,stimData,eccEdges,nBoot);
[responseFigure,responseComponents] = plotKResponseDecomposition(T);
disp(T)
disp(responseComponents)

fileName = sprintf('KResponseFigure_V%d.mat',compileOpts.ROI);
save(fileName,'T','responseComponents','compileOpts','eccEdges','nBoot')
