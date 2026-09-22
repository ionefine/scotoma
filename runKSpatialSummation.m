% runKSpatialSummation  Test whether fixed-pRF k changes under CSS.
%
% Change ROI and rerun for V1, V2, and V3. The Kay et al. (2013)
% exponents are approximate medians read from their Figure 7B.

clear subjectData stimData
compileOpts.ROI = 3;                    % 1=V1, 2=V2, 3=V3
compileOpts.radRange = [0 4];
compileOpts.minvexpl = 0.2;
compileOpts.minSigma = 0.1;
compileOpts.subList = 1:10;
compileOpts.edgeSigma = 0;
[subjectData,stimData] = compileStimAndSubData(compileOpts);

kayN = [0.33 0.10 0.065];
nGrid = [1 kayN(compileOpts.ROI)];      % use [1 .75 .5 .33 .2 .1 .065] for a curve
eccEdges = 1.25:0.25:3;
cssOpts = struct('TR',1.2,'kBounds',[0 1],'tolK',0.001, ...
                 'nCoarse',21,'minVox',20,'verbose',true);
CSSK = fitKSpatialSummation(subjectData,stimData,nGrid,eccEdges,cssOpts);
disp(CSSK.summary)
% 
% figure; clf
% C = lines(numel(nGrid));
% tiledlayout(2,1,'TileSpacing','compact');
% nexttile; hold on
% for j = 1:numel(nGrid)
%     errorbar(CSSK.ecc,CSSK.meanK(:,j),CSSK.semK(:,j),'o-', ...
%              'Color',C(j,:),'MarkerFaceColor',C(j,:));
% end
% yline(0,':k'); xlabel('pRF eccentricity (deg)'); ylabel('k');
% legend(compose('n = %.3g',nGrid),'Location','best'); grid on
% title(sprintf('V%d: fixed-pRF k with compressive spatial summation',compileOpts.ROI));
% nexttile; hold on
% for j = 1:numel(nGrid)
%     if nGrid(j) == 1, continue, end
%     errorbar(CSSK.ecc,CSSK.meanDeltaK(:,j),CSSK.semDeltaK(:,j),'o-', ...
%              'Color',C(j,:),'MarkerFaceColor',C(j,:));
% end
% yline(0,':k'); xlabel('pRF eccentricity (deg)'); ylabel('\Delta k from n=1');
% grid on

fileName = sprintf('KSpatialSummation_V%d.mat',compileOpts.ROI);
save(fileName,'CSSK','compileOpts','cssOpts','nGrid','eccEdges')
