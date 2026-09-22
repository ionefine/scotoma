% runInferCSSKCurves
%
% For each ROI, use the observed fixed-pRF eccentricity-by-k curve from
% runKResponseFigure as the target. Simulate CSS responses using a Gaussian
% narrowed by sqrt(n), then find the CSS k in each bin that reproduces the
% target when the responses are fitted with the original linear pRF.

clear CSSInference
kayN = [0.33 0.10 0.065];
effectiveSizeScale = [1 1 1]; % CSS effective size / fitted linear pRF size
kGrid = 0:0.005:1;
simOpts = struct('TR',1.2,'maxVoxPerBin',25,'minVox',5, ...
                 'rngSeed',1,'tolK',0.001,'verbose',true);
CSSInference = cell(3,1);
summaryBlocks = cell(3,1);

for roi = 1:3
    resultFile = sprintf('KResponseFigure_V%d.mat',roi);
    if ~isfile(resultFile)
        error('runInferCSSKCurves:missingTarget', ...
              'Run runKResponseFigure for V%d first; %s is missing.',roi,resultFile);
    end
    D = load(resultFile,'T','eccEdges');
    if ~isfield(D,'T') || ~istable(D.T) || ...
       ~all(ismember({'ecc','k'},D.T.Properties.VariableNames)) || ...
       ~isfield(D,'eccEdges')
        error('runInferCSSKCurves:badTarget','%s lacks T.ecc, T.k, or eccEdges.',resultFile);
    end
    targetK = double(D.T.k(:));
    eccEdges = double(D.eccEdges(:).');
    expectedEcc = ((eccEdges(1:end-1)+eccEdges(2:end))/2).';
    if numel(targetK) ~= numel(expectedEcc) || ...
       any(abs(double(D.T.ecc(:))-expectedEcc) > 1e-10)
        error('runInferCSSKCurves:badTarget', ...
              'The rows of T do not match the saved eccentricity bins in %s.',resultFile);
    end

    compileOpts.ROI = roi;
    compileOpts.radRange = [0 4];
    compileOpts.minvexpl = 0.2;
    compileOpts.minSigma = 0.1;
    compileOpts.subList = 1:10;
    compileOpts.edgeSigma = 0;
    [subjectData,stimData] = compileStimAndSubData(compileOpts);
    simOpts.effectiveSizeScale = effectiveSizeScale(roi);
    CSSInference{roi} = inferCSSKFromLinearCurve(subjectData,stimData, ...
                              kayN(roi),eccEdges,targetK,kGrid,simOpts);
    R = CSSInference{roi}.table;
    R.ROI = repmat(roi,height(R),1);
    R = movevars(R,'ROI','Before','ecc');
    summaryBlocks{roi} = R;
    disp(R)
end

CSSInferenceTable = vertcat(summaryBlocks{:});
figure('Color','w','Position',[100 100 1120 360]);
tiledlayout(1,3,'TileSpacing','compact');
for roi = 1:3
    R = CSSInference{roi};
    nexttile; hold on
    plot(R.ecc,R.targetLinearK,'ko-','MarkerFaceColor','k','LineWidth',1.5);
    plot(R.ecc,R.inferredCSSK,'o-','Color',[0.84 0.37 0], ...
         'MarkerFaceColor',[0.84 0.37 0],'LineWidth',1.5);
    plot(R.ecc,R.verifiedLinearK,'--','Color',[0 0.45 0.70],'LineWidth',1.2);
    xline(2,'k--'); yline(0,'k:');
    xlabel('pRF eccentricity (deg)'); ylabel('k');
    title(sprintf('V%d, n = %.3g, size = %.2g x linear', ...
                  roi,kayN(roi),effectiveSizeScale(roi)));
    ylim([0 1]); grid on; set(gca,'TickDir','out','Box','off');
    if roi == 1
        legend('observed linear k','inferred CSS k', ...
               'linear refit of simulation','Location','best','Box','off');
    end
end

save('InferredCSSKCurves.mat','CSSInference','CSSInferenceTable', ...
     'kayN','effectiveSizeScale','kGrid','simOpts')
