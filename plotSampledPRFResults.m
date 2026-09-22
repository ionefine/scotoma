% plotSampledPRFResults
% Load saved CSS and linear-shift results and recreate all figures.
% No models are refitted.

resultsDir = '.';
roiList = 1:3;
writeSummaryTables = false;
CSS = cell(1,max(roiList));
Shifts = cell(1,max(roiList));
CSSSummary = table();
LinearShiftSummary = table();
for roi = roiList
    cssFile = fullfile(resultsDir,sprintf('SampledPRFCSS_V%d.mat',roi));
    shiftFile = fullfile(resultsDir,sprintf('SampledLinearShifts_V%d.mat',roi));
    if ~isfile(cssFile), error('Missing file: %s',cssFile); end
    if ~isfile(shiftFile), error('Missing file: %s',shiftFile); end
    tmp = load(cssFile,'result');
    if ~isfield(tmp,'result') || ~isfield(tmp.result,'summaryTable')
        error('%s does not contain result.summaryTable.',cssFile);
    end
    CSS{roi} = tmp.result;
    if isfield(CSS{roi},'ROI') && isfinite(CSS{roi}.ROI) && CSS{roi}.ROI ~= roi
        error('%s contains ROI %g, not ROI %g.',cssFile,CSS{roi}.ROI,roi);
    end
    tmp = load(shiftFile,'result');
    if ~isfield(tmp,'result') || ~isfield(tmp.result,'summaryTable')
        error('%s does not contain result.summaryTable.',shiftFile);
    end
    Shifts{roi} = tmp.result;
    if isfield(Shifts{roi},'ROI') && isfinite(Shifts{roi}.ROI) && Shifts{roi}.ROI ~= roi
        error('%s contains ROI %g, not ROI %g.',shiftFile,Shifts{roi}.ROI,roi);
    end
    T = CSS{roi}.summaryTable;
    if ~ismember('ROI',T.Properties.VariableNames)
        T = addvars(T,repmat(roi,height(T),1),'Before',1,'NewVariableNames','ROI');
    end
    CSSSummary = [CSSSummary;T]; %#ok<AGROW>
    T = Shifts{roi}.summaryTable;
    if ismember('ROI',T.Properties.VariableNames)
        if any(T.ROI ~= roi)
            error('%s has an inconsistent ROI column.',shiftFile);
        end
    else
        T = addvars(T,repmat(roi,height(T),1),'Before',1,'NewVariableNames','ROI');
    end
    LinearShiftSummary = [LinearShiftSummary;T]; %#ok<AGROW>
    plotted = T.ecc >= 0 & T.ecc <= 5 & T.nSubjects > 0;
    if any(T.nNoFillFailed(plotted) > 0 | T.nWithKFailed(plotted) > 0)
        warning('plotSampledPRFResults:nonconvergence', ...
            'One or more plotted V%d bins contain non-converged subject fits.',roi);
    end
    if any(T.nFixedAtBoundary(plotted) > 0 | T.nNoFillAtBoundary(plotted) > 0 | ...
           T.nWithKAtBoundary(plotted) > 0)
        warning('plotSampledPRFResults:shiftAtBoundary', ...
            ['One or more plotted V%d shift bins have subject estimates pinned at a ' ...
             'parameter bound; Figures 3 and 4 understate their uncertainty.'],roi);
    end
    % The same checks for the CSS side, which Figures 1 and 2 draw. Without
    % these a bin where every subject''s k pinned at 0 or 1 plots silently.
    C = CSS{roi}.summaryTable;
    cssPlotted = C.ecc >= 0 & C.ecc <= 5 & C.nSubjects > 0;
    if any(C.nLinearKAtBoundary(cssPlotted) > 0 | C.nCSSKAtBoundary(cssPlotted) > 0)
        warning('plotSampledPRFResults:cssKAtBoundary', ...
            ['One or more plotted V%d CSS bins have subject k estimates pinned at a ' ...
             'k bound; Figure 1 understates their uncertainty.'],roi);
    end
end

if writeSummaryTables
    disp(CSSSummary)
    disp(LinearShiftSummary)
    writetable(CSSSummary,fullfile(resultsDir,'SampledPRFCSS_summary.csv'));
    writetable(LinearShiftSummary,fullfile(resultsDir,'SampledLinearShifts_summary.csv'));
end

% Figure 1: k estimated with independently fitted linear and CSS pRFs.
figure('Color','w','Name','Linear versus CSS k','Position',[80 80 360*numel(roiList) 360]);
tiledlayout(1,numel(roiList),'TileSpacing','compact','Padding','compact');
for roi = roiList
    T = CSS{roi}.summaryTable;
    nexttile; hold on
    q = T.ecc >= 0 & T.ecc <= 5;
    h1 = errorbar(T.ecc(q),T.kLinear(q),T.kLinear_se(q), ...
        'ko-','MarkerFaceColor','k','LineWidth',1.3);
    h2 = errorbar(T.ecc(q),T.kCSS(q),T.kCSS_se(q), ...
        'o-','Color',[0.84 0.37 0],'MarkerFaceColor',[0.84 0.37 0],'LineWidth',1.3);
    xline(2,'k--'); yline(0,'k:');
    xlabel('original pRF eccentricity (deg)'); ylabel('k');
    title(sprintf('V%d: linear versus CSS',roi)); xlim([0 5]); ylim([0 1]);
    grid on; set(gca,'TickDir','out','Box','off');
    if roi == roiList(1), legend([h1 h2],{'linear pRF','CSS pRF'},'Location','best','Box','off'); end
end

% Figure 5: original fixed-pRF k estimates in the three visual areas.
roiColors = [1 0 0;0 0.72 0;0 0 1];
figure('Color','w','Name','Fixed pRF k by visual area','Position',[100 100 760 500]);
hold on
h = gobjects(numel(roiList),1);
for j = 1:numel(roiList)
    roi = roiList(j);
    T = Shifts{roi}.summaryTable;
    q = T.ecc >= 0 & T.ecc <= 5 & isfinite(T.k_fixedPRF);
    h(j) = errorbar(T.ecc(q),T.k_fixedPRF(q),T.k_fixedPRF_se(q),'o-', ...
        'Color',roiColors(roi,:),'MarkerFaceColor',roiColors(roi,:), ...
        'MarkerSize',6,'LineWidth',1.3,'CapSize',5);
end
xline(2,'k--'); yline(0,'k:'); yline(1,'k:');
xlabel('pRF eccentricity (deg)'); ylabel('filling-in factor k');
xlim([0 5]);
roiLabels = arrayfun(@(r) sprintf('V%d',r),roiList,'UniformOutput',false);
legend(h,roiLabels,'Location','best','Box','off');
grid on; set(gca,'TickDir','out','Box','off');

% Figure 6: bottom-up and predictive components of the scotoma response.
% Each panel uses the same ROI color as the k graph. Bottom-up drive is the
% darker band; the predictive addition is the lighter band stacked above it.
figure('Color','w','Name','Response decomposition by visual area', ...
    'Position',[80 100 1120 380]);
tiledlayout(1,numel(roiList),'TileSpacing','compact','Padding','compact');
for j = 1:numel(roiList)
    roi = roiList(j);
    T = Shifts{roi}.summaryTable;
    q = T.ecc >= 0 & T.ecc <= 5 & isfinite(T.massIn_median) & ...
        isfinite(T.k_fixedPRF);
    ecc = T.ecc(q);
    m = T.massIn_median(q);
    k = T.k_fixedPRF(q);
    bottomUp = 1-m;
    predictive = k.*m;
    total = bottomUp+predictive;
    baseColor = roiColors(roi,:);
    darkColor = 0.55*baseColor;
    lightColor = 0.65+0.35*baseColor;
    ax = nexttile; hold(ax,'on');
    hBottom = addBand(ax,ecc,zeros(size(bottomUp)),bottomUp,darkColor,0.85);
    hPredictive = addBand(ax,ecc,bottomUp,total,lightColor,0.90);
    hTotal = plot(ax,ecc,total,'o-','Color',baseColor, ...
        'MarkerFaceColor',baseColor,'MarkerSize',5,'LineWidth',1.6);
    yline(ax,1,'k-','LineWidth',1);
    xline(ax,2,'k--','LineWidth',1);
    xlabel(ax,'pRF eccentricity (deg)');
    if j == 1, ylabel(ax,'response / full-field response'); end
    title(ax,sprintf('V%d',roi));
    xlim(ax,[0 5]);
    grid(ax,'on'); set(ax,'TickDir','out','Box','off','Layer','top');
    if j == 1
        legend(ax,[hBottom hPredictive hTotal], ...
            {'bottom-up','predictive','total scotoma response'}, ...
            'Location','best','Box','off');
    end
end

% Figure 2: fitted CSS exponent and effective-size change.
figure('Color','w','Name','CSS diagnostics','Position',[80 480 360*numel(roiList) 360]);
tiledlayout(1,numel(roiList),'TileSpacing','compact','Padding','compact');
for roi = roiList
    T = CSS{roi}.summaryTable;
    nexttile; hold on
    yyaxis left
    q = T.ecc >= 0 & T.ecc <= 5;
    errorbar(T.ecc(q),T.medianCssN(q),T.medianCssN_se(q),'o-', ...
        'Color',[0.49 0.18 0.56], ...
        'MarkerFaceColor',[0.49 0.18 0.56],'LineWidth',1.5);
    ylabel('subject-median CSS n'); ylim([0 1.35]);
    yyaxis right
    errorbar(T.ecc(q),T.medianEffectiveSizeRatio(q), ...
        T.medianEffectiveSizeRatio_se(q),'s-','Color',[0 0.45 0.70], ...
        'MarkerFaceColor',[0 0.45 0.70],'LineWidth',1.5);
    ylabel('effective CSS size / linear size');
    xline(2,'k--'); yline(1,'k:');
    xlabel('original pRF eccentricity (deg)'); title(sprintf('V%d: CSS diagnostics',roi));
    xlim([0 5]);
    grid on; set(gca,'TickDir','out','Box','off');
end

% Figure 3: linear-pRF shifts and size changes, with and without k.
colNo = [0.84 0.37 0]; colK = [0 0.45 0.70]; colNull = [0.55 0.55 0.55];
figure('Color','w','Name','Linear pRF shifts','Position',[80 80 360*numel(roiList) 700]);
tiledlayout(2,numel(roiList),'TileSpacing','compact','Padding','compact');
for roi = roiList
    T = Shifts{roi}.summaryTable;
    nexttile; hold on
    yline(0,'Color',[0.82 0.82 0.82]); xline(2,'k--');
    q = T.ecc >= 0 & T.ecc <= 5;
    h0 = errorbar(T.ecc(q),T.dtheta_noFill(q),T.dtheta_noFill_se(q), ...
        '-','Color',colNull,'LineWidth',0.8);
    h1 = errorbar(T.ecc(q),T.dr_noFill(q),T.dr_noFill_se(q), ...
        'o-','Color',colNo,'MarkerFaceColor',colNo,'LineWidth',1.4);
    h2 = errorbar(T.ecc(q),T.dr_withK(q),T.dr_withK_se(q), ...
        'o-','Color',colK,'MarkerFaceColor',colK,'LineWidth',1.4);
    title(sprintf('V%d radial shift',roi)); ylabel('\Delta r (deg; - toward scotoma)');
    xlabel('linear-pRF eccentricity (deg)'); grid on; set(gca,'TickDir','out','Box','off');
    xlim([0 5]);
    if roi == roiList(1)
        legend([h1 h2 h0],{'no filling','with k','tangential null'},'Location','best','Box','off');
    end
end
for roi = roiList
    T = Shifts{roi}.summaryTable;
    nexttile; hold on
    yline(0,'Color',[0.82 0.82 0.82]); xline(2,'k--');
    q = T.ecc >= 0 & T.ecc <= 5;
    errorbar(T.ecc(q),T.dsigma_noFill(q),T.dsigma_noFill_se(q), ...
        'o-','Color',colNo, ...
        'MarkerFaceColor',colNo,'LineWidth',1.4);
    errorbar(T.ecc(q),T.dsigma_withK(q),T.dsigma_withK_se(q), ...
        'o-','Color',colK, ...
        'MarkerFaceColor',colK,'LineWidth',1.4);
    title(sprintf('V%d size change',roi)); ylabel('\Delta\sigma (deg; + larger)');
    xlabel('linear-pRF eccentricity (deg)'); grid on; set(gca,'TickDir','out','Box','off');
    xlim([0 5]);
end

% Figure 4: how allowing shifts changes the linear-model estimate of k.
figure('Color','w','Name','Fixed versus joint k','Position',[80 80 360*numel(roiList) 360]);
tiledlayout(1,numel(roiList),'TileSpacing','compact','Padding','compact');
for roi = roiList
    T = Shifts{roi}.summaryTable;
    nexttile; hold on
    q = T.ecc >= 0 & T.ecc <= 5;
    h1 = errorbar(T.ecc(q),T.k_fixedPRF(q),T.k_fixedPRF_se(q), ...
        'ko-','MarkerFaceColor','k','LineWidth',1.3);
    h2 = errorbar(T.ecc(q),T.k_joint(q),T.k_joint_se(q), ...
        'o-','Color',colK,'MarkerFaceColor',colK,'LineWidth',1.3);
    xline(2,'k--'); yline(0,'k:');
    xlabel('linear-pRF eccentricity (deg)'); ylabel('k');
    title(sprintf('V%d: effect of allowing shifts',roi));
    xlim([0 5]);
    grid on; set(gca,'TickDir','out','Box','off');
    if roi == roiList(1), legend([h1 h2],{'fixed pRF','joint shift + k'},'Location','best','Box','off'); end
end

function h = addBand(ax,x,y1,y2,color,alpha)
valid = isfinite(x) & isfinite(y1) & isfinite(y2);
starts = find(valid & [true;~valid(1:end-1)]);
stops = find(valid & [~valid(2:end);true]);
h = patch(ax,nan,nan,color,'EdgeColor','none','FaceAlpha',alpha);
found = false;
for i = 1:numel(starts)
    idx = starts(i):stops(i);
    if numel(idx) < 2, continue, end
    p = patch(ax,[x(idx);flipud(x(idx))],[y1(idx);flipud(y2(idx))],color, ...
        'EdgeColor','none','FaceAlpha',alpha);
    if ~found, delete(h); h = p; found = true; end
end
end
