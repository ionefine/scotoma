function T = prfShiftFigure(src, eccRange, dR, nBoot, iterOpts)
% prfShiftFigure  Apparent pRF shift and size change, with and without
% filling-in in the model.
%
% Two panels compare three models plus a tangential null:
%   no filling term   fixed-beta apparent pRF shift/expansion
%   with k            joint shift + missing-stimulus model
%   with k2           joint shift + full-field-addition model
%   tangential null   must sit at zero for a trustworthy bin
%
% A gap between the first two is the change in estimated pRF geometry when
% the specified filling term is added; its interpretation depends on the
% fixed-beta and additive-component assumptions in prfShiftFit.
%
% RANGE.  Default 1.25 to 3 deg. Inside about 1.25 deg the scotoma leaves
% almost no bottom-up drive, so there is no amplitude to anchor the
% linearisation and the tangential null departs from zero. Beyond about
% 3 deg too little pRF mass falls inside the scotoma for k to be identified.
% Both limits are visible in the figure: the null is plotted, and the
% intervals widen.
%
% ITERATION.  A single linearisation saturates dsigma near 0.45 deg. Pass
% src as a struct {subjectData, stimData, x, y} to re-linearise and refit
% (prfShiftIterate), which removes that ceiling. Pass a cell array of
% prfShiftFit outputs instead to use the single-pass estimate.
%
% INPUTS
%   src       either a cell array of prfShiftFit outputs (single pass), or
%             struct('subjectData',{{...}},'stimData',{{...}},'x',x,'y',y)
%             to iterate; the two cell arrays must have matching subjects
%   eccRange  [lo hi] degrees      (default [1.25 3])
%   dR        bin width, degrees   (default 0.25)
%   nBoot     bootstrap replicates (default 2000)
%   iterOpts  passed to prfShiftIterate; .maxIter = 0 gives the original
%             single-pass linearisation
%
% OUTPUT
%   T   table of plotted values plus fixed-pRF k and k2 (primary definitions)
%       and k_joint/k2_joint from the geometry models. The fixed-pRF k2 is
%       the quantity that reduces to g-1 outside the scotoma.

    if nargin < 2 || isempty(eccRange), eccRange = [1.25 3]; end
    if nargin < 3 || isempty(dR),       dR       = 0.25;     end
    if nargin < 4 || isempty(nBoot),    nBoot    = 2000;     end
    if nargin < 5 || isempty(iterOpts), iterOpts = struct();  end
    if numel(eccRange) ~= 2 || any(~isfinite(eccRange)) || eccRange(2) <= eccRange(1)
        error('prfShiftFigure:badRange','eccRange must be [low high] with high > low.');
    end
    if ~isscalar(dR) || ~isfinite(dR) || dR <= 0
        error('prfShiftFigure:badBinWidth','dR must be positive.');
    end
    if ~isscalar(nBoot) || nBoot < 0 || nBoot ~= round(nBoot)
        error('prfShiftFigure:badBootstrap','nBoot must be a nonnegative integer.');
    end
    nB = round(diff(eccRange)/dR);
    if abs(nB*dR-diff(eccRange)) > 1e-10*max(1,diff(eccRange))
        error('prfShiftFigure:unevenBins','dR must divide the eccentricity range exactly.');
    end

    edges = linspace(eccRange(1),eccRange(2),nB+1);
    ctr   = (edges(1:end-1) + dR/2).';     % column, to match the estimates

    % Either use supplied single-pass fits, or iterate the linearisation.
    % Each model is iterated separately because its pRF estimate depends on
    % which filling term is present.
    isRaw = isstruct(src) && all(isfield(src,{'subjectData','stimData','x','y'}));
    if isRaw
        optsK = iterOpts; optsK.fillMeasure = 'k';
        optsK2 = iterOpts; optsK2.fillMeasure = 'k2';
        [fitFixed,~,~] = prfShiftIterate(src.subjectData,src.stimData, ...
                                        src.x,src.y,edges,4,optsK);
        % All models start at the same zero-offset pRFs. Reuse that costly
        % first fit rather than recomputing it three more times.
        optsK.initialFits = fitFixed;
        optsK2.initialFits = fitFixed;
        [fitNo,  dNo ] = prfShiftIterate(src.subjectData, src.stimData, ...
                                         src.x, src.y, edges, 1:3, optsK);
        [fitK, dK] = prfShiftIterate(src.subjectData, src.stimData, ...
                                    src.x, src.y, edges, 1:4, optsK);
        [fitK2, dK2] = prfShiftIterate(src.subjectData, src.stimData, ...
                                      src.x, src.y, edges, 1:4, optsK2);
    else
        fitNo = src; if ~iscell(fitNo), fitNo = {fitNo}; end
        fitFixed = fitNo; fitK = fitNo; fitK2 = fitNo;
        dNo = zeros(nB,3); dK = zeros(nB,3); dK2 = zeros(nB,3);
    end

    z = nan(nB,1);
    [drNo,drNoLo,drNoHi, dsNo,dsNoLo,dsNoHi] = deal(z);
    [drK,drKLo,drKHi,dsK,dsKLo,dsKHi,kJoint,kJointLo,kJointHi] = deal(z);
    [drK2,drK2Lo,drK2Hi,dsK2,dsK2Lo,dsK2Hi,k2Joint,k2JointLo,k2JointHi] = deal(z);
    [kFixed,kFixedLo,kFixedHi,k2Fixed,k2FixedLo,k2FixedHi] = deal(z);
    [dt,dtLo,dtHi,nFixed,nNo,nK,nK2,massMedian,massRMS] = deal(z);

    for i = 1:nB
        bN = cellfun(@(F) F.ecc > edges(i) & F.ecc <= edges(i+1), fitNo, ...
                     'UniformOutput', false);
        bKbin = cellfun(@(F) F.ecc > edges(i) & F.ecc <= edges(i+1), fitK, ...
                     'UniformOutput', false);
        bK2bin = cellfun(@(F) F.ecc > edges(i) & F.ecc <= edges(i+1), fitK2, ...
                     'UniformOutput', false);
        bF = cellfun(@(F) F.ecc > edges(i) & F.ecc <= edges(i+1), fitFixed, ...
                     'UniformOutput', false);

        % Match the voxel population across the three geometry models.
        mShift = cell(size(fitNo)); mFixed = cell(size(fitFixed));
        for s = 1:numel(fitNo)
            mShift{s} = bN{s} & bKbin{s} & bK2bin{s} & ...
                        fitNo{s}.ok & fitK{s}.ok & fitK2{s}.ok2;
            mFixed{s} = bF{s} & fitFixed{s}.ok & fitFixed{s}.ok2;
        end

        mv = [];
        for s = 1:numel(fitNo)
            if isfield(fitNo{s},'massIn')
                qMass = mShift{s} & isfinite(fitNo{s}.massIn);
                mv = [mv; fitNo{s}.massIn(qMass)]; %#ok<AGROW>
            end
        end
        if ~isempty(mv)
            massMedian(i) = median(mv);
            massRMS(i) = sqrt(mean(mv.^2));
        end

        [bFixed,cFixed,nFixed(i)] = prfShiftPool(fitFixed,4,mFixed,nBoot,0.05,'k');
        [b2Fixed,c2Fixed] = prfShiftPool(fitFixed,4,mFixed,nBoot,0.05,'k2');
        [bNo,cNo_,nNo(i)] = prfShiftPool(fitNo,1:3,mShift,nBoot,0.05,'k');
        [bK,cK,nK(i)] = prfShiftPool(fitK,1:4,mShift,nBoot,0.05,'k');
        [bK2,cK2,nK2(i)] = prfShiftPool(fitK2,1:4,mShift,nBoot,0.05,'k2');

        if nFixed(i) >= 20 && isfinite(bFixed) && isfinite(b2Fixed)
            kFixed(i)=bFixed; kFixedLo(i)=cFixed(1); kFixedHi(i)=cFixed(2);
            k2Fixed(i)=b2Fixed; k2FixedLo(i)=c2Fixed(1); k2FixedHi(i)=c2Fixed(2);
        end

        if nNo(i) >= 20 && all(isfinite(bNo))
            bNo(1:3) = bNo(1:3) + dNo(i,:).';
            cNo_(1:3,:) = cNo_(1:3,:) + dNo(i,:).';
            drNo(i)=bNo(1); drNoLo(i)=cNo_(1,1); drNoHi(i)=cNo_(1,2);
            dt(i)=bNo(2); dtLo(i)=cNo_(2,1); dtHi(i)=cNo_(2,2);
            dsNo(i)=bNo(3); dsNoLo(i)=cNo_(3,1); dsNoHi(i)=cNo_(3,2);
        end
        if nK(i) >= 20 && all(isfinite(bK))
            bK(1:3) = bK(1:3) + dK(i,:).';
            cK(1:3,:) = cK(1:3,:) + dK(i,:).';
            drK(i)=bK(1); drKLo(i)=cK(1,1); drKHi(i)=cK(1,2);
            dsK(i)=bK(3); dsKLo(i)=cK(3,1); dsKHi(i)=cK(3,2);
            kJoint(i)=bK(4); kJointLo(i)=cK(4,1); kJointHi(i)=cK(4,2);
        end
        if nK2(i) >= 20 && all(isfinite(bK2))
            bK2(1:3) = bK2(1:3) + dK2(i,:).';
            cK2(1:3,:) = cK2(1:3,:) + dK2(i,:).';
            drK2(i)=bK2(1); drK2Lo(i)=cK2(1,1); drK2Hi(i)=cK2(1,2);
            dsK2(i)=bK2(3); dsK2Lo(i)=cK2(3,1); dsK2Hi(i)=cK2(3,2);
            k2Joint(i)=bK2(4); k2JointLo(i)=cK2(4,1); k2JointHi(i)=cK2(4,2);
        end
    end

    T = table(ctr,massMedian,massRMS,nFixed,nNo,nK,nK2,drNo,drNoLo,drNoHi,drK,drKLo,drKHi, ...
              drK2,drK2Lo,drK2Hi,dsNo,dsNoLo,dsNoHi,dsK,dsKLo,dsKHi, ...
              dsK2,dsK2Lo,dsK2Hi,kFixed,kFixedLo,kFixedHi,kJoint,kJointLo,kJointHi, ...
              k2Fixed,k2FixedLo,k2FixedHi,k2Joint,k2JointLo,k2JointHi, ...
              dt,dtLo,dtHi, ...
        'VariableNames', {'ecc','massIn_median','massIn_rms', ...
        'nVox_fixedPRF','nVox_noFill','nVox_k','nVox_k2', ...
        'dr_noFill','dr_noFill_lo','dr_noFill_hi', ...
        'dr_k','dr_k_lo','dr_k_hi','dr_k2','dr_k2_lo','dr_k2_hi', ...
        'dsig_noFill','dsig_noFill_lo','dsig_noFill_hi', ...
        'dsig_k','dsig_k_lo','dsig_k_hi','dsig_k2','dsig_k2_lo','dsig_k2_hi', ...
        'k','k_lo','k_hi','k_joint','k_joint_lo','k_joint_hi', ...
        'k2','k2_lo','k2_hi','k2_joint','k2_joint_lo','k2_joint_hi', ...
        'dtheta_null','dtheta_null_lo','dtheta_null_hi'});

    % --- draw -------------------------------------------------------------
    colNo   = [0.84 0.37 0.00];    % vermillion, Okabe-Ito
    colK    = [0.00 0.45 0.70];    % blue
    colK2   = [0.00 0.62 0.45];    % bluish green
    colNull = [0.60 0.60 0.60];

    figure('Color','w','Position',[100 100 880 380]); clf

    panel = {{drNo,drNoLo,drNoHi,drK,drKLo,drKHi,drK2,drK2Lo,drK2Hi, ...
              'Apparent pRF shift','\delta r  (deg, + outward)',true}, ...
             {dsNo,dsNoLo,dsNoHi,dsK,dsKLo,dsKHi,dsK2,dsK2Lo,dsK2Hi, ...
              'Apparent pRF size change','\delta\sigma  (deg)',false}};

    for k = 1:2
        p = panel{k};
        ax = subplot(1,2,k); hold(ax,'on');
        yline(0, '-', 'Color', [0.85 0.85 0.85]);
        xline(2, '--', 'Color', [0.45 0.45 0.45], 'LineWidth', 1);

        h = gobjects(0); lbl = {};
        if p{12}     % the null belongs in one panel only
            hN = errorbar(ctr, dt, dt-dtLo, dtHi-dt, '-', 'Color', colNull, ...
                          'LineWidth', 0.75, 'CapSize', 0);
            h(end+1) = hN; lbl{end+1} = 'tangential (null)'; %#ok<AGROW>
        end
        h1 = errorbar(ctr, p{1}, p{1}-p{2}, p{3}-p{1}, 'o-', 'Color', colNo, ...
                      'MarkerFaceColor', colNo, 'MarkerSize', 5, ...
                      'LineWidth', 1.6, 'CapSize', 0);
        h2 = errorbar(ctr, p{4}, p{4}-p{5}, p{6}-p{4}, 'o-', 'Color', colK, ...
                      'MarkerFaceColor', colK, 'MarkerSize', 5, ...
                      'LineWidth', 1.6, 'CapSize', 0);
        h3 = errorbar(ctr, p{7}, p{7}-p{8}, p{9}-p{7}, 's-', 'Color', colK2, ...
                      'MarkerFaceColor', colK2, 'MarkerSize', 5, ...
                      'LineWidth', 1.4, 'CapSize', 0);
        h = [h1, h2, h3, h];
        lbl = [{'fixed \beta, no filling','joint model with k','joint model with k2'}, lbl];

        xlabel('pRF eccentricity (deg)');
        ylabel(p{11});
        title(p{10}, 'FontWeight', 'normal');
        xlim(eccRange);
        set(ax, 'TickDir','out', 'Box','off', 'Layer','top');

        v = [p{2}; p{3}; p{5}; p{6}; p{8}; p{9}];
        if p{12}, v = [v; dtLo; dtHi]; end
        v = v(isfinite(v));
        if ~isempty(v)
            lim = max(0.1, 1.1*max(abs(v)));
            ylim([-lim lim]);
        end

        yl = ylim;
        text(2, yl(1) + 0.04*diff(yl), ' scotoma edge', ...
             'FontSize', 8, 'Color', [0.45 0.45 0.45]);

        if k == 1
            legend(h, lbl, 'Location','southwest', 'Box','off', 'FontSize', 8);
        end
    end
end
