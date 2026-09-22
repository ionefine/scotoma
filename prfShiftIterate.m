function [fits, delta, hist] = prfShiftIterate(subjectData, stimData, x, y, ...
                                               edges, cols, opts)
% prfShiftIterate  Gauss-Newton refinement of the pRF shift estimate.
%
% WHY.  prfShiftFit linearises the pRF about its full-field parameters, so a
% single pass estimates only a first-order displacement. Measured against
% ground truth on V1, 1.25-3 deg, noiseless, the single pass gives:
%
%      true dsigma   0.05  0.10  0.20  0.30  0.50  0.75  1.00
%      recovered     0.05  0.10  0.18  0.25  0.36  0.42  0.45
%      true dr       0.05  0.10  0.20  0.30  0.50  0.75  1.00
%      recovered     0.05  0.11  0.22  0.34  0.59  0.90  1.18
%
% dr is well behaved -- within 10% below 0.2 deg, inflating to about +20%
% beyond, never saturating. dsigma SATURATES near 0.45 deg: above a true
% 0.5 deg the single pass can say an expansion exists but cannot size it.
%
% Re-linearising about the current estimate and refitting removes that
% compression, because each pass only ever has to estimate a small
% remaining increment. This is Gauss-Newton; at convergence it is the full
% nonlinear fit, at a fraction of the cost of a search.
%
% HOW.  The shift parameters are shared within an eccentricity bin, so each
% bin carries its own cumulative offset and every voxel is re-linearised
% about its own bin's current estimate. Each committed update costs one
% prfShiftFit per subject, followed by one final fit at the updated point.
%
% INPUTS
%   subjectData  cell array, one entry per subject
%   stimData     cell array matching subjectData; each entry contains
%                .SfullRaw{run}, .SscotRaw{run}, .hrfParams, .tStim, .TR
%   x, y         visual field grids
%   edges        eccentricity bin edges, e.g. 1.25:0.25:3
%   cols         model, as in prfShiftPool: 1:3 (no filling-in) or 1:4
%                (both). Passing 4 alone returns immediately -- the filling
%                term enters exactly, so it needs no iteration.
%   opts         .maxIter   committed Gauss-Newton updates; 0 gives the
%                           original single-pass estimate   (3)
%                .tol       stop when every bin's step is below this, deg
%                                                            (0.01)
%                .damp      step multiplier, < 1 to slow it  (1.0)
%                .maxStep   cap on one pass's step, deg      (0.5)
%                .fitOpts   passed through to prfShiftFit
%                .fillMeasure 'k' (default) or 'k2'; passed to prfShiftPool
%                .initialFits optional zero-offset fits to reuse on the
%                           first pass; avoids repeating an identical fit
%                .useParallel fit subjects with parfor       (false)
%
% OUTPUTS
%   fits   cell array of prfShiftFit outputs at the CONVERGED linearisation
%          point. Pool these to get the remaining increment and add that
%          increment to delta for the total.
%   delta  [nBin x 3] cumulative [dr dtheta dsigma] per bin, degrees
%   hist   [nBin x 3 x maxIter] the cumulative estimate after each update,
%          so convergence (or oscillation) is visible

    if nargin < 7 || isempty(opts), opts = struct(); end
    if ~isfield(opts,'maxIter'), opts.maxIter = 3;    end
    if ~isfield(opts,'tol'),     opts.tol     = 0.01; end
    if ~isfield(opts,'damp'),    opts.damp    = 1.0;  end
    if ~isfield(opts,'maxStep'), opts.maxStep = 0.5;  end
    if ~isfield(opts,'fitOpts'), opts.fitOpts = struct('minEcc',1e-6,'verbose',false); end
    if ~isfield(opts,'fillMeasure'), opts.fillMeasure = 'k'; end
    if ~isfield(opts,'initialFits'), opts.initialFits = []; end
    if ~isfield(opts,'useParallel'), opts.useParallel = false; end
    if numel(edges) < 2 || any(~isfinite(edges)) || any(diff(edges) <= 0)
        error('prfShiftIterate:badEdges','edges must be finite and strictly increasing.');
    end
    if isempty(cols) || any(~ismember(cols,1:4)) || numel(unique(cols)) ~= numel(cols)
        error('prfShiftIterate:badColumns','cols must contain unique values from 1:4.');
    end
    if ~any(strcmpi(opts.fillMeasure,{'k','k2'}))
        error('prfShiftIterate:badFillMeasure','opts.fillMeasure must be ''k'' or ''k2''.');
    end
    if ~isscalar(opts.maxIter) || opts.maxIter < 0 || opts.maxIter ~= round(opts.maxIter)
        error('prfShiftIterate:badIterations','opts.maxIter must be a nonnegative integer.');
    end
    if ~isscalar(opts.tol) || ~isfinite(opts.tol) || opts.tol < 0 || ...
       ~isscalar(opts.damp) || ~isfinite(opts.damp) || opts.damp <= 0 || ...
       ~isscalar(opts.maxStep) || ~isfinite(opts.maxStep) || opts.maxStep <= 0
        error('prfShiftIterate:badOptions','tol, damp, and maxStep are invalid.');
    end
    if ~(islogical(opts.useParallel) || isnumeric(opts.useParallel)) || ...
       ~isscalar(opts.useParallel) || ~ismember(opts.useParallel,[0 1])
        error('prfShiftIterate:badParallel','useParallel must be true or false.');
    end

    if ~iscell(subjectData), subjectData = {subjectData}; end
    nSub = numel(subjectData);
    if ~iscell(stimData) || numel(stimData) ~= nSub
        error('prfShiftIterate:stimulusCount', ...
              'stimData must contain one entry per subject.');
    end
    if ~isempty(opts.initialFits) && ...
       (~iscell(opts.initialFits) || numel(opts.initialFits) ~= nSub || ...
        ~all(cellfun(@(F) isstruct(F) && isscalar(F),opts.initialFits)))
        error('prfShiftIterate:badInitialFits', ...
              'initialFits must contain one fit struct per subject.');
    end
    nBin = numel(edges) - 1;

    % which entries of cols are the three pRF parameters
    [~, pIdx] = ismember(1:3, cols);
    if all(pIdx == 0)
        if isempty(opts.initialFits)
            fits = fitAll(subjectData,stimData,x,y,opts.fitOpts,[],opts.useParallel);
        else
            fits = opts.initialFits;
        end
        delta = zeros(nBin,3);
        hist  = zeros(nBin,3,1);
        return
    end

    delta  = zeros(nBin,3);
    hist   = nan(nBin,3,opts.maxIter);
    offset = cell(nSub,1);
    for s = 1:nSub, offset{s} = zeros(numel(subjectData{s}.sigma), 3); end

    if opts.maxIter == 0 && ~isempty(opts.initialFits)
        fits = opts.initialFits;
        return
    end

    for it = 1:opts.maxIter
        if it == 1 && ~isempty(opts.initialFits)
            fits = opts.initialFits;
        else
            fits = fitAll(subjectData,stimData,x,y,opts.fitOpts,offset,opts.useParallel);
        end

        step = zeros(nBin,3);
        for i = 1:nBin
            m = cellfun(@(F) F.ecc > edges(i) & F.ecc <= edges(i+1), fits, ...
                        'UniformOutput', false);
            [b,n] = prfShiftPool(fits,cols,m,opts.fillMeasure);
            if n < 20 || any(~isfinite(b)), continue, end
            for j = 1:3
                if pIdx(j) > 0, step(i,j) = b(pIdx(j)); end
            end
        end

        step = opts.damp * step;
        step = max(min(step, opts.maxStep), -opts.maxStep);
        step(~isfinite(step)) = 0;

        delta = delta + step;
        hist(:,:,it) = delta;

        for s = 1:nSub
            bin = discretize(sqrt(subjectData{s}.prfXY(:,1).^2 + ...
                                  subjectData{s}.prfXY(:,2).^2), edges);
            ok  = ~isnan(bin);
            offset{s}(:) = 0;
            offset{s}(ok,:) = delta(bin(ok),:);
        end

        if max(abs(step(:))) < opts.tol, break, end
    end

    % final pass at the converged point, so the returned cross-products
    % match the reported delta
    fits = fitAll(subjectData,stimData,x,y,opts.fitOpts,offset,opts.useParallel);
end

% =====================================================================
function fits = fitAll(subjectData,stimData,x,y,fitOpts,offset,useParallel)
    nSub = numel(subjectData);
    fits = cell(nSub,1);
    canParallel = useParallel && ~isempty(ver('parallel')) && ...
                  license('test','Distrib_Computing_Toolbox');
    if canParallel
        parfor s = 1:nSub
            o = fitOpts;
            if ~isempty(offset), o.offset = offset{s}; end
            fits{s} = prfShiftFit(subjectData{s},stimData{s},x,y,o);
        end
    else
        persistent warnedNoParallel
        if useParallel && isempty(warnedNoParallel)
            warning('prfShiftIterate:noParallelToolbox', ...
                    'Parallel toolbox unavailable; using a serial subject loop.');
            warnedNoParallel = true;
        end
        for s = 1:nSub
            o = fitOpts;
            if ~isempty(offset), o.offset = offset{s}; end
            fits{s} = prfShiftFit(subjectData{s},stimData{s},x,y,o);
        end
    end
end
