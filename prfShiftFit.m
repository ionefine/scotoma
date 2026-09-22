function out = prfShiftFit(subjectData, stimData, x, y, opts)
% prfShiftFit  Per-voxel cross-products for the pRF shift / filling-in test.
%
% MODELS. The scotoma-condition timecourse is represented in two ways.
% The headline missing-stimulus coefficient k uses beta estimated from the
% full-field condition and holds it fixed in the scotoma condition:
%
%   Y ~ beta*[P + D*delta + k*K]
%
% The full-field-addition control k2 fixes beta at its value from the
% full-field condition:
%
%   Y ~ beta*[P + D*delta + k2*F]
%
%   P  = Sscot*G(p+offset)             surviving bottom-up drive
%   D  = Sscot*dG/dp                   local pRF derivatives
%   K  = (Sfull-Sscot)*G(p)            removed-stimulus drive
%   F  = Sfull*G(p)                     full-field drive
%
% G(p) in K and F is always the ORIGINAL full-field pRF. It must not move
% with the Gauss-Newton offset. k is identified only where K is appreciable. k2
% is identified everywhere F is nonzero, but outside the scotoma it measures
% a gain difference rather than filling-in.
%
% WHY CROSS-PRODUCTS.  P and K are strongly correlated (~0.95 near the
% scotoma border), so within one voxel's fit their coefficients trade off
% and their errors are anticorrelated. Estimating a ratio per voxel and
% averaging therefore measures that covariance and cancels to zero. The
% parameters must be estimated jointly across voxels. Both models subtract
% beta*P from the data. Projecting P out would instead make its amplitude
% free and would no longer estimate the k defined above; for k2 it would
% also destroy the outside-scotoma gain control.
%
% Selecting column subsets of the same cross-products gives the models
% without refitting anything:
%   columns 1:3  pRF free to shift and resize, no filling-in
%   column  4    pRF fixed, filling-in only
%   columns 1:4  both
%
% INPUTS
%   subjectData  .Yfull{run}, .Yscot{run}  [nT x nVox]
%                .hemIdx [nVox x 1]  1 = left hemisphere, 2 = right
%                .prfXY [nVox x 2], .sigma [nVox x 1]
%   stimData     .SfullRaw{run}, .SscotRaw{run}  [nT x nPix], UNCONVOLVED
%                .hrfParams (1x2, {L,R}), .tStim{run}, .TR
%   x, y         [ny x nx] visual field coordinate grids (funcOf.x/.y)
%   opts         .offset   [nVox x 3] of [dr dtheta dsigma] to linearise
%                          AROUND, instead of the original pRF. Zeros (the
%                          default) is the single-pass estimate. Supplying
%                          the previous pass's estimate and refitting is a
%                          Gauss-Newton step, which removes the compression
%                          the first-order approximation shows for large
%                          changes -- see prfShiftIterate.
%                .minEcc   skip only effectively zero eccentricity; polar
%                          angle and the radial/tangential split are
%                          undefined at the exact foveal centre         (1e-6)
%                .hPos     position step, fraction of sigma          (0.10)
%                .hSig     sigma step, fraction of sigma             (0.10)
%                .verbose                                            (true)
%
% OUTPUT  out, per voxel. Missing-stimulus model (k), beta fixed:
%   WtW, WtZ, ZtZ, dof, ok     columns [dr dtheta dsigma k]
% Full-field-addition model (k2), with beta fixed:
%   WtW2, WtZ2, ZtZ2, dof2, ok2  columns [dr dtheta dsigma k2]
%   beta [nVox x 1]       full-field amplitude
%   massIn [nVox x 1]     fraction of original pRF inside the scotoma
%   ecc, sigma, ok

    if nargin < 5 || isempty(opts), opts = struct(); end
    if ~isfield(opts,'minEcc'),  opts.minEcc  = 1e-6; end
    if ~isfield(opts,'hPos'),    opts.hPos    = 0.10; end
    if ~isfield(opts,'hSig'),    opts.hSig    = 0.10; end
    if ~isfield(opts,'verbose'), opts.verbose = true; end
    if ~isscalar(opts.minEcc) || ~isfinite(opts.minEcc) || opts.minEcc < 0 || ...
       ~isscalar(opts.hPos) || ~isfinite(opts.hPos) || opts.hPos <= 0 || ...
       ~isscalar(opts.hSig) || ~isfinite(opts.hSig) || opts.hSig <= 0
        error('prfShiftFit:badOptions','minEcc, hPos, and hSig are invalid.');
    end

    if ~all(isfield(stimData,{'SfullRaw','SscotRaw','hrfParams','tStim','TR'}))
        error('prfShiftFit:oldStimData', ...
              ['stimData must carry SfullRaw/SscotRaw, hrfParams, tStim and TR. ', ...
               'Pre-convolved Sfull/Sscot are no longer used: the HRF differs ', ...
               'between hemispheres, so the design is convolved per voxel here.']);
    end
    if ~isfield(subjectData,'hemIdx')
        error('prfShiftFit:missingHemisphere', ...
              'subjectData.hemIdx is required to pick each voxel''s HRF.');
    end
    nRun = numel(stimData.SfullRaw);
    if nRun == 0 || numel(stimData.SscotRaw) ~= nRun || ...
       numel(subjectData.Yfull) ~= nRun || numel(subjectData.Yscot) ~= nRun
        error('prfShiftFit:runMismatch','Stimulus and data run counts must agree.');
    end
    nVox = size(subjectData.Yscot{1}, 2);
    if ~isequal(size(x),size(y)) || numel(x) ~= size(stimData.SfullRaw{1},2)
        error('prfShiftFit:gridMismatch','x/y and stimulus pixel dimensions disagree.');
    end
    if size(subjectData.prfXY,1) ~= nVox || numel(subjectData.sigma) ~= nVox || ...
       numel(subjectData.hemIdx) ~= nVox
        error('prfShiftFit:voxelMismatch','pRF parameters and BOLD data have different voxel counts.');
    end
    hemIdx = double(subjectData.hemIdx(:));
    nHem = numel(stimData.hrfParams);
    if any(hemIdx < 1 | hemIdx > nHem)
        error('prfShiftFit:badHemIdx','hemIdx indexes outside stimData.hrfParams.');
    end
    TR = stimData.TR;

    if ~isfield(opts,'offset') || isempty(opts.offset)
        opts.offset = zeros(nVox, 3);
    end
    if ~isequal(size(opts.offset),[nVox,3]) || any(~isfinite(opts.offset(:)))
        error('prfShiftFit:badOffset','opts.offset must be a finite nVox-by-3 matrix.');
    end

    for r = 1:nRun
        if ~isequal(size(stimData.SfullRaw{r}),size(stimData.SscotRaw{r})) || ...
           size(stimData.SfullRaw{r},2) ~= numel(x) || ...
           size(stimData.SfullRaw{r},1) ~= size(subjectData.Yfull{r},1) || ...
           size(subjectData.Yfull{r},2) ~= nVox || ...
           ~isequal(size(subjectData.Yfull{r}),size(subjectData.Yscot{r}))
            error('prfShiftFit:dimensionMismatch','Dimensions disagree in run %d.',r);
        end
    end

    % --- concatenate runs, centring each run separately -------------------
    Yscot = catCentred(subjectData.Yscot);
    Yfull = catCentred(subjectData.Yfull);

    % Raw, unconvolved designs. Each voxel's regressors are convolved below
    % with that voxel's own hemisphere HRF, after projection onto its pRF.
    Ss = cell(nRun,1); Sd = cell(nRun,1); Sf = cell(nRun,1);
    dPix = false(1,size(stimData.SfullRaw{1},2));
    for r = 1:nRun
        Sf{r} = double(stimData.SfullRaw{r});
        Ss{r} = double(stimData.SscotRaw{r});
        Sd{r} = Sf{r} - Ss{r};
        dPix = dPix | any(abs(Sd{r}) > 1e-8,1);
    end

    % hrf{h}{r}: hemisphere h, run r.
    hrf = cell(nHem,1);
    for h = 1:nHem
        hrf{h} = cell(nRun,1);
        for r = 1:nRun
            hh = double(hrf_twogamma(stimData.hrfParams(h),stimData.tStim{r}));
            hrf{h}{r} = hh(:);
            if any(~isfinite(hrf{h}{r}))
                error('prfShiftFit:badHRF','HRF for hemisphere %d, run %d is not finite.',h,r);
            end
        end
    end

    % --- sampling limit --------------------------------------------------
    % A displacement finer than one pixel cannot be resolved, and a step
    % larger than sigma would evaluate the Gaussian at a negative width.
    dx = abs(diff(x,1,2)); dy = abs(diff(y,1,1));
    spacing = [dx(:);dy(:)];
    spacing = spacing(isfinite(spacing) & spacing > 0);
    if isempty(spacing), error('prfShiftFit:badGrid','Cannot determine pixel spacing.'); end
    pix = min(spacing);

    ecc   = sqrt(subjectData.prfXY(:,1).^2 + subjectData.prfXY(:,2).^2);
    sigma = subjectData.sigma(:);

    out = struct('ecc',ecc, 'sigma',sigma, ...
                 'WtW',nan(nVox,4,4), 'WtZ',nan(nVox,4), ...
                 'ZtZ',nan(nVox,1), 'dof',nan(nVox,1), 'ok',false(nVox,1), ...
                 'WtW2',nan(nVox,4,4), 'WtZ2',nan(nVox,4), ...
                 'ZtZ2',nan(nVox,1), 'dof2',nan(nVox,1), 'ok2',false(nVox,1), ...
                 'beta',nan(nVox,1), 'massIn',nan(nVox,1));

    for v = 1:nVox
        if opts.verbose && mod(v,500) == 0
            fprintf('    voxel %d of %d\n', v, nVox);
        end

        sg = sigma(v);
        if ecc(v) < opts.minEcc, continue, end
        if ~isfinite(sg) || sg < pix, continue, end   % pRF below one pixel

        x0 = subjectData.prfXY(v,1);
        y0 = subjectData.prfXY(v,2);

        % The radial / tangential frame is fixed by the ORIGINAL pRF centre
        % and never re-derived, so that an accumulated offset always means
        % displacement along the same axes.
        th = atan2(y0, x0);
        ur = [ cos(th), sin(th)];       % radial, outward
        ut = [-sin(th), cos(th)];       % tangential

        % linearisation point: the original pRF, displaced by any offset
        cx  = x0 + opts.offset(v,1)*ur(1) + opts.offset(v,2)*ut(1);
        cy  = y0 + opts.offset(v,1)*ur(2) + opts.offset(v,2)*ut(2);
        sgc = sg + opts.offset(v,3);
        if ~isfinite(sgc) || sgc < pix, continue, end

        hp = min(max(opts.hPos*sgc, 0.75*pix), 0.5*sgc);
        hs = min(max(opts.hSig*sgc, 0.75*pix), 0.5*sgc);

        % central differences of the AREA-NORMALISED Gaussian, so the
        % weights come out in degrees
        G  = normGauss(cx, cy, sgc, x, y);
        B  = [ (normGauss(cx+hp*ur(1), cy+hp*ur(2), sgc, x, y) - ...
                normGauss(cx-hp*ur(1), cy-hp*ur(2), sgc, x, y)) / (2*hp), ...
               (normGauss(cx+hp*ut(1), cy+hp*ut(2), sgc, x, y) - ...
                normGauss(cx-hp*ut(1), cy-hp*ut(2), sgc, x, y)) / (2*hp), ...
               (normGauss(cx, cy, sgc+hs, x, y) - ...
                normGauss(cx, cy, sgc-hs, x, y)) / (2*hs) ];

        % beta is the amplitude of the ORIGINAL full-field pRF, so it does
        % not drift as the linearisation point moves
        if any(opts.offset(v,:) ~= 0)
            G0 = normGauss(x0, y0, sg, x, y);
        else
            G0 = G;
        end
        out.massIn(v) = sum(G0(dPix));

        % Project onto the pRF first, then convolve that single time course
        % with THIS voxel's hemisphere HRF. Convolution along time commutes
        % with the spatial projection, so this matches convolving the whole
        % design, without needing one convolved design per hemisphere.
        hv = hrf{hemIdx(v)};
        P = []; D = []; K = []; F = [];
        for r = 1:nRun
            P  = [P;  centre(convHRF(Ss{r} * G, hv{r},TR))];   %#ok<AGROW>
            D  = [D;  centre(convHRF(Ss{r} * B, hv{r},TR))];   %#ok<AGROW>
            K  = [K;  centre(convHRF(Sd{r} * G0,hv{r},TR))];   %#ok<AGROW>
            F  = [F;  centre(convHRF(Sf{r} * G0,hv{r},TR))];   %#ok<AGROW>
        end
        Pf = F;   % beta predictor and the k2 regressor are the same quantity

        Ys = Yscot(:,v);
        Yf = Yfull(:,v);
        good = isfinite(Ys) & isfinite(Yf) & isfinite(P) & isfinite(K) & ...
               isfinite(F) & isfinite(Pf) & all(isfinite(D), 2);
        n = nnz(good);
        if n < 20, continue, end

        Ys = Ys(good); Yf = Yf(good);
        P = P(good); D = D(good,:); K = K(good); F = F(good); Pf = Pf(good);

        if sum(Pf.^2) < eps, continue, end

        % amplitude from the FULL-FIELD condition, where feedforward is
        % complete and the estimate is clean. Using it for the scotoma
        % condition assumes the gain is unchanged outside the scotoma.
        beta = (Pf.' * Yf) / (Pf.' * Pf);
        if ~isfinite(beta) || beta <= 0, continue, end

        % k: beta is fixed and K contains only the removed stimulus. Outside
        % the scotoma K approaches zero, so those voxels carry little information
        % about k even though they can still constrain the pRF terms.
        Z = Ys - beta*P;
        W = beta * [D, K];
        if all(isfinite(W(:))) && all(isfinite(Z))
            out.WtW(v,:,:) = W.' * W;
            out.WtZ(v,:) = W.' * Z;
            out.ZtZ(v) = Z.' * Z;
            out.dof(v) = max(n-nRun,1); % one mean removed per run
            out.ok(v) = true;
        end

        % k2: beta is fixed. Do not residualise P; outside the scotoma F=P,
        % so projecting P out would make the intended gain control vanish.
        Z2 = Ys - beta*P;
        W2 = beta * [D, F];
        if all(isfinite(W2(:))) && all(isfinite(Z2))
            out.WtW2(v,:,:) = W2.' * W2;
            out.WtZ2(v,:) = W2.' * Z2;
            out.ZtZ2(v) = Z2.' * Z2;
            out.dof2(v) = max(n-nRun,1);
            out.ok2(v) = true;
        end
        out.beta(v) = beta;
    end

    if opts.verbose
        fprintf('    fitted %d k voxels and %d k2 voxels of %d\n', ...
                nnz(out.ok), nnz(out.ok2), nVox);
    end
end

% =====================================================================
function g = normGauss(x0, y0, sg, x, y)
% Unit-area Gaussian on the stimulus grid, using the project's Gauss.m.
    pRF.center = [x0, y0];
    pRF.sig    = sg;
    pRF.ar     = 1;
    g = Gauss(pRF, x, y, 1);
    s = sum(g);
    if ~isfinite(s) || s <= 0
        g = zeros(size(g));
    else
        g = g / s;
    end
end

function Z = centre(Z)
    Z = Z - mean(Z, 1, 'omitnan');
end

function Z = catCentred(C)
% Concatenate runs, removing each run's own mean first.
    Z = [];
    for r = 1:numel(C)
        Zr = double(C{r});
        Z  = [Z; Zr - mean(Zr, 1, 'omitnan')];   %#ok<AGROW>
    end
end
