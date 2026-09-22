function out = inferCSSKFromLinearCurve(subjectData,stimData,trueN, ...
                                         eccEdges,targetLinearK,kGrid,opts)
% inferCSSKFromLinearCurve  Infer CSS k from an observed linear-model curve.
%
% Noiseless CSS responses are simulated over kGrid using real pRF geometries
% and subject HRFs. For every eccentricity bin, the function finds the CSS
% k whose simulated response produces targetLinearK when fitted incorrectly
% with n=1. It then simulates that bin-specific CSS k curve and verifies it
% by fitting both the linear and correct-CSS models.
%
% The simulation uses separate pRFs for the two models. The misspecified
% linear fit uses each voxel's original fitted Gaussian, Glinear. The CSS
% generator and correct CSS refit use a narrower Gaussian, Gcss, with
%
%   sigma(Gcss) = effectiveSizeScale*sqrt(trueN)*sigma(Glinear).
%
% Because a CSS point-response has effective width sigma(Gcss)/sqrt(trueN),
% effectiveSizeScale=1 preserves the original effective pRF size instead of
% accidentally enlarging it by 1/sqrt(trueN). Values below one explicitly
% test smaller effective CSS pRFs.
%
% The generative model is
%   full: (Afull*Gcss)^trueN
%   scot: [Ascot*Gcss + k*(Afull-Ascot)*Gcss]^trueN.
% Each fitting model estimates beta from its own full-field predictor before
% beta is held fixed for the scotoma fit.
%
% opts.TR               repetition time                 (default 1.2)
% opts.maxVoxPerBin     sampled voxels/subject/bin      (default 25)
% opts.minVox           minimum voxels/subject/bin      (default 5)
% opts.rngSeed          sampling seed                   (default 1)
% opts.tolK             CSS verification tolerance     (default 0.001)
% opts.effectiveSizeScale CSS effective size / linear size (default 1)
% opts.verbose                                             (default true)

if nargin < 7 || isempty(opts), opts = struct(); end
if ~isfield(opts,'TR'), opts.TR = 1.2; end
if ~isfield(opts,'maxVoxPerBin'), opts.maxVoxPerBin = 25; end
if ~isfield(opts,'minVox'), opts.minVox = 5; end
if ~isfield(opts,'rngSeed'), opts.rngSeed = 1; end
if ~isfield(opts,'tolK'), opts.tolK = 0.001; end
if ~isfield(opts,'effectiveSizeScale'), opts.effectiveSizeScale = 1; end
if ~isfield(opts,'verbose'), opts.verbose = true; end
if ~iscell(subjectData), subjectData = {subjectData}; end
if ~iscell(stimData) || numel(stimData) ~= numel(subjectData)
    error('inferCSSKFromLinearCurve:subjectCount', ...
          'stimData must contain one entry per subject.');
end
eccEdges = double(eccEdges(:).');
targetLinearK = double(targetLinearK(:));
kGrid = double(kGrid(:).');
nBin = numel(eccEdges)-1;
if numel(eccEdges) < 2 || any(~isfinite(eccEdges)) || any(diff(eccEdges) <= 0)
    error('inferCSSKFromLinearCurve:badEdges','eccEdges must increase strictly.');
end
if numel(targetLinearK) ~= nBin || any(isinf(targetLinearK))
    error('inferCSSKFromLinearCurve:badTarget', ...
          'targetLinearK must have one finite value or NaN per eccentricity bin.');
end
if ~isscalar(trueN) || ~isfinite(trueN) || trueN <= 0
    error('inferCSSKFromLinearCurve:badExponent','trueN must be positive.');
end
if numel(kGrid) < 2 || any(~isfinite(kGrid)) || any(diff(kGrid) <= 0) || ...
   kGrid(1) < 0 || kGrid(end) > 1
    error('inferCSSKFromLinearCurve:badGrid', ...
          'kGrid must increase strictly within [0,1].');
end
if ~isscalar(opts.TR) || opts.TR <= 0 || ~isfinite(opts.TR) || ...
   ~isscalar(opts.maxVoxPerBin) || opts.maxVoxPerBin < 1 || ...
   opts.maxVoxPerBin ~= round(opts.maxVoxPerBin) || ...
   ~isscalar(opts.minVox) || opts.minVox < 1 || opts.minVox ~= round(opts.minVox) || ...
   ~isscalar(opts.tolK) || opts.tolK <= 0 || ~isfinite(opts.tolK) || ...
   ~isscalar(opts.effectiveSizeScale) || ~isfinite(opts.effectiveSizeScale) || ...
   opts.effectiveSizeScale <= 0
    error('inferCSSKFromLinearCurve:badOptions','Invalid option value.');
end
if ~(islogical(opts.verbose) || isnumeric(opts.verbose)) || ...
   ~isscalar(opts.verbose) || ~ismember(opts.verbose,[0 1])
    error('inferCSSKFromLinearCurve:badVerbose','verbose must be true or false.');
end

nSub = numel(subjectData);
nRun = numel(subjectData{1}.Yfull);
[Afull,Ascot,stimTime,xGrid,yGrid] = loadRawStimuli(nRun,opts.TR);
oldRng = rng;
restoreRng = onCleanup(@() rng(oldRng)); %#ok<NASGU>
rng(opts.rngSeed);

P = cell(nSub,1);
nVox = zeros(nSub,nBin);
nSubpixel = zeros(nSub,1);
for s = 1:nSub
    S = subjectData{s}; D = stimData{s};
    required = {'Gprf','prfXY','sigma','Yfull','hemIdx'};
    if ~isstruct(S) || ~all(isfield(S,required)) || ...
       ~isstruct(D) || ~isfield(D,'hrfParams') || numel(S.Yfull) ~= nRun
        error('inferCSSKFromLinearCurve:badSubject','Subject %d is incomplete.',s);
    end
    eccAll = hypot(S.prfXY(:,1),S.prfXY(:,2));
    binAll = discretize(eccAll,eccEdges);
    keep = false(size(eccAll));
    for b = 1:nBin
        idx = find(binAll == b);
        if numel(idx) > opts.maxVoxPerBin
            idx = idx(randperm(numel(idx),opts.maxVoxPerBin));
        end
        keep(idx) = true;
        nVox(s,b) = numel(idx);
    end
    Glinear = double(S.Gprf(:,keep));
    gMass = sum(Glinear,1);
    if isempty(Glinear) || any(~isfinite(Glinear(:))) || any(~isfinite(gMass) | gMass <= 0)
        error('inferCSSKFromLinearCurve:badPRF','Invalid sampled pRFs for subject %d.',s);
    end
    Glinear = Glinear./gMass;
    xy = double(S.prfXY(keep,:));
    sigmaLinear = double(S.sigma(keep));
    [Gcss,sigmaCSS,nSubpixel(s)] = makeCSSGaussians(xGrid,yGrid,xy, ...
                                      sigmaLinear,trueN,opts.effectiveSizeScale);
    ecc = eccAll(keep);
    bin = binAll(keep);
    % One HRF per hemisphere; hemIdxKeep selects it for each sampled voxel.
    nHem = numel(D.hrfParams);
    hemIdxKeep = double(S.hemIdx(keep));
    hemIdxKeep = hemIdxKeep(:).';
    if any(hemIdxKeep < 1 | hemIdxKeep > nHem)
        error('inferCSSKFromLinearCurve:badHemIdx', ...
              'hemIdx indexes outside hrfParams for subject %d.',s);
    end
    hrf = cell(nHem,1);
    for h = 1:nHem
        hrf{h} = cell(nRun,1);
        for r = 1:nRun
            hh = double(hrf_twogamma(D.hrfParams(h),stimTime{r}));
            hrf{h}{r} = hh(:);
        end
    end
    driveFCSS = cell(nRun,1); driveSCSS = cell(nRun,1); driveDCSS = cell(nRun,1);
    fullTrue = cell(nRun,1); fullLinear = cell(nRun,1);
    scotLinear = cell(nRun,1); missingLinear = cell(nRun,1);
    for r = 1:nRun
        if size(Afull{r},2) ~= size(Glinear,1)
            error('inferCSSKFromLinearCurve:gridMismatch', ...
                  'Stimulus and pRF grids disagree for subject %d.',s);
        end
        driveFCSS{r} = max(Afull{r}*Gcss,0);
        driveSCSS{r} = max(Ascot{r}*Gcss,0);
        driveDCSS{r} = max((Afull{r}-Ascot{r})*Gcss,0);
        driveFLinear = max(Afull{r}*Glinear,0);
        driveSLinear = max(Ascot{r}*Glinear,0);
        driveDLinear = max((Afull{r}-Ascot{r})*Glinear,0);
        fullTrue{r} = cssPredict(driveFCSS{r},trueN,hrf,opts.TR,hemIdxKeep,r);
        fullLinear{r} = cssPredict(driveFLinear,1,hrf,opts.TR,hemIdxKeep,r);
        scotLinear{r} = cssPredict(driveSLinear,1,hrf,opts.TR,hemIdxKeep,r);
        missingLinear{r} = cssPredict(driveDLinear,1,hrf,opts.TR,hemIdxKeep,r);
    end
    betaLinear = fitBeta(fullLinear,fullTrue);
    betaCSS = fitBeta(fullTrue,fullTrue);
    valid = isfinite(betaLinear) & betaLinear > 0 & ...
            isfinite(betaCSS) & betaCSS > 0;
    P{s} = struct('bin',bin,'valid',valid(:),'driveSCSS',{driveSCSS}, ...
                  'driveDCSS',{driveDCSS},'hrf',{hrf},'hemIdx',hemIdxKeep(:),'fullTrue',{fullTrue}, ...
                  'scotLinear',{scotLinear},'missingLinear',{missingLinear}, ...
                  'betaLinear',betaLinear(:),'betaCSS',betaCSS(:), ...
                  'sigmaLinear',sigmaLinear(:),'sigmaCSS',sigmaCSS(:));
end
if any(nSubpixel > 0)
    warning('inferCSSKFromLinearCurve:subpixelPRF', ...
            ['%d sampled CSS Gaussians have sigma below the stimulus-pixel spacing. ' ...
             'Their sampled shapes may be grid-sensitive; inspect out.nSubpixelGaussian.'], ...
            sum(nSubpixel));
end

% Forward map: true CSS k -> k recovered by the misspecified linear model.
kLinearBySubject = nan(nSub,nBin,numel(kGrid));
for g = 1:numel(kGrid)
    for s = 1:nSub
        Yscot = simulateScot(P{s},kGrid(g),trueN,opts.TR);
        for b = 1:nBin
            q = P{s}.bin == b & P{s}.valid;
            if nnz(q) < opts.minVox, continue, end
            kLinearBySubject(s,b,g) = fitLinearK(P{s},Yscot,q);
        end
    end
    if opts.verbose
        fprintf('CSS exponent %.3g: simulated k grid %d/%d\n', ...
                trueN,g,numel(kGrid));
    end
end
linearMap = reshape(mean(kLinearBySubject,1,'omitnan'),nBin,numel(kGrid));

% Invert each bin's forward map at the observed linear-model k.
inferredCSSK = nan(nBin,1);
atLimit = false(nBin,1);
attainableMin = min(linearMap,[],2,'omitnan');
attainableMax = max(linearMap,[],2,'omitnan');
for b = 1:nBin
    if ~isfinite(targetLinearK(b)), continue, end
    x = linearMap(b,:);
    good = isfinite(x);
    if nnz(good) < 2, continue, end
    x = x(good); kg = kGrid(good); target = targetLinearK(b);
    exact = find(abs(x-target) <= 10*eps(max(1,abs(target))),1);
    crossing = find((x(1:end-1)-target).*(x(2:end)-target) < 0,1);
    if ~isempty(exact)
        inferredCSSK(b) = kg(exact);
    elseif ~isempty(crossing)
        inferredCSSK(b) = kg(crossing)+(target-x(crossing))* ...
            (kg(crossing+1)-kg(crossing))/(x(crossing+1)-x(crossing));
    else
        [~,ii] = min(abs(x-target));
        inferredCSSK(b) = kg(ii);
        atLimit(b) = target < min(x) || target > max(x);
        if ~atLimit(b)
            warning('inferCSSKFromLinearCurve:nonmonotonicMap', ...
                    ['No grid segment crossed the target in eccentricity bin %d ' ...
                     'although it lies within the map range; using the closest grid value.'],b);
        end
    end
end

% Simulate the inferred eccentricity-dependent CSS curve and refit it both
% ways. The linear refit should reproduce targetLinearK; the CSS refit
% should reproduce inferredCSSK.
verifiedLinearBySubject = nan(nSub,nBin);
recoveredCSSBySubject = nan(nSub,nBin);
for s = 1:nSub
    kVoxel = inferredCSSK(P{s}.bin);
    Yscot = simulateScot(P{s},kVoxel,trueN,opts.TR);
    for b = 1:nBin
        q = P{s}.bin == b & P{s}.valid;
        if nnz(q) < opts.minVox, continue, end
        verifiedLinearBySubject(s,b) = fitLinearK(P{s},Yscot,q);
        recoveredCSSBySubject(s,b) = fitCSSK(P{s},Yscot,q,trueN,opts.TR,opts.tolK);
    end
end

out = struct();
out.exponent = trueN;
out.effectiveSizeScale = opts.effectiveSizeScale;
out.cssGaussianScale = opts.effectiveSizeScale*sqrt(trueN);
out.eccEdges = eccEdges;
out.ecc = ((eccEdges(1:end-1)+eccEdges(2:end))/2).';
out.kGrid = kGrid;
out.targetLinearK = targetLinearK;
out.linearMap = linearMap;
out.kLinearBySubject = kLinearBySubject;
out.inferredCSSK = inferredCSSK;
out.verifiedLinearK = mean(verifiedLinearBySubject,1,'omitnan').';
out.recoveredCSSK = mean(recoveredCSSBySubject,1,'omitnan').';
out.verifiedLinearBySubject = verifiedLinearBySubject;
out.recoveredCSSBySubject = recoveredCSSBySubject;
out.attainableLinearMin = attainableMin;
out.attainableLinearMax = attainableMax;
out.atLimit = atLimit;
out.nVox = nVox;
out.nSubpixelGaussian = nSubpixel;
qCheck = ~atLimit & isfinite(targetLinearK) & isfinite(out.verifiedLinearK);
if any(qCheck & abs(out.verifiedLinearK-targetLinearK) > 0.02)
    warning('inferCSSKFromLinearCurve:inversionError', ...
            'A linear-model verification differs from its target by more than 0.02; refine kGrid.');
end
qCheck = isfinite(inferredCSSK) & isfinite(out.recoveredCSSK);
if any(qCheck & abs(out.recoveredCSSK-inferredCSSK) > 0.01)
    warning('inferCSSKFromLinearCurve:cssRecoveryError', ...
            'The correct CSS refit differs from the simulated k by more than 0.01.');
end
out.table = table(out.ecc,repmat(trueN,nBin,1), ...
                  repmat(opts.effectiveSizeScale,nBin,1), ...
                  repmat(out.cssGaussianScale,nBin,1),targetLinearK, ...
                  inferredCSSK,out.verifiedLinearK,out.recoveredCSSK, ...
                  attainableMin,attainableMax,atLimit,sum(nVox,1).', ...
    'VariableNames',{'ecc','exponent','effective_size_scale','css_gaussian_scale', ...
                     'target_linear_k','inferred_css_k', ...
                     'verified_linear_k','recovered_css_k', ...
                     'attainable_linear_min','attainable_linear_max', ...
                     'at_limit','nVox'});
end

function [Afull,Ascot,stimTime,xGrid,yGrid] = loadRawStimuli(nRun,TR)
Afull = cell(nRun,1); Ascot = cell(nRun,1); stimTime = cell(nRun,1);
xGrid = []; yGrid = [];
for r = 1:nRun
    [imgF,funcF] = loadScotomaStimuli(r,'logbar',TR);
    [imgS,funcS] = loadScotomaStimuli(r,'scotoma',TR);
    if ~isequal(size(imgF),size(imgS)) || ...
       ~isequal(size(funcF.x),size(funcS.x)) || ...
       ~isequal(size(funcF.t),size(funcS.t)) || ...
       max(abs(funcF.x(:)-funcS.x(:))) > 1e-10 || ...
       max(abs(funcF.y(:)-funcS.y(:))) > 1e-10 || ...
       max(abs(funcF.t(:)-funcS.t(:))) > 1e-10
        error('inferCSSKFromLinearCurve:stimulusMismatch', ...
              'Full and scotoma stimuli disagree in run %d.',r);
    end
    nt = numel(funcF.t); nPix = numel(funcF.x);
    Afull{r} = reshape(double(imgF),[nPix,nt]).';
    Ascot{r} = reshape(double(imgS),[nPix,nt]).';
    stimTime{r} = double(funcF.t(:).');
    if r == 1
        xGrid = double(funcF.x); yGrid = double(funcF.y);
    elseif ~isequal(size(xGrid),size(funcF.x)) || ...
           ~isequal(size(yGrid),size(funcF.y)) || ...
           max(abs(xGrid(:)-funcF.x(:))) > 1e-10 || ...
           max(abs(yGrid(:)-funcF.y(:))) > 1e-10
        error('inferCSSKFromLinearCurve:gridMismatch', ...
              'The stimulus grid differs across runs.');
    end
end
end

function [G,sigmaCSS,nSubpixel] = makeCSSGaussians(x,y,xy,sigmaLinear,n,sizeScale)
sigmaCSS = sigmaLinear(:)*sqrt(n)*sizeScale;
nVox = numel(sigmaCSS);
if size(xy,1) ~= nVox || size(xy,2) ~= 2 || any(~isfinite(xy(:))) || ...
   any(~isfinite(sigmaCSS) | sigmaCSS <= 0)
    error('inferCSSKFromLinearCurve:badCSSPRF','Invalid CSS pRF parameters.');
end
G = zeros(numel(x),nVox);
for v = 1:nVox
    pRF = struct('center',xy(v,:),'sig',sigmaCSS(v),'ar',1);
    G(:,v) = Gauss(pRF,x,y,true);
end
area = sum(G,1);
if any(~isfinite(area) | area <= 0)
    error('inferCSSKFromLinearCurve:badCSSPRF', ...
          'A sampled CSS Gaussian has zero or invalid mass; use a finer stimulus grid.');
end
G = G./area;
ux = unique(double(x(:))); uy = unique(double(y(:)));
dx = diff(sort(ux)); dy = diff(sort(uy));
dxy = [dx(dx > 0);dy(dy > 0)];
if isempty(dxy)
    nSubpixel = 0;
else
    pixelSpacing = min(dxy);
    nSubpixel = nnz(sigmaCSS < pixelSpacing);
end
end

function beta = fitBeta(P,Y)
nVox = size(P{1},2);
num = zeros(1,nVox); den = zeros(1,nVox);
for r = 1:numel(P)
    good = isfinite(P{r}) & isfinite(Y{r});
    Pr = P{r}; Yr = Y{r}; Pr(~good) = 0; Yr(~good) = 0;
    num = num+sum(Pr.*Yr,1); den = den+sum(Pr.^2,1);
end
beta = num./den;
beta(den <= eps | ~isfinite(beta)) = NaN;
end

function Y = simulateScot(P,k,trueN,TR)
Y = cell(numel(P.driveSCSS),1);
if isscalar(k)
    for r = 1:numel(Y)
        Y{r} = cssPredict(P.driveSCSS{r}+k*P.driveDCSS{r},trueN,P.hrf,TR,P.hemIdx,r);
    end
else
    k = k(:).';
    for r = 1:numel(Y)
        Y{r} = cssPredict(P.driveSCSS{r}+P.driveDCSS{r}.*k,trueN,P.hrf,TR,P.hemIdx,r);
    end
end
end

function k = fitLinearK(P,Y,q)
A = 0; c = 0;
b = P.betaLinear(q).';
for r = 1:numel(Y)
    Z = Y{r}(:,q)-P.scotLinear{r}(:,q).*b;
    W = P.missingLinear{r}(:,q).*b;
    good = isfinite(Z) & isfinite(W);
    A = A+sum(W(good).^2); c = c+sum(W(good).*Z(good));
end
if A <= eps || ~isfinite(A) || ~isfinite(c), k = NaN; return, end
k = min(max(c/A,0),1);
end

function k = fitCSSK(P,Y,q,n,TR,tolK)
b = P.betaCSS(q).';
objective = @(v) cssSSE(v,P,Y,q,b,n,TR);
grid = linspace(0,1,21);
sse = arrayfun(objective,grid);
[~,ii] = min(sse);
k = grid(ii);
if ii > 1 && ii < numel(grid)
    o = optimset('Display','off','TolX',tolK);
    [v,f] = fminbnd(objective,grid(ii-1),grid(ii+1),o);
    if f < sse(ii), k = v; end
end
end

function sse = cssSSE(k,P,Y,q,b,n,TR)
sse = 0;
for r = 1:numel(Y)
    pred = cssPredict(P.driveSCSS{r}(:,q)+k*P.driveDCSS{r}(:,q),n,P.hrf,TR,P.hemIdx(q),r);
    R = Y{r}(:,q)-pred.*b;
    R = R(isfinite(R));
    sse = sse+sum(R.^2);
end
if ~isfinite(sse), sse = Inf; end
end

function Y = cssPredict(drive,n,hrf,TR,hemIdx,r)
% hrf{h}{r} is hemisphere h on run r; the exponent precedes convolution.
Y = centre(convByHemisphere(max(drive,0).^n,hrf,TR,hemIdx,r));
end
