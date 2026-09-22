function out = fitSampledPRFCSS(subjectData,stimData,x,y,roi,opts)
% fitSampledPRFCSS  Fit linear and CSS pRFs with independent locations.
%
% Independent study-1 pRF files supply x and y. The three current full-field
% runs estimate sigma and beta for the linear model (n=1), and sigma, n, and
% beta for the CSS model. Those parameters are then held fixed while k is
% fitted only to the three scotoma runs. Group error bars are standard errors
% across subject-level estimates; voxels are never treated as subjects.

if nargin < 6 || isempty(opts), opts = struct(); end
opts = defaults(opts);
if ~iscell(subjectData), subjectData = {subjectData}; end
if ~iscell(stimData) || numel(stimData) ~= numel(subjectData)
    error('fitSampledPRFCSS:subjectCount','stimData must contain one entry per subject.');
end
if exist('fminsearchcon','file') ~= 2
    error('fitSampledPRFCSS:missingOptimizer','fminsearchcon must be on the MATLAB path.');
end
validateattributes(x,{'numeric'},{'real','finite','nonempty'},mfilename,'x');
validateattributes(y,{'numeric'},{'real','finite','size',size(x)},mfilename,'y');
validateattributes(roi,{'numeric'},{'scalar','integer','positive','finite'},mfilename,'roi');
nSub = numel(subjectData); nBin = numel(opts.eccEdges)-1;
kLinearBySubject = nan(nSub,nBin); kCSSBySubject = nan(nSub,nBin);
kLinearAtBoundary = false(nSub,nBin); kCSSAtBoundary = false(nSub,nBin);
nVoxBySubject = zeros(nSub,nBin); medianNBySubject = nan(nSub,nBin);
medianLinearSigmaBySubject = nan(nSub,nBin);
medianCssEffectiveSigmaBySubject = nan(nSub,nBin);
voxelTables = cell(nSub,1);
oldRng = rng; cleanupRng = onCleanup(@() rng(oldRng)); %#ok<NASGU>
rng(opts.rngSeed+1000*roi);
for s = 1:nSub
    S = subjectData{s}; D = stimData{s};
    validateSubject(S,D,x,y,s);
    nRun = numel(S.Yfull);
    eccAll = hypot(S.prfXY(:,1),S.prfXY(:,2));
    binAll = discretize(eccAll,opts.eccEdges);
    sampled = sampleByBin(binAll,nBin,opts.nPerBin);
    if isempty(sampled)
        warning('fitSampledPRFCSS:noVoxels','Subject %d has no voxels in range.',s);
        voxelTables{s} = table(); continue
    end
    sampledBin = binAll(sampled); nSample = numel(sampled);
    Yfull = selectColumns(S.Yfull,sampled); Yscot = selectColumns(S.Yscot,sampled);
    xy = double(S.prfXY(sampled,:)); hemIdx = double(S.hemIdx(sampled));
    hrf = makeHRFs(D,nRun);
    args = struct('Afull',{D.SfullRaw},'hrf',{hrf},'TR',D.TR,'x',double(x), ...
        'y',double(y),'sigmaBounds',opts.sigmaBounds,'nBounds',opts.nBounds, ...
        'nStarts',opts.cssStarts,'optim',opts.optim);
    [linearSigma,linearBeta,linearSSE,linearR2,linearCorr,linearExitflag, ...
     cssSigma,cssN,cssBeta,cssSSE,cssR2,cssCorr,cssExitflag] = deal(nan(nSample,1));
    useParfor = opts.useParallel && license('test','Distrib_Computing_Toolbox');
    if useParfor
        parfor v = 1:nSample
            yv = selectVoxel(Yfull,v);
            [linearSigma(v),linearBeta(v),linearSSE(v),linearR2(v),linearCorr(v), ...
             linearExitflag(v),cssSigma(v),cssN(v),cssBeta(v),cssSSE(v),cssR2(v), ...
             cssCorr(v),cssExitflag(v)] = fitVoxel(yv,xy(v,:), ...
                double(S.sigma(sampled(v))),args,hemIdx(v));
        end
    else
        for v = 1:nSample
            yv = selectVoxel(Yfull,v);
            [linearSigma(v),linearBeta(v),linearSSE(v),linearR2(v),linearCorr(v), ...
             linearExitflag(v),cssSigma(v),cssN(v),cssBeta(v),cssSSE(v),cssR2(v), ...
             cssCorr(v),cssExitflag(v)] = fitVoxel(yv,xy(v,:), ...
                double(S.sigma(sampled(v))),args,hemIdx(v));
            if opts.verbose && (mod(v,25) == 0 || v == nSample)
                fprintf('V%d subject %d/%d: fitted %d/%d voxels\n',roi,s,nSub,v,nSample);
            end
        end
    end
    validLinear = isfinite(linearSigma) & linearSigma > 0 & isfinite(linearBeta) & ...
        linearBeta > 0 & isfinite(linearR2) & linearR2 >= opts.minFullR2 & linearExitflag > 0;
    validCSS = isfinite(cssSigma) & cssSigma > 0 & isfinite(cssN) & ...
        isfinite(cssBeta) & cssBeta > 0 & isfinite(cssR2) & ...
        cssR2 >= opts.minFullR2 & cssExitflag > 0;
    validForK = validLinear & validCSS;
    cssEffectiveSigma = cssSigma./sqrt(cssN);
    linearSigmaAtBoundary = atBounds(linearSigma,opts.sigmaBounds,opts.boundaryTol);
    cssSigmaAtBoundary = atBounds(cssSigma,opts.sigmaBounds,opts.boundaryTol);
    cssNAtBoundary = atBounds(cssN,opts.nBounds,opts.boundaryTol);
    for b = 1:nBin
        q = sampledBin == b & validForK;
        nVoxBySubject(s,b) = nnz(q);
        if nnz(q) >= opts.minVoxPerBin
            [kLinearBySubject(s,b),kLinearAtBoundary(s,b)] = fitK(xy(q,:), ...
                linearSigma(q),ones(nnz(q),1),linearBeta(q),selectVoxels(Yscot,q), ...
                D,hrf,x,y,hemIdx(q),opts);
            [kCSSBySubject(s,b),kCSSAtBoundary(s,b)] = fitK(xy(q,:),cssSigma(q), ...
                cssN(q),cssBeta(q),selectVoxels(Yscot,q),D,hrf,x,y,hemIdx(q),opts);
        end
        qCss = sampledBin == b & validCSS;
        if nnz(qCss) >= opts.minVoxPerBin
            medianNBySubject(s,b) = median(cssN(qCss),'omitnan');
        end
        if nnz(q) >= opts.minVoxPerBin
            medianLinearSigmaBySubject(s,b) = median(linearSigma(q),'omitnan');
            medianCssEffectiveSigmaBySubject(s,b) = median(cssEffectiveSigma(q),'omitnan');
        end
    end
    subNum = s; if isfield(S,'subNum'), subNum = S.subNum; end
    sourceVoxelIndex = sampled;
    if isfield(S,'sourceVoxelIndex'), sourceVoxelIndex = S.sourceVoxelIndex(sampled); end
    voxelTables{s} = table(repmat(subNum,nSample,1),repmat(roi,nSample,1), ...
        reshape(sampled,[],1),reshape(sourceVoxelIndex,[],1),reshape(sampledBin,[],1), ...
        reshape(eccAll(sampled),[],1),xy(:,1),xy(:,2),reshape(double(S.sigma(sampled)),[],1), ...
        reshape(double(S.w_vox(sampled)),[],1),linearSigma,linearBeta,linearSSE,linearR2, ...
        linearCorr,linearExitflag,validLinear,linearSigmaAtBoundary,cssSigma,cssN, ...
        cssEffectiveSigma,cssBeta,cssSSE,cssR2,cssCorr,cssExitflag,validCSS, ...
        cssSigmaAtBoundary,cssNAtBoundary,validForK, ...
        'VariableNames',{'subNum','ROI','voxelIndex','sourceVoxelIndex','eccBin', ...
        'originalEcc','independentX','independentY','originalSigma','originalVexpl', ...
        'linearSigma','linearBeta','linearSSE','linearR2','linearCorr','linearExitflag', ...
        'validLinearPRF','linearSigmaAtBoundary','cssSigma','cssN','cssEffectiveSigma', ...
        'cssBeta','cssSSE','cssR2','cssCorr','cssExitflag','validCSSPRF', ...
        'cssSigmaAtBoundary','cssNAtBoundary','validForK'});
end
hasRows = cellfun(@(T) istable(T) && height(T) > 0,voxelTables);
if any(hasRows), voxelTable = vertcat(voxelTables{hasRows}); else, voxelTable = table(); end
summaryTable = summarize(kLinearBySubject,kCSSBySubject,nVoxBySubject, ...
    medianNBySubject,medianLinearSigmaBySubject,medianCssEffectiveSigmaBySubject, ...
    kLinearAtBoundary,kCSSAtBoundary,roi,opts);
out = struct('ROI',roi,'options',opts,'voxelTable',voxelTable,'summaryTable',summaryTable, ...
    'kLinearBySubject',kLinearBySubject,'kCSSBySubject',kCSSBySubject, ...
    'kLinearAtBoundary',kLinearAtBoundary,'kCSSAtBoundary',kCSSAtBoundary, ...
    'nVoxBySubject',nVoxBySubject,'medianNBySubject',medianNBySubject, ...
    'medianLinearSigmaBySubject',medianLinearSigmaBySubject, ...
    'medianCssEffectiveSigmaBySubject',medianCssEffectiveSigmaBySubject);
end

function opts = defaults(opts)
if ~isfield(opts,'eccEdges'), opts.eccEdges = 0:0.25:5; end
if ~isfield(opts,'nPerBin'), opts.nPerBin = Inf; end
if ~isfield(opts,'minVoxPerBin'), opts.minVoxPerBin = 5; end
if ~isfield(opts,'rngSeed'), opts.rngSeed = 1; end
if ~isfield(opts,'sigmaBounds'), opts.sigmaBounds = [0.05 8]; end
if ~isfield(opts,'nBounds'), opts.nBounds = [0.05 1.35]; end
if ~isfield(opts,'cssStarts'), opts.cssStarts = [0.12 0.33 0.75 1 1.2]; end
if ~isfield(opts,'minFullR2'), opts.minFullR2 = -Inf; end
if ~isfield(opts,'kBounds'), opts.kBounds = [0 1]; end
if ~isfield(opts,'kGridStep'), opts.kGridStep = 0.05; end
if ~isfield(opts,'boundaryTol'), opts.boundaryTol = 0.01; end
if ~isfield(opts,'useParallel'), opts.useParallel = true; end
if ~isfield(opts,'verbose'), opts.verbose = true; end
if ~isfield(opts,'optim')
    opts.optim = optimset('Display','off','MaxIter',400,'MaxFunEvals',2000, ...
        'TolX',1e-4,'TolFun',1e-7);
end
opts.eccEdges = double(opts.eccEdges(:).'); opts.sigmaBounds = double(opts.sigmaBounds(:).');
opts.nBounds = double(opts.nBounds(:).'); opts.kBounds = double(opts.kBounds(:).');
opts.cssStarts = double(opts.cssStarts(:).');
if numel(opts.eccEdges) < 2 || opts.eccEdges(1) < 0 || any(~isfinite(opts.eccEdges)) || ...
        any(diff(opts.eccEdges) <= 0)
    error('fitSampledPRFCSS:badEccEdges','eccEdges must increase from a nonnegative value.');
end
if numel(opts.sigmaBounds) ~= 2 || opts.sigmaBounds(1) <= 0 || ...
        opts.sigmaBounds(2) <= opts.sigmaBounds(1)
    error('fitSampledPRFCSS:badSigmaBounds','sigmaBounds must be increasing and positive.');
end
if numel(opts.nBounds) ~= 2 || opts.nBounds(1) <= 0 || opts.nBounds(2) <= opts.nBounds(1) || ...
        any(opts.cssStarts < opts.nBounds(1) | opts.cssStarts > opts.nBounds(2))
    error('fitSampledPRFCSS:badNBounds','n bounds and starting values are inconsistent.');
end
if numel(opts.kBounds) ~= 2 || opts.kBounds(1) < 0 || opts.kBounds(2) > 1 || ...
        opts.kBounds(2) <= opts.kBounds(1)
    error('fitSampledPRFCSS:badKBounds','kBounds must increase within [0,1].');
end
end

function validateSubject(S,D,x,y,s)
requiredS = {'Yfull','Yscot','prfXY','sigma','w_vox','hemIdx'};
requiredD = {'SfullRaw','SscotRaw','hrfParams','tStim','TR'};
if ~isstruct(S) || ~all(isfield(S,requiredS)) || ~isstruct(D) || ~all(isfield(D,requiredD))
    error('fitSampledPRFCSS:badSubject','Subject %d is incomplete.',s);
end
nRun = numel(S.Yfull); nVox = size(S.prfXY,1);
if nRun == 0 || numel(S.Yscot) ~= nRun || numel(D.SfullRaw) ~= nRun || ...
        numel(D.SscotRaw) ~= nRun || numel(D.tStim) ~= nRun || size(S.prfXY,2) ~= 2 || ...
        numel(S.sigma) ~= nVox || numel(S.w_vox) ~= nVox || numel(S.hemIdx) ~= nVox
    error('fitSampledPRFCSS:dimensionMismatch','Subject %d has inconsistent dimensions.',s);
end
for r = 1:nRun
    if ~isequal(size(S.Yfull{r}),size(S.Yscot{r})) || size(S.Yfull{r},2) ~= nVox || ...
            ~isequal(size(D.SfullRaw{r}),size(D.SscotRaw{r})) || ...
            size(D.SfullRaw{r},1) ~= size(S.Yfull{r},1) || size(D.SfullRaw{r},2) ~= numel(x)
        error('fitSampledPRFCSS:runDimensionMismatch', ...
              'Subject %d, run %d has mismatched BOLD/stimulus dimensions.',s,r);
    end
end
if ~isequal(size(x),size(y)) || any(~isfinite(S.prfXY(:))) || ...
        any(~isfinite(S.sigma(:)) | S.sigma(:) <= 0) || ...
        any(double(S.hemIdx(:)) < 1 | double(S.hemIdx(:)) > numel(D.hrfParams))
    error('fitSampledPRFCSS:invalidValues','Subject %d has invalid pRF or hemisphere values.',s);
end
end

function sampled = sampleByBin(bin,nBin,nPerBin)
sampled = [];
for b = 1:nBin
    idx = find(bin == b);
    if isfinite(nPerBin) && numel(idx) > nPerBin, idx = idx(randperm(numel(idx),nPerBin)); end
    sampled = [sampled;idx(:)]; %#ok<AGROW>
end
end

function hrf = makeHRFs(D,nRun)
nHem = numel(D.hrfParams); hrf = cell(nHem,1);
for h = 1:nHem
    hrf{h} = cell(nRun,1);
    for r = 1:nRun
        hh = double(hrf_twogamma(D.hrfParams(h),D.tStim{r}));
        if isempty(hh) || any(~isfinite(hh)), error('fitSampledPRFCSS:badHRF','Invalid HRF.'); end
        hrf{h}{r} = hh(:);
    end
end
end

function [sLin,bLin,sseLin,r2Lin,corrLin,exitLin,sCss,nCss,bCss,sseCss,r2Css,corrCss,exitCss] = ...
    fitVoxel(Y,xy,sigma0,args,hem)
lb = args.sigmaBounds(1); ub = args.sigmaBounds(2);
startSigma = min(max(sigma0,lb+1e-6),ub-1e-6);
[pLin,~,exitLin] = fminsearchcon(@(p) fullObjective(xy,p(1),1,Y,args,hem), ...
    startSigma,lb,ub,[],[],[],args.optim);
sLin = pLin(1); [sseLin,bLin,r2Lin,corrLin] = fullObjective(xy,sLin,1,Y,args,hem);
best = Inf; [sCss,nCss,bCss,sseCss,r2Css,corrCss,exitCss] = deal(NaN);
for n0 = args.nStarts
    start = [min(max(sLin*sqrt(n0),lb+1e-6),ub-1e-6),n0];
    [candidate,f,flag] = fminsearchcon(@(p) fullObjective(xy,p(1),p(2),Y,args,hem), ...
        start,[lb args.nBounds(1)],[ub args.nBounds(2)],[],[],[],args.optim);
    if isfinite(f) && f < best
        best = f; sCss = candidate(1); nCss = candidate(2); exitCss = flag;
    end
end
if isfinite(sCss), [sseCss,bCss,r2Css,corrCss] = fullObjective(xy,sCss,nCss,Y,args,hem); end
end

function [sse,beta,r2,rModel] = fullObjective(xy,sigma,n,Y,args,hem)
G = gaussian(xy,sigma,args.x,args.y); P = cell(numel(Y),1);
for r = 1:numel(Y)
    P{r} = cssPredict(nonnegative(args.Afull{r}*G),n,args.hrf,args.TR,hem,r);
end
[beta,sse,r2,rModel] = fitAmplitude(P,Y);
if ~isfinite(sse), sse = realmax/1e100; end
end

function [beta,sse,r2,rModel] = fitAmplitude(P,Y)
num = 0; den = 0; sst = 0; allP = []; allY = [];
for r = 1:numel(P)
    good = isfinite(P{r}) & isfinite(Y{r}); pr = P{r}(good); yr = Y{r}(good);
    num = num+sum(pr.*yr); den = den+sum(pr.^2); sst = sst+sum(yr.^2);
end
if den <= eps || sst <= eps, beta = NaN; sse = Inf; r2 = NaN; rModel = NaN; return, end
beta = num/den;
if ~isfinite(beta) || beta <= 0, beta = NaN; sse = Inf; r2 = NaN; rModel = NaN; return, end
sse = 0;
for r = 1:numel(P)
    good = isfinite(P{r}) & isfinite(Y{r}); pr = beta*P{r}(good); yr = Y{r}(good);
    sse = sse+sum((yr-pr).^2); allP = [allP;pr(:)]; allY = [allY;yr(:)]; %#ok<AGROW>
end
r2 = 1-sse/sst;
if numel(allP) > 1 && std(allP) > 0 && std(allY) > 0
    C = corrcoef(allP,allY); rModel = C(1,2);
else
    rModel = NaN;
end
end

function [k,atBoundary] = fitK(xy,sigma,n,beta,Y,D,hrf,x,y,hemIdx,opts)
G = gaussianMatrix(xy,sigma,x,y); driveS = cell(numel(Y),1); driveD = cell(numel(Y),1);
for r = 1:numel(Y)
    driveS{r} = nonnegative(double(D.SscotRaw{r})*G);
    driveD{r} = nonnegative((double(D.SfullRaw{r})-double(D.SscotRaw{r}))*G);
end
objective = @(v) kObjective(v,driveS,driveD,n,beta,Y,hrf,D.TR,hemIdx);
grid = opts.kBounds(1):opts.kGridStep:opts.kBounds(2);
if grid(end) < opts.kBounds(2), grid = [grid opts.kBounds(2)]; end
sse = arrayfun(objective,grid); [best,ii] = min(sse); k = grid(ii);
if ii > 1 && ii < numel(grid)
    [candidate,f] = fminbnd(objective,grid(ii-1),grid(ii+1),optimset('Display','off','TolX',1e-4));
    if f < best, k = candidate; best = f; end
end
if ~isfinite(best), k = NaN; end
atBoundary = isfinite(k) && min(abs(k-opts.kBounds)) <= opts.boundaryTol;
end

function sse = kObjective(k,driveS,driveD,n,beta,Y,hrf,TR,hemIdx)
sse = 0; nGood = 0;
for r = 1:numel(Y)
    P = cssPredict(driveS{r}+k*driveD{r},n,hrf,TR,hemIdx,r);
    R = Y{r}-P.*beta(:).'; good = isfinite(R);
    sse = sse+sum(R(good).^2); nGood = nGood+nnz(good);
end
if nGood == 0 || ~isfinite(sse), sse = Inf; end
end

function T = summarize(kLin,kCss,nVox,medianN,medianLinSigma,medianCssEff,linBoundary,cssBoundary,roi,opts)
nBin = size(kLin,2); ecc = ((opts.eccEdges(1:end-1)+opts.eccEdges(2:end))/2).';
[meanLin,seLin,meanCss,seCss,delta,seDelta,meanN,seN,meanLinSigma,seLinSigma, ...
 meanCssEff,seCssEff,meanRatio,seRatio] = deal(nan(nBin,1));
[nSubjects,nVoxTotal,nLinearBoundary,nCssBoundary] = deal(zeros(nBin,1));
for b = 1:nBin
    q = isfinite(kLin(:,b)) & isfinite(kCss(:,b));
    nSubjects(b) = nnz(q); nVoxTotal(b) = sum(nVox(q,b));
    [meanLin(b),seLin(b)] = meanAndSE(kLin(q,b));
    [meanCss(b),seCss(b)] = meanAndSE(kCss(q,b));
    [delta(b),seDelta(b)] = meanAndSE(kCss(q,b)-kLin(q,b));
    nLinearBoundary(b) = nnz(linBoundary(q,b)); nCssBoundary(b) = nnz(cssBoundary(q,b));
    [meanN(b),seN(b)] = meanAndSE(medianN(:,b));
    [meanLinSigma(b),seLinSigma(b)] = meanAndSE(medianLinSigma(:,b));
    [meanCssEff(b),seCssEff(b)] = meanAndSE(medianCssEff(:,b));
    [meanRatio(b),seRatio(b)] = meanAndSE(medianCssEff(:,b)./medianLinSigma(:,b));
end
T = table(repmat(roi,nBin,1),ecc,nSubjects,nVoxTotal,meanLin,seLin,meanCss,seCss, ...
    delta,seDelta,meanN,seN,meanLinSigma,seLinSigma,meanCssEff,seCssEff,meanRatio,seRatio, ...
    nLinearBoundary,nCssBoundary,'VariableNames',{'ROI','ecc','nSubjects','nVox', ...
    'kLinear','kLinear_se','kCSS','kCSS_se','deltaK','deltaK_se','medianCssN', ...
    'medianCssN_se','medianLinearSigma','medianLinearSigma_se','medianCssEffectiveSigma', ...
    'medianCssEffectiveSigma_se','medianEffectiveSizeRatio','medianEffectiveSizeRatio_se', ...
    'nLinearKAtBoundary','nCSSKAtBoundary'});
end

function q = atBounds(v,bounds,tol)
q = isfinite(v) & (abs(v-bounds(1)) <= tol | abs(v-bounds(2)) <= tol);
end

function G = gaussian(xy,sigma,x,y)
if numel(xy) ~= 2 || ~isfinite(sigma) || sigma <= 0, G = nan(numel(x),1); return, end
z = -((double(x)-xy(1)).^2+(double(y)-xy(2)).^2)/(2*sigma^2); z = z-max(z(:));
G = exp(z(:)); mass = sum(G);
if ~isfinite(mass) || mass <= 0, G(:) = NaN; else, G = G/mass; end
end

function G = gaussianMatrix(xy,sigma,x,y)
G = nan(numel(x),size(xy,1));
for v = 1:size(xy,1), G(:,v) = gaussian(xy(v,:),sigma(v),x,y); end
end

function drive = nonnegative(drive)
drive(isfinite(drive) & drive < 0) = 0;
end

function Y = cssPredict(drive,n,hrf,TR,hemIdx,r)
Y = centre(convByHemisphere(bsxfun(@power,nonnegative(drive),double(n(:).')), ...
                            hrf,TR,hemIdx,r));
end

function Y = selectColumns(C,idx)
Y = cell(size(C)); for r = 1:numel(C), Y{r} = centre(double(C{r}(:,idx))); end
end

function Y = selectVoxel(C,v)
Y = cell(size(C)); for r = 1:numel(C), Y{r} = C{r}(:,v); end
end

function Y = selectVoxels(C,q)
Y = cell(size(C)); for r = 1:numel(C), Y{r} = C{r}(:,q); end
end

function Z = centre(Z)
Z = Z-mean(Z,1,'omitnan');
end

function [m,se] = meanAndSE(v)
v = v(isfinite(v));
if isempty(v), m = NaN; se = NaN; return, end
m = mean(v); if numel(v) < 2, se = NaN; else, se = std(v,0)/sqrt(numel(v)); end
end
