function out = fitSampledPRFCSS(subjectData,stimData,x,y,roi,opts)
% fitSampledPRFCSS  Refit sampled voxels with linear and CSS pRF models.
%
% Both pRF models are fitted only to the three independent mapping runs
% stored in Yprf. Their fitted geometry is then held fixed. Beta is
% re-estimated from the separate three full-field comparison runs, and k is
% fitted to the three scotoma runs. The same voxels are used throughout.
%
% The CSS model is
%   neural(t) = (stimulus(t)*G(x0,y0,sigma))^n
% followed by HRF convolution. The scotoma model replaces stimulus with
%   Sscot + k*(Sfull-Sscot).
%
% Voxel handling:
%   * compileStimAndSubData applies the original-pRF inclusion criteria.
%   * Both mapping models must fit successfully before a voxel contributes
%     to the direct linear-versus-CSS k comparison.
%   * CSS n diagnostics use successful CSS mapping fits and therefore do
%     not depend on the separate full/scotoma comparison runs.
%   * Bounds are flagged in the output, not silently discarded.
%
% Required external functions: loadScotomaStimuli, hrf_twogamma, and
% fminsearchcon (included in pRF-master/matlab/externals).

if nargin < 6 || isempty(opts), opts = struct(); end
opts = defaults(opts);
if ~iscell(subjectData), subjectData = {subjectData}; end
if ~iscell(stimData) || numel(stimData) ~= numel(subjectData)
    error('fitSampledPRFCSS:subjectCount','stimData must contain one entry per subject.');
end
if exist('fminsearchcon','file') ~= 2
    error('fitSampledPRFCSS:missingOptimizer', ...
          'Add pRF-master and its subfolders to the MATLAB path so fminsearchcon is available.');
end
validateattributes(x,{'numeric'},{'real','finite','nonempty'},mfilename,'x');
validateattributes(y,{'numeric'},{'real','finite','size',size(x)},mfilename,'y');
validateattributes(roi,{'numeric'},{'real','finite','scalar','integer','>=',1},mfilename,'roi');
nSub = numel(subjectData);
nRun = numel(subjectData{1}.Yfull);
[Afull,Ascot,stimTime] = loadRawStimuli(nRun,opts.TR,x,y);
oldRng = rng;
cleanupRng = onCleanup(@() rng(oldRng)); %#ok<NASGU>
rng(opts.rngSeed+1000*roi);
nBin = numel(opts.eccEdges)-1;
kLinearBySubject = nan(nSub,nBin);
kCSSBySubject = nan(nSub,nBin);
nVoxBySubject = zeros(nSub,nBin);
medianNBySubject = nan(nSub,nBin);
medianLinearSigmaBySubject = nan(nSub,nBin);
medianCssEffectiveSigmaBySubject = nan(nSub,nBin);
voxelTables = cell(nSub,1);
for s = 1:nSub
    S = subjectData{s}; D = stimData{s};
    required = {'Yprf','Yfull','Yscot','prfXY','sigma','w_vox','hemIdx'};
    if ~isstruct(S) || ~all(isfield(S,required)) || ~isstruct(D) || ...
            ~all(isfield(D,{'hrfParams','Aprf','tPrf'}))
        error('fitSampledPRFCSS:badSubject','Subject %d is incomplete.',s);
    end
    if numel(D.hrfParams) < 1 || numel(D.hrfParams) > 2
        error('fitSampledPRFCSS:badHRF', ...
              'stimData{%d}.hrfParams must hold one struct per hemisphere.',s);
    end
    if any(double(S.hemIdx(:)) < 1 | double(S.hemIdx(:)) > numel(D.hrfParams))
        error('fitSampledPRFCSS:badHemIdx', ...
              'hemIdx for subject %d indexes outside hrfParams.',s);
    end
    if numel(S.Yprf) ~= nRun || numel(S.Yfull) ~= nRun || numel(S.Yscot) ~= nRun
        error('fitSampledPRFCSS:runCount','Run counts disagree for subject %d.',s);
    end
    nAll = size(S.prfXY,1);
    if numel(S.sigma) ~= nAll || numel(S.w_vox) ~= nAll || numel(S.hemIdx) ~= nAll
        error('fitSampledPRFCSS:voxelCount','pRF fields disagree for subject %d.',s);
    end
    for r = 1:nRun
        if size(S.Yprf{r},2) ~= nAll || size(S.Yfull{r},2) ~= nAll || ...
                ~isequal(size(S.Yfull{r}),size(S.Yscot{r})) || ...
                size(S.Yprf{r},1) ~= size(D.Aprf{r},1)
            error('fitSampledPRFCSS:dataSize','BOLD dimensions disagree for subject %d.',s);
        end
    end
    eccAll = hypot(S.prfXY(:,1),S.prfXY(:,2));
    binAll = discretize(eccAll,opts.eccEdges);
    sampled = [];
    sampledBin = [];
    for b = 1:nBin
        candidates = find(binAll == b);
        if numel(candidates) > opts.nPerBin
            candidates = candidates(randperm(numel(candidates),opts.nPerBin));
        end
        sampled = [sampled;candidates(:)]; %#ok<AGROW>
        sampledBin = [sampledBin;repmat(b,numel(candidates),1)]; %#ok<AGROW>
    end
    if isempty(sampled)
        warning('fitSampledPRFCSS:noVoxels','Subject %d has no voxels in the requested bins.',s);
        voxelTables{s} = table();
        continue
    end
    % One HRF per hemisphere, so hrf{h}{r} and hrfPrf{h}{r}. hemIdx selects
    % the row of D.hrfParams that belongs to each sampled voxel.
    nHem = numel(D.hrfParams);
    hrf = cell(nHem,1);
    hrfPrf = cell(nHem,1);
    Yprf = cell(nRun,1);
    Yfull = cell(nRun,1);
    Yscot = cell(nRun,1);
    for h = 1:nHem
        hrf{h} = cell(nRun,1);
        hrfPrf{h} = cell(nRun,1);
        for r = 1:nRun
            hrf{h}{r} = double(hrf_twogamma(D.hrfParams(h),stimTime{r}));
            hrf{h}{r} = hrf{h}{r}(:);
            hrfPrf{h}{r} = double(hrf_twogamma(D.hrfParams(h),D.tPrf{r}));
            hrfPrf{h}{r} = hrfPrf{h}{r}(:);
        end
    end
    for r = 1:nRun
        Yprf{r} = centre(double(S.Yprf{r}(:,sampled)));
        Yfull{r} = centre(double(S.Yfull{r}(:,sampled)));
        Yscot{r} = centre(double(S.Yscot{r}(:,sampled)));
    end
    hemIdx = double(S.hemIdx(sampled));
    hemIdx = hemIdx(:);
    nSample = numel(sampled);
    pLinear = nan(nSample,3); betaLinear = nan(nSample,1);
    sseLinear = nan(nSample,1); r2Linear = nan(nSample,1); corrLinear = nan(nSample,1);
    pCSS = nan(nSample,4); betaCSS = nan(nSample,1);
    sseCSS = nan(nSample,1); r2CSS = nan(nSample,1); corrCSS = nan(nSample,1);
    p0 = [double(S.prfXY(sampled,:)),double(S.sigma(sampled))];
    fitArgs = struct('Afull',{D.Aprf},'hrf',{hrfPrf},'TR',opts.TR,'x',double(x), ...
                     'y',double(y),'centreLimit',opts.centreLimit, ...
                     'sigmaBounds',opts.sigmaBounds,'nBounds',opts.nBounds, ...
                     'cssStarts',opts.cssStarts,'optim',opts.optim);
    useParfor = opts.useParallel && license('test','Distrib_Computing_Toolbox');
    if useParfor
        parfor v = 1:nSample
            yv = selectVoxel(Yprf,v);
            [pLinear(v,:),betaLinear(v),sseLinear(v),r2Linear(v),corrLinear(v), ...
             pCSS(v,:),betaCSS(v),sseCSS(v),r2CSS(v),corrCSS(v)] = ...
                fitVoxelModels(yv,p0(v,:),fitArgs,hemIdx(v));
        end
    else
        for v = 1:nSample
            yv = selectVoxel(Yprf,v);
            [pLinear(v,:),betaLinear(v),sseLinear(v),r2Linear(v),corrLinear(v), ...
             pCSS(v,:),betaCSS(v),sseCSS(v),r2CSS(v),corrCSS(v)] = ...
                fitVoxelModels(yv,p0(v,:),fitArgs,hemIdx(v));
            if opts.verbose && (mod(v,10) == 0 || v == nSample)
                fprintf('V%d subject %d/%d: fitted %d/%d sampled voxels\n', ...
                        roi,s,nSub,v,nSample);
            end
        end
    end
    mappingBetaLinear = betaLinear;
    mappingBetaCSS = betaCSS;
    betaLinear = refitBeta(pLinear,ones(nSample,1),Yfull,Afull,hrf,x,y,opts.TR,hemIdx);
    betaCSS = refitBeta(pCSS(:,1:3),pCSS(:,4),Yfull,Afull,hrf,x,y,opts.TR,hemIdx);
    validLinearPRF = all(isfinite(pLinear),2) & isfinite(mappingBetaLinear) & ...
                     mappingBetaLinear > 0 & isfinite(r2Linear) & ...
                     r2Linear >= opts.minMappingR2;
    validCSSPRF = all(isfinite(pCSS),2) & isfinite(mappingBetaCSS) & ...
                  mappingBetaCSS > 0 & isfinite(r2CSS) & ...
                  r2CSS >= opts.minMappingR2;
    validForK = validLinearPRF & validCSSPRF & ...
                isfinite(betaLinear) & betaLinear > 0 & ...
                isfinite(betaCSS) & betaCSS > 0;
    cssEffectiveSigma = pCSS(:,3)./sqrt(pCSS(:,4));
    for b = 1:nBin
        q = sampledBin == b & validForK;
        nVoxBySubject(s,b) = nnz(q);
        if nnz(q) < opts.minVoxPerBin, continue, end
        kLinearBySubject(s,b) = fitK(pLinear(q,:),ones(nnz(q),1),betaLinear(q), ...
                                      selectVoxels(Yscot,q),Afull,Ascot,hrf,x,y,opts,hemIdx(q));
        kCSSBySubject(s,b) = fitK(pCSS(q,1:3),pCSS(q,4),betaCSS(q), ...
                                  selectVoxels(Yscot,q),Afull,Ascot,hrf,x,y,opts,hemIdx(q));
    end
    for b = 1:nBin
        qCss = sampledBin == b & validCSSPRF;
        qBoth = sampledBin == b & validLinearPRF & validCSSPRF;
        if nnz(qCss) >= opts.minVoxPerBin
            medianNBySubject(s,b) = median(pCSS(qCss,4),'omitnan');
        end
        if nnz(qBoth) >= opts.minVoxPerBin
            medianLinearSigmaBySubject(s,b) = median(pLinear(qBoth,3),'omitnan');
            medianCssEffectiveSigmaBySubject(s,b) = median(cssEffectiveSigma(qBoth),'omitnan');
        end
    end
    subNum = s;
    if isfield(S,'subNum'), subNum = S.subNum; end
    nAtBoundary = abs(pCSS(:,4)-opts.nBounds(1)) <= opts.boundaryTol | ...
                  abs(pCSS(:,4)-opts.nBounds(2)) <= opts.boundaryTol;
    % MATLAB preserves the orientation of some indexed vectors. In
    % particular, w_vox is often stored as 1-by-n, which previously made
    % that table variable a row while every other variable had nSample
    % rows. Force every scalar-per-voxel field to an nSample-by-1 column.
    subColumn = repmat(subNum,nSample,1);
    roiColumn = repmat(roi,nSample,1);
    voxelIndex = reshape(sampled,[],1);
    if isfield(S,'sourceVoxelIndex')
        sourceVoxelIndex = reshape(S.sourceVoxelIndex(sampled),[],1);
    else
        sourceVoxelIndex = voxelIndex;
    end
    eccBin = reshape(sampledBin,[],1);
    originalEcc = reshape(eccAll(sampled),[],1);
    originalX = reshape(double(S.prfXY(sampled,1)),[],1);
    originalY = reshape(double(S.prfXY(sampled,2)),[],1);
    originalSigma = reshape(double(S.sigma(sampled)),[],1);
    originalVexpl = reshape(double(S.w_vox(sampled)),[],1);
    voxelTables{s} = table(subColumn,roiColumn,voxelIndex,sourceVoxelIndex,eccBin,originalEcc, ...
        originalX,originalY,originalSigma,originalVexpl, ...
        pLinear(:,1),pLinear(:,2),pLinear(:,3),mappingBetaLinear,betaLinear, ...
        sseLinear,r2Linear,corrLinear,validLinearPRF, ...
        pCSS(:,1),pCSS(:,2),pCSS(:,3),pCSS(:,4),cssEffectiveSigma, ...
        mappingBetaCSS,betaCSS,sseCSS,r2CSS,corrCSS,validCSSPRF,validForK,nAtBoundary, ...
        'VariableNames',{'subNum','ROI','voxelIndex','sourceVoxelIndex','eccBin','originalEcc', ...
        'originalX','originalY','originalSigma','originalVexpl', ...
        'linearX','linearY','linearSigma','linearMappingBeta','linearBeta', ...
        'linearSSE','linearR2','linearCorr','validLinearPRF', ...
        'cssX','cssY','cssSigma','cssN','cssEffectiveSigma','cssMappingBeta', ...
        'cssBeta','cssSSE','cssR2','cssCorr','validCSSPRF','validForK','cssNAtBoundary'});
    if opts.verbose
        fprintf('V%d subject %d/%d complete: %d/%d fits usable for k\n', ...
                roi,s,nSub,nnz(validForK),nSample);
    end
end
hasVariables = cellfun(@(T) istable(T) && width(T) > 0,voxelTables);
if any(hasVariables)
    voxelTable = vertcat(voxelTables{hasVariables});
else
    voxelTable = table();
end
summaryTable = summarizeResults(kLinearBySubject,kCSSBySubject,nVoxBySubject, ...
    medianNBySubject,medianLinearSigmaBySubject, ...
    medianCssEffectiveSigmaBySubject,roi,opts);
out = struct();
out.ROI = roi;
out.options = opts;
out.voxelTable = voxelTable;
out.summaryTable = summaryTable;
out.kLinearBySubject = kLinearBySubject;
out.kCSSBySubject = kCSSBySubject;
out.nVoxBySubject = nVoxBySubject;
out.medianNBySubject = medianNBySubject;
out.medianLinearSigmaBySubject = medianLinearSigmaBySubject;
out.medianCssEffectiveSigmaBySubject = medianCssEffectiveSigmaBySubject;
end

function opts = defaults(opts)
if ~isfield(opts,'TR'), opts.TR = 1.2; end
if ~isfield(opts,'eccEdges'), opts.eccEdges = 0:0.25:8; end
if ~isfield(opts,'nPerBin'), opts.nPerBin = 10; end
if ~isfield(opts,'minVoxPerBin'), opts.minVoxPerBin = 5; end
if ~isfield(opts,'rngSeed'), opts.rngSeed = 1; end
if ~isfield(opts,'centreLimit'), opts.centreLimit = 8.5; end
if ~isfield(opts,'sigmaBounds'), opts.sigmaBounds = [0.05 8]; end
if ~isfield(opts,'nBounds'), opts.nBounds = [0.05 1.35]; end
if ~isfield(opts,'cssStarts'), opts.cssStarts = [0.12 0.33 0.75 1.15]; end
if ~isfield(opts,'minMappingR2')
    if isfield(opts,'minFullR2')
        opts.minMappingR2 = opts.minFullR2; % compatibility with older scripts
    else
        opts.minMappingR2 = 0;
    end
end
if ~isfield(opts,'kBounds'), opts.kBounds = [0 1]; end
if ~isfield(opts,'kGridStep'), opts.kGridStep = 0.05; end
if ~isfield(opts,'boundaryTol'), opts.boundaryTol = 0.01; end
if ~isfield(opts,'useParallel'), opts.useParallel = true; end
if ~isfield(opts,'verbose'), opts.verbose = true; end
if ~isfield(opts,'optim')
    opts.optim = optimset('Display','off','MaxIter',300,'MaxFunEvals',1500, ...
                          'TolX',1e-3,'TolFun',1e-6);
end

opts.eccEdges = double(opts.eccEdges(:).');
opts.sigmaBounds = double(opts.sigmaBounds(:).');
opts.nBounds = double(opts.nBounds(:).');
opts.kBounds = double(opts.kBounds(:).');
opts.cssStarts = double(opts.cssStarts(:).');
if numel(opts.eccEdges) < 2 || any(~isfinite(opts.eccEdges)) || any(diff(opts.eccEdges) <= 0)
    error('fitSampledPRFCSS:badEccEdges','eccEdges must increase strictly.');
end
if numel(opts.sigmaBounds) ~= 2 || opts.sigmaBounds(1) <= 0 || opts.sigmaBounds(2) <= opts.sigmaBounds(1)
    error('fitSampledPRFCSS:badSigmaBounds','sigmaBounds must be [positive lower, larger upper].');
end
if numel(opts.nBounds) ~= 2 || any(~isfinite(opts.nBounds)) || ...
        opts.nBounds(1) <= 0 || opts.nBounds(2) <= opts.nBounds(1)
    error('fitSampledPRFCSS:badNBounds','nBounds must contain two increasing positive values.');
end
if ~isscalar(opts.minMappingR2) || ~isfinite(opts.minMappingR2)
    error('fitSampledPRFCSS:badMappingR2','minMappingR2 must be a finite scalar.');
end
if any(opts.cssStarts < opts.nBounds(1) | opts.cssStarts > opts.nBounds(2))
    error('fitSampledPRFCSS:badCSSStarts','Every CSS start must lie within nBounds.');
end
if numel(opts.kBounds) ~= 2 || opts.kBounds(1) < 0 || opts.kBounds(2) > 1 || opts.kBounds(2) <= opts.kBounds(1)
    error('fitSampledPRFCSS:badKBounds','kBounds must increase within [0,1].');
end
end

function beta = refitBeta(p,n,Y,Afull,hrf,x,y,TR,hemIdx)
G = gaussianMatrix(p,x,y);
nVox = size(p,1);
num = zeros(nVox,1);
den = zeros(nVox,1);
for r = 1:numel(Y)
    drive = max(Afull{r}*G,0);
    P = cssPredict(drive,n,hrf,TR,hemIdx,r);
    good = isfinite(P) & isfinite(Y{r});
    P(~good) = 0;
    Yr = Y{r}; Yr(~good) = 0;
    num = num+sum(P.*Yr,1).';
    den = den+sum(P.^2,1).';
end
beta = num./den;
beta(~isfinite(beta) | den <= eps | beta <= 0) = NaN;
end

function [pLin,bLin,sseLin,r2Lin,corrLin,pCss,bCss,sseCss,r2Css,corrCss] = ...
    fitVoxelModels(Y,p0,args,h)
centreBounds = [-args.centreLimit args.centreLimit];
lbLin = [centreBounds(1) centreBounds(1) args.sigmaBounds(1)];
ubLin = [centreBounds(2) centreBounds(2) args.sigmaBounds(2)];
p0 = min(max(double(p0),lbLin+1e-6),ubLin-1e-6);
objLin = @(p) fullObjective(p,1,Y,args,h);
[pLin,sseLin] = fminsearchcon(objLin,p0,lbLin,ubLin,[],[],[],args.optim);
pLin = pLin(:).';
[sseLin,bLin,r2Lin,corrLin] = fullObjective(pLin,1,Y,args,h);
lbCss = [lbLin args.nBounds(1)];
ubCss = [ubLin args.nBounds(2)];
[pCss,bCss,sseCss,r2Css,corrCss] = deal(nan);
best = Inf;
for n0 = args.cssStarts
    sigma0 = min(max(pLin(3)*sqrt(n0),args.sigmaBounds(1)+1e-6),args.sigmaBounds(2)-1e-6);
    start = [pLin(1:2) sigma0 n0];
    objCss = @(p) fullObjective(p(1:3),p(4),Y,args,h);
    [candidate,f] = fminsearchcon(objCss,start,lbCss,ubCss,[],[],[],args.optim);
    candidate = candidate(:).';
    if isfinite(f) && f < best
        best = f;
        pCss = candidate;
    end
end
if all(isfinite(pCss))
    [sseCss,bCss,r2Css,corrCss] = fullObjective(pCss(1:3),pCss(4),Y,args,h);
else
    pCss = nan(1,4);
end
end

function [sse,beta,r2,rModel] = fullObjective(p,n,Y,args,h)
if hypot(p(1),p(2)) > args.centreLimit
    excess = hypot(p(1),p(2))-args.centreLimit;
    sse = 1e12*(1+excess^2); beta = NaN; r2 = NaN; rModel = NaN;
    return
end
G = gaussian(p,args.x,args.y);
P = cell(numel(args.Afull),1);
for r = 1:numel(P)
    drive = max(args.Afull{r}*G,0);
    P{r} = cssPredict(drive,n,args.hrf,args.TR,h,r);
end
[beta,sse,r2,rModel] = fitAmplitude(P,Y);
if ~isfinite(sse), sse = 1e12; end
end

function [beta,sse,r2,rModel] = fitAmplitude(P,Y)
num = 0; den = 0; sst = 0;
for r = 1:numel(P)
    good = isfinite(P{r}) & isfinite(Y{r});
    pr = P{r}(good); yr = Y{r}(good);
    num = num+sum(pr.*yr); den = den+sum(pr.^2); sst = sst+sum(yr.^2);
end
if den <= eps || sst <= eps
    beta = NaN; sse = Inf; r2 = NaN; rModel = NaN;
    return
end
beta = max(num/den,0);
sse = 0; allP = []; allY = [];
for r = 1:numel(P)
    good = isfinite(P{r}) & isfinite(Y{r});
    pr = beta*P{r}(good); yr = Y{r}(good);
    sse = sse+sum((yr-pr).^2);
    allP = [allP;pr(:)]; allY = [allY;yr(:)]; %#ok<AGROW>
end
r2 = 1-sse/sst;
if numel(allP) > 1 && std(allP) > 0 && std(allY) > 0
    C = corrcoef(allP,allY); rModel = C(1,2);
else
    rModel = NaN;
end
end

function k = fitK(p,n,beta,Y,Afull,Ascot,hrf,x,y,opts,hemIdx)
G = gaussianMatrix(p,x,y);
driveS = cell(numel(Afull),1); driveD = cell(numel(Afull),1);
for r = 1:numel(Afull)
    driveS{r} = max(Ascot{r}*G,0);
    driveD{r} = max((Afull{r}-Ascot{r})*G,0);
end
objective = @(v) kObjective(v,driveS,driveD,n,beta,Y,hrf,opts.TR,hemIdx);
grid = opts.kBounds(1):opts.kGridStep:opts.kBounds(2);
if grid(end) < opts.kBounds(2), grid = [grid opts.kBounds(2)]; end
sse = arrayfun(objective,grid);
[~,ii] = min(sse);
k = grid(ii);
if ii > 1 && ii < numel(grid)
    o = optimset('Display','off','TolX',1e-3);
    [candidate,f] = fminbnd(objective,grid(ii-1),grid(ii+1),o);
    if f < sse(ii), k = candidate; end
end
end

function sse = kObjective(k,driveS,driveD,n,beta,Y,hrf,TR,hemIdx)
sse = 0;
for r = 1:numel(Y)
    P = cssPredict(driveS{r}+k*driveD{r},n,hrf,TR,hemIdx,r);
    R = Y{r}-P.*beta(:).';
    R = R(isfinite(R));
    sse = sse+sum(R.^2);
end
if ~isfinite(sse), sse = Inf; end
end

function T = summarizeResults(kLin,kCss,nVox,medianNBySubject, ...
    medianLinearSigmaBySubject,medianCssEffectiveSigmaBySubject,roi,opts)
nBin = size(kLin,2);
ecc = ((opts.eccEdges(1:end-1)+opts.eccEdges(2:end))/2).';
[meanLin,seLin,meanCss,seCss,delta,seDelta] = deal(nan(nBin,1));
nSubjects = zeros(nBin,1); nVoxTotal = zeros(nBin,1);
[meanN,seN,meanLinearSigma,seLinearSigma,meanCssEffectiveSigma,seCssEffectiveSigma, ...
 meanEffectiveSizeRatio,seEffectiveSizeRatio] = ...
    deal(nan(nBin,1));
for b = 1:nBin
    q = isfinite(kLin(:,b)) & isfinite(kCss(:,b));
    nSubjects(b) = nnz(q); nVoxTotal(b) = sum(nVox(q,b));
    if any(q)
        lin = kLin(q,b); css = kCss(q,b); dif = css-lin;
        meanLin(b) = mean(lin); meanCss(b) = mean(css); delta(b) = mean(dif);
        seLin(b) = standardError(lin);
        seCss(b) = standardError(css);
        seDelta(b) = standardError(dif);
    end
    [meanN(b),seN(b)] = meanAndSE(medianNBySubject(:,b));
    [meanLinearSigma(b),seLinearSigma(b)] = meanAndSE(medianLinearSigmaBySubject(:,b));
    [meanCssEffectiveSigma(b),seCssEffectiveSigma(b)] = ...
        meanAndSE(medianCssEffectiveSigmaBySubject(:,b));
    ratio = medianCssEffectiveSigmaBySubject(:,b)./medianLinearSigmaBySubject(:,b);
    [meanEffectiveSizeRatio(b),seEffectiveSizeRatio(b)] = meanAndSE(ratio);
end
T = table(repmat(roi,nBin,1),ecc,nSubjects,nVoxTotal, ...
    meanLin,seLin,meanCss,seCss,delta,seDelta,meanN,seN, ...
    meanLinearSigma,seLinearSigma,meanCssEffectiveSigma,seCssEffectiveSigma, ...
    meanEffectiveSizeRatio,seEffectiveSizeRatio, ...
    'VariableNames',{'ROI','ecc','nSubjects','nVox','kLinear','kLinear_se', ...
    'kCSS','kCSS_se','deltaK','deltaK_se','medianCssN','medianCssN_se', ...
    'medianLinearSigma','medianLinearSigma_se','medianCssEffectiveSigma', ...
    'medianCssEffectiveSigma_se','medianEffectiveSizeRatio', ...
    'medianEffectiveSizeRatio_se'});
end

function [Afull,Ascot,stimTime] = loadRawStimuli(nRun,TR,x,y)
Afull = cell(nRun,1); Ascot = cell(nRun,1); stimTime = cell(nRun,1);
for r = 1:nRun
    [imgF,funcF] = loadScotomaStimuli(r,'logbar',TR);
    [imgS,funcS] = loadScotomaStimuli(r,'scotoma',TR);
    same = isequal(size(imgF),size(imgS)) && isequal(size(funcF.x),size(x)) && ...
           isequal(size(funcF.y),size(y)) && max(abs(funcF.x(:)-x(:))) <= 1e-10 && ...
           max(abs(funcF.y(:)-y(:))) <= 1e-10 && ...
           isequal(size(funcF.x),size(funcS.x)) && isequal(size(funcF.y),size(funcS.y)) && ...
           max(abs(funcF.x(:)-funcS.x(:))) <= 1e-10 && ...
           max(abs(funcF.y(:)-funcS.y(:))) <= 1e-10 && ...
           isequal(size(funcF.t),size(funcS.t)) && ...
           max(abs(funcF.t(:)-funcS.t(:))) <= 1e-10;
    if ~same
        error('fitSampledPRFCSS:stimulusMismatch','Stimulus or grid mismatch in run %d.',r);
    end
    nt = numel(funcF.t); nPix = numel(funcF.x);
    Afull{r} = reshape(double(imgF),[nPix nt]).';
    Ascot{r} = reshape(double(imgS),[nPix nt]).';
    stimTime{r} = double(funcF.t(:).');
end
end

function G = gaussian(p,x,y)
G = exp(-((double(x)-p(1)).^2+(double(y)-p(2)).^2)/(2*p(3)^2));
G = G(:);
mass = sum(G);
if ~isfinite(mass) || mass <= 0, G(:) = NaN; else, G = G/mass; end
end

function G = gaussianMatrix(p,x,y)
G = zeros(numel(x),size(p,1));
for v = 1:size(p,1), G(:,v) = gaussian(p(v,:),x,y); end
end

function Y = cssPredict(drive,n,hrf,TR,hemIdx,r)
% drive is [nT x nVox]. hrf{h}{r} is the HRF for hemisphere h on run r, and
% hemIdx says which hemisphere each column belongs to (scalar when every
% column shares one). The exponent is applied BEFORE convolution, so this
% cannot be folded into a pre-convolved design.
neural = bsxfun(@power,max(drive,0),double(n(:).'));
Y = centre(convByHemisphere(neural,hrf,TR,hemIdx,r));
end

function Y = selectVoxel(Yall,v)
Y = cell(size(Yall));
for r = 1:numel(Yall), Y{r} = Yall{r}(:,v); end
end

function Y = selectVoxels(Yall,q)
Y = cell(size(Yall));
for r = 1:numel(Yall), Y{r} = Yall{r}(:,q); end
end

function Z = centre(Z)
Z = Z-mean(Z,1,'omitnan');
end

function [m,se] = meanAndSE(x)
x = x(isfinite(x));
if isempty(x), m = NaN; se = NaN; return, end
m = mean(x);
se = standardError(x);
end

function se = standardError(x)
x = x(isfinite(x));
if numel(x) < 2, se = NaN; else, se = std(x,0)/sqrt(numel(x)); end
end
