function out = fitKSpatialSummation(subjectData,stimData,nGrid,eccEdges,opts)
% fitKSpatialSummation  Re-estimate fixed-pRF k with CSS spatial summation.
%
% The pRF centers, sizes, and Gprf matrices are held fixed. For exponent n,
% the neural drive is
%
%   full:  (Afull*G)^n
%   scot:  [Ascot*G + k*(Afull-Ascot)*G]^n
%
% where A is the raw, unconvolved stimulus. The exponent is therefore
% applied after spatial pooling and before HRF convolution. Beta is
% re-estimated from the full-field runs for every voxel and n, then held
% fixed while k is fitted to the scotoma runs. One k is fitted per subject
% and eccentricity bin by pooled least squares across its voxels. All
% exponents use the same valid voxels. The n=1 fit is the existing fixed-pRF
% model with the same kBounds constraint.
%
% INPUTS
%   subjectData  cell array from compileStimAndSubData
%   stimData     matching cell array; each entry needs .hrfParams (1x2,
%                {L,R}) and .SfullRaw/.SscotRaw. subjectData needs .hemIdx.
%   nGrid        positive CSS exponents; include 1 for the linear baseline
%   eccEdges     eccentricity bin edges, using bins (lower,upper]
%   opts.TR       repetition time                         (default 1.2)
%       .kBounds  bounds for k                            (default [0 1])
%       .tolK     fminbnd tolerance                       (default 0.001)
%       .nCoarse  points in global k scan                 (default 21)
%       .minVox   minimum voxels per subject/bin          (default 20)
%       .verbose                                            (default true)
%
% OUTPUT
%   out.k          [subject x bin x exponent]
%   out.meanK      [bin x exponent], mean across subjects
%   out.semK       [bin x exponent]
%   out.deltaK     subject estimates minus the n=1 estimate
%   out.atBoundary estimates within tolK of a k bound
%   out.summary    long-format table of means and differences

if nargin < 5 || isempty(opts), opts = struct(); end
if ~isstruct(opts) || ~isscalar(opts)
    error('fitKSpatialSummation:badOptions','opts must be a scalar struct.');
end
if ~isfield(opts,'TR'),      opts.TR = 1.2; end
if ~isfield(opts,'kBounds'), opts.kBounds = [0 1]; end
if ~isfield(opts,'tolK'),    opts.tolK = 0.001; end
if ~isfield(opts,'nCoarse'), opts.nCoarse = 21; end
if ~isfield(opts,'minVox'),  opts.minVox = 20; end
if ~isfield(opts,'verbose'), opts.verbose = true; end
if ~iscell(subjectData), subjectData = {subjectData}; end
if isempty(subjectData) || ~isstruct(subjectData{1}) || ...
   ~isscalar(subjectData{1}) || ~isfield(subjectData{1},'Yfull') || ...
   ~iscell(subjectData{1}.Yfull)
    error('fitKSpatialSummation:noSubjects','subjectData is empty or invalid.');
end
if ~iscell(stimData) || numel(stimData) ~= numel(subjectData)
    error('fitKSpatialSummation:subjectCount', ...
          'stimData must contain one entry per subject.');
end
nGrid = nGrid(:).';
eccEdges = eccEdges(:).';
if ~isnumeric(nGrid) || ~isreal(nGrid) || isempty(nGrid) || ...
   any(~isfinite(nGrid) | nGrid <= 0)
    error('fitKSpatialSummation:badExponent','nGrid must contain positive finite values.');
end
if ~isnumeric(eccEdges) || ~isreal(eccEdges) || numel(eccEdges) < 2 || ...
   any(~isfinite(eccEdges)) || any(diff(eccEdges) <= 0)
    error('fitKSpatialSummation:badEdges','eccEdges must increase strictly.');
end
if ~isnumeric(opts.kBounds) || ~isreal(opts.kBounds) || ...
   numel(opts.kBounds) ~= 2 || any(~isfinite(opts.kBounds)) || ...
   opts.kBounds(2) <= opts.kBounds(1) || opts.kBounds(1) < 0
    error('fitKSpatialSummation:badBounds','kBounds must be [low high] with 0 <= low < high.');
end
if ~isnumeric(opts.TR) || ~isreal(opts.TR) || ~isscalar(opts.TR) || ...
   ~isfinite(opts.TR) || opts.TR <= 0 || ...
   ~isnumeric(opts.tolK) || ~isreal(opts.tolK) || ~isscalar(opts.tolK) || ...
   ~isfinite(opts.tolK) || opts.tolK <= 0 || ...
   ~isnumeric(opts.nCoarse) || ~isreal(opts.nCoarse) || ...
   ~isscalar(opts.nCoarse) || opts.nCoarse < 3 || opts.nCoarse ~= round(opts.nCoarse) || ...
   ~isnumeric(opts.minVox) || ~isreal(opts.minVox) || ...
   ~isscalar(opts.minVox) || opts.minVox < 1 || opts.minVox ~= round(opts.minVox)
    error('fitKSpatialSummation:badOptions', ...
          'TR, tolK, nCoarse, or minVox is invalid.');
end
if ~(islogical(opts.verbose) || isnumeric(opts.verbose)) || ...
   ~isscalar(opts.verbose) || ~ismember(opts.verbose,[0 1])
    error('fitKSpatialSummation:badVerbose','verbose must be a scalar logical value.');
end

nSub = numel(subjectData);
nBin = numel(eccEdges)-1;
nN = numel(nGrid);
nRun = numel(subjectData{1}.Yfull);
if nRun == 0
    error('fitKSpatialSummation:noRuns','No runs were found.');
end

% Raw stimuli are common across subjects. Keep them unconvolved here.
Afull = cell(nRun,1); Ascot = cell(nRun,1); stimTime = cell(nRun,1);
dPix = [];
gridX = []; gridY = [];
for r = 1:nRun
    [imgF,funcF] = loadScotomaStimuli(r,'logbar',opts.TR);
    [imgS,funcS] = loadScotomaStimuli(r,'scotoma',opts.TR);
    sameGrid = isequal(size(funcF.x),size(funcS.x)) && ...
               max(abs(funcF.x(:)-funcS.x(:))) < 1e-10 && ...
               max(abs(funcF.y(:)-funcS.y(:))) < 1e-10;
    sameTime = isequal(size(funcF.t),size(funcS.t)) && ...
               max(abs(funcF.t(:)-funcS.t(:))) < 1e-10;
    if ~sameGrid || ~sameTime || ~isequal(size(imgF),size(imgS)) || ...
       any(imgF(:) < -1e-10) || any(imgS(:) < -1e-10) || ...
       any(imgS(:) > imgF(:)+1e-10)
        error('fitKSpatialSummation:stimulusMismatch', ...
              'Raw full and scotoma stimuli disagree in run %d.',r);
    end
    if r == 1
        gridX = funcF.x; gridY = funcF.y;
    elseif ~isequal(size(gridX),size(funcF.x)) || ...
           max(abs(gridX(:)-funcF.x(:))) > 1e-10 || ...
           max(abs(gridY(:)-funcF.y(:))) > 1e-10
        error('fitKSpatialSummation:gridMismatch', ...
              'The raw stimulus grid differs in run %d.',r);
    end
    nt = numel(funcF.t);
    nPix = numel(funcF.x);
    Afull{r} = reshape(double(imgF),[nPix,nt]).';
    Ascot{r} = reshape(double(imgS),[nPix,nt]).';
    stimTime{r} = double(funcF.t(:).');
    changed = any(abs(Afull{r}-Ascot{r}) > 1e-10,1);
    if isempty(dPix), dPix = changed; else, dPix = dPix | changed; end
end

out = struct();
out.nGrid = nGrid;
out.eccEdges = eccEdges;
out.ecc = (eccEdges(1:end-1)+eccEdges(2:end))/2;
out.kBounds = opts.kBounds;
out.subjectIDs = (1:nSub).';
if all(cellfun(@(S) isstruct(S) && isscalar(S) && isfield(S,'subNum'),subjectData))
    out.subjectIDs = cellfun(@(S) S.subNum,subjectData(:));
end
out.k = nan(nSub,nBin,nN);
out.sse = nan(nSub,nBin,nN);
out.nVox = zeros(nSub,nBin,nN);
out.atBoundary = false(nSub,nBin,nN);
out.massMedian = nan(nSub,nBin);
out.beta = cell(nSub,nN);
out.commonVoxels = cell(nSub,1);
optimOpts = optimset('Display','off','TolX',opts.tolK);

for s = 1:nSub
    S = subjectData{s};
    required = {'Yfull','Yscot','Gprf','prfXY'};
    if ~isstruct(S) || ~isscalar(S) || ~all(isfield(S,required)) || ...
       ~iscell(S.Yfull) || ~iscell(S.Yscot) || ...
       numel(S.Yfull) ~= nRun || numel(S.Yscot) ~= nRun
        error('fitKSpatialSummation:badSubject','Subject %d has incomplete data.',s);
    end
    if ~isstruct(stimData{s}) || ~isscalar(stimData{s}) || ...
       ~all(isfield(stimData{s},{'hrfParams','SfullRaw','SscotRaw'})) || ...
       ~isfield(S,'hemIdx') || ...
       numel(stimData{s}.SfullRaw) ~= nRun || numel(stimData{s}.SscotRaw) ~= nRun
        error('fitKSpatialSummation:missingHRF', ...
              'stimData{%d} lacks hrfParams or its run predictors.',s);
    end
    G = double(S.Gprf);
    if size(G,1) ~= size(Afull{1},2) || size(G,2) ~= size(S.prfXY,1) || ...
       size(S.prfXY,2) ~= 2 || any(~isfinite(G(:)))
        error('fitKSpatialSummation:badPRF','Subject %d has invalid Gprf or prfXY.',s);
    end
    gMass = sum(G,1);
    if any(~isfinite(gMass) | gMass <= 0)
        error('fitKSpatialSummation:badPRF','Subject %d has invalid pRF mass.',s);
    end
    G = G./gMass;
    nVox = size(G,2);
    if nVox == 0
        error('fitKSpatialSummation:noVoxels','Subject %d contains no voxels.',s);
    end
    ecc = hypot(S.prfXY(:,1),S.prfXY(:,2));
    massIn = sum(G(dPix,:),1).';
    bin = nan(nVox,1);
    for b = 1:nBin
        q = ecc > eccEdges(b) & ecc <= eccEdges(b+1);
        bin(q) = b;
    end

    Yfull = cell(nRun,1); Yscot = cell(nRun,1);
    driveF = cell(nRun,1); driveS = cell(nRun,1); driveD = cell(nRun,1);
    % hrf{h}{r}: one HRF per hemisphere, per run.
    nHem = numel(stimData{s}.hrfParams);
    hemIdx = double(S.hemIdx(:));
    if numel(hemIdx) ~= nVox || any(hemIdx < 1 | hemIdx > nHem)
        error('fitKSpatialSummation:badHemIdx', ...
              'hemIdx for subject %d is missing or indexes outside hrfParams.',s);
    end
    hrf = cell(nHem,1);
    for h = 1:nHem
        hrf{h} = cell(nRun,1);
        for r = 1:nRun
            hh = double(hrf_twogamma(stimData{s}.hrfParams(h),stimTime{r}));
            hrf{h}{r} = hh(:);
            if any(~isfinite(hrf{h}{r}))
                error('fitKSpatialSummation:badHRF', ...
                      'HRF is invalid for subject %d, hemisphere %d, run %d.',s,h,r);
            end
        end
    end
    for r = 1:nRun
        if size(S.Yfull{r},1) ~= size(Afull{r},1) || ...
           ~isequal(size(S.Yfull{r}),size(S.Yscot{r})) || ...
           size(S.Yfull{r},2) ~= nVox
            error('fitKSpatialSummation:dimensionMismatch', ...
                  'Dimensions disagree for subject %d, run %d.',s,r);
        end
        Yfull{r} = centre(double(S.Yfull{r}));
        Yscot{r} = centre(double(S.Yscot{r}));
        driveF{r} = max(Afull{r}*G,0);
        driveS{r} = max(Ascot{r}*G,0);
        driveD{r} = max(driveF{r}-driveS{r},0);
        % The locally loaded raw stimulus must match the one
        % compileStimAndSubData stored. This catches a wrong reshape or a
        % stimulus/grid mismatch between the two independent load paths,
        % which is what the old stored-versus-raw n=1 check caught before
        % the stored designs stopped being pre-convolved.
        storedF = double(stimData{s}.SfullRaw{r});
        storedS = double(stimData{s}.SscotRaw{r});
        scale = max([1,max(abs(storedF(:))),max(abs(storedS(:)))]);
        if ~isequal(size(Afull{r}),size(storedF)) || ...
           ~isequal(size(Ascot{r}),size(storedS)) || ...
           max(abs(Afull{r}(:)-storedF(:))) > 1e-8*scale || ...
           max(abs(Ascot{r}(:)-storedS(:))) > 1e-8*scale
            error('fitKSpatialSummation:rawMismatch', ...
                  'Locally loaded and stored raw stimuli differ for subject %d, run %d.',s,r);
        end
    end

    % Estimate beta for every exponent first, then use their intersection
    % so differences in k cannot be caused by differences in voxel selection.
    betaByN = nan(nVox,nN);
    for j = 1:nN
        n = nGrid(j);
        num = zeros(1,nVox); den = zeros(1,nVox); nPair = zeros(1,nVox);
        for r = 1:nRun
            F = cssPredict(driveF{r},n,hrf,opts.TR,hemIdx,r);
            good = isfinite(F) & isfinite(Yfull{r});
            F(~good) = 0;
            Y = Yfull{r}; Y(~good) = 0;
            num = num + sum(F.*Y,1);
            den = den + sum(F.^2,1);
            nPair = nPair + sum(good,1);
        end
        beta = (num./den).';
        beta(den(:) <= eps | nPair(:) < 20 | ~isfinite(beta) | beta <= 0) = NaN;
        betaByN(:,j) = beta;
        out.beta{s,j} = beta;
    end
    commonBeta = all(isfinite(betaByN),2);
    out.commonVoxels{s} = commonBeta;
    for b = 1:nBin
        q = bin == b & commonBeta;
        if any(q), out.massMedian(s,b) = median(massIn(q),'omitnan'); end
    end

    for j = 1:nN
        n = nGrid(j);
        beta = betaByN(:,j);
        for b = 1:nBin
            q = bin == b & commonBeta;
            out.nVox(s,b,j) = nnz(q);
            if nnz(q) < opts.minVox, continue, end
            if abs(n-1) <= 10*eps
                % Exact reduction to the existing fixed-pRF linear model,
                % subject to the common kBounds constraint.
                [kBest,sseBest] = linearFit(driveS,driveD,Yscot,hrf, ...
                                            beta,q,opts.TR,opts.kBounds,hemIdx);
            else
                objective = @(k) scotSSE(k,driveS,driveD,Yscot,hrf,beta,q,n,opts.TR,hemIdx);
                coarseK = linspace(opts.kBounds(1),opts.kBounds(2),opts.nCoarse);
                coarseSSE = arrayfun(objective,coarseK);
                [sseBest,ii] = min(coarseSSE);
                kBest = coarseK(ii);
                if ii > 1 && ii < numel(coarseK)
                    [kLocal,sseLocal] = fminbnd(objective,coarseK(ii-1), ...
                                                coarseK(ii+1),optimOpts);
                    if sseLocal < sseBest, kBest = kLocal; sseBest = sseLocal; end
                end
                if ~isfinite(sseBest), kBest = NaN; end
            end
            out.k(s,b,j) = kBest;
            out.sse(s,b,j) = sseBest;
            out.atBoundary(s,b,j) = min(abs(kBest-opts.kBounds)) <= opts.tolK;
        end
        if opts.verbose
            fprintf('subject %d/%d, n = %.4g: fitted %d/%d bins\n', ...
                    s,nSub,n,nnz(isfinite(out.k(s,:,j))),nBin);
        end
    end
end

out.meanK = reshape(mean(out.k,1,'omitnan'),nBin,nN);
out.semK = reshape(std(out.k,0,1,'omitnan')./sqrt(sum(isfinite(out.k),1)),nBin,nN);
out.nSubjects = reshape(sum(isfinite(out.k),1),nBin,nN);
out.nAtBoundary = reshape(sum(out.atBoundary,1),nBin,nN);
if opts.verbose && any(out.atBoundary(:))
    warning('fitKSpatialSummation:kBoundary', ...
            '%d subject/bin estimates reached a k bound; inspect nAtBoundary.', ...
            nnz(out.atBoundary));
end
out.deltaK = nan(size(out.k));
iLinear = find(abs(nGrid-1) <= 10*eps,1);
if ~isempty(iLinear)
    for j = 1:nN
        out.deltaK(:,:,j) = out.k(:,:,j)-out.k(:,:,iLinear);
    end
    out.meanDeltaK = reshape(mean(out.deltaK,1,'omitnan'),nBin,nN);
    out.semDeltaK = reshape(std(out.deltaK,0,1,'omitnan')./ ...
                         sqrt(sum(isfinite(out.deltaK),1)),nBin,nN);
else
    out.meanDeltaK = nan(nBin,nN);
    out.semDeltaK = nan(nBin,nN);
end
[E,N] = ndgrid(out.ecc,nGrid);
out.summary = table(E(:),N(:),out.meanK(:),out.semK(:), ...
                    out.meanDeltaK(:),out.semDeltaK(:),out.nSubjects(:), ...
                    out.nAtBoundary(:), ...
    'VariableNames',{'ecc','exponent','k','k_sem','delta_k','delta_k_sem', ...
                     'nSubjects','nAtBoundary'});
end

function Y = cssPredict(drive,n,hrf,TR,hemIdx,r)
% hrf{h}{r} is hemisphere h on run r; hemIdx says which hemisphere each
% column of drive belongs to. The exponent is applied before convolution.
neural = bsxfun(@power,max(drive,0),double(n));
Y = centre(convByHemisphere(neural,hrf,TR,hemIdx,r));
end


function sse = scotSSE(k,driveS,driveD,Yscot,hrf,beta,q,n,TR,hemIdx)
% Pooled SSE for one subject, eccentricity bin, exponent, and candidate k.
sse = 0;
b = beta(q).';
nGood = 0;
for r = 1:numel(driveS)
    P = cssPredict(driveS{r}(:,q)+k*driveD{r}(:,q),n,hrf,TR,hemIdx(q),r);
    R = Yscot{r}(:,q)-P.*b;
    good = isfinite(R);
    sse = sse+sum(R(good).^2);
    nGood = nGood+nnz(good);
end
if nGood < 20*nnz(q) || ~isfinite(sse), sse = Inf; end
end

function [k,sse] = linearFit(driveS,driveD,Yscot,hrf,beta,q,TR,kBounds,hemIdx)
% Closed-form n=1 fit; algebraically identical to the fixed-pRF k model.
A = 0; c = 0; zz = 0;
b = beta(q).';
for r = 1:numel(driveS)
    P = cssPredict(driveS{r}(:,q),1,hrf,TR,hemIdx(q),r);
    K = cssPredict(driveD{r}(:,q),1,hrf,TR,hemIdx(q),r);
    Z = Yscot{r}(:,q)-P.*b;
    W = K.*b;
    good = isfinite(Z) & isfinite(W);
    A = A+sum(W(good).^2);
    c = c+sum(W(good).*Z(good));
    zz = zz+sum(Z(good).^2);
end
if ~isfinite(A) || A <= eps
    k = NaN; sse = Inf; return
end
k = min(max(c/A,kBounds(1)),kBounds(2));
sse = max(zz-2*k*c+k^2*A,0);
end

function Z = centre(Z)
Z = Z-mean(Z,1,'omitnan');
end
