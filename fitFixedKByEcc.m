function T = fitFixedKByEcc(subjectData,stimData,eccEdges,nBoot,opts)
% fitFixedKByEcc  Fast fixed-pRF estimates of k and k2 by eccentricity.
%
% This computes only the fixed-pRF coefficients. It uses the Gprf matrices
% already made by compileStimAndSubData, so it avoids the Gaussian
% derivatives and repeated re-linearisation required by the geometry fit.
% Beta is estimated from the full-field condition and held fixed.
%
% INPUTS
%   subjectData, stimData  cell arrays from compileStimAndSubData
%   eccEdges               eccentricity bin edges
%   nBoot                  subject bootstrap replicates (default 2000)
%   opts.alpha             interval alpha             (default 0.05)
%       .minEcc            minimum pRF eccentricity   (default 0.1)
%       .minVox            minimum voxels per bin     (default 20)
%
% OUTPUT
%   T  one row per eccentricity bin, directly usable by
%      plotKResponseDecomposition

if nargin < 4 || isempty(nBoot), nBoot = 2000; end
if nargin < 5 || isempty(opts), opts = struct(); end
if ~isfield(opts,'alpha'), opts.alpha = 0.05; end
if ~isfield(opts,'minEcc'), opts.minEcc = 0.1; end
if ~isfield(opts,'minVox'), opts.minVox = 20; end
if ~iscell(subjectData), subjectData = {subjectData}; end
if ~iscell(stimData) || numel(stimData) ~= numel(subjectData)
    error('fitFixedKByEcc:subjectCount', ...
          'stimData must contain one entry per subject.');
end
eccEdges = double(eccEdges(:).');
if numel(eccEdges) < 2 || any(~isfinite(eccEdges)) || any(diff(eccEdges) <= 0)
    error('fitFixedKByEcc:badEdges','eccEdges must increase strictly.');
end
if ~isscalar(nBoot) || nBoot < 0 || nBoot ~= round(nBoot)
    error('fitFixedKByEcc:badBootstrap','nBoot must be a nonnegative integer.');
end
if ~isscalar(opts.alpha) || opts.alpha <= 0 || opts.alpha >= 1 || ...
   ~isscalar(opts.minEcc) || ~isfinite(opts.minEcc) || opts.minEcc < 0 || ...
   ~isscalar(opts.minVox) || opts.minVox < 1 || opts.minVox ~= round(opts.minVox)
    error('fitFixedKByEcc:badOptions','alpha, minEcc, or minVox is invalid.');
end

nSub = numel(subjectData);
nBin = numel(eccEdges)-1;
[Ak,ck,zzk,dfk,Ak2,ck2,zzk2,dfk2,nVox] = deal(zeros(nSub,nBin));
massByBin = cell(nBin,1);

for s = 1:nSub
    S = subjectData{s}; D = stimData{s};
    requiredS = {'Yfull','Yscot','Gprf','prfXY','sigma','hemIdx'};
    if ~isstruct(S) || ~all(isfield(S,requiredS)) || ...
       ~isstruct(D) || ~all(isfield(D,{'SfullRaw','SscotRaw','hrfParams','tStim','TR'}))
        error('fitFixedKByEcc:badSubject','Subject %d is incomplete.',s);
    end
    nRun = numel(S.Yfull);
    if nRun == 0 || numel(S.Yscot) ~= nRun || ...
       numel(D.SfullRaw) ~= nRun || numel(D.SscotRaw) ~= nRun
        error('fitFixedKByEcc:runCount','Run counts disagree for subject %d.',s);
    end
    % One HRF per hemisphere. The designs arrive unconvolved, so each
    % voxel's projected time course is convolved with its own HRF.
    nHem = numel(D.hrfParams);
    hemIdx = double(S.hemIdx(:).');
    if any(hemIdx < 1 | hemIdx > nHem)
        error('fitFixedKByEcc:badHemIdx','hemIdx indexes outside hrfParams for subject %d.',s);
    end
    hrf = cell(nHem,1);
    for h = 1:nHem
        hrf{h} = cell(nRun,1);
        for r = 1:nRun
            hh = double(hrf_twogamma(D.hrfParams(h),D.tStim{r}));
            hrf{h}{r} = hh(:);
        end
    end
    G = double(S.Gprf);
    gMass = sum(G,1);
    if any(~isfinite(G(:))) || any(~isfinite(gMass) | gMass <= 0)
        error('fitFixedKByEcc:badPRF','Invalid Gprf for subject %d.',s);
    end
    G = G./gMass;
    n = size(G,2);
    ecc = hypot(S.prfXY(:,1),S.prfXY(:,2));
    if numel(ecc) ~= n || numel(S.sigma) ~= n
        error('fitFixedKByEcc:voxelCount','Voxel dimensions disagree for subject %d.',s);
    end

    F = cell(nRun,1); P = cell(nRun,1); good0 = cell(nRun,1);
    dPix = false(1,size(G,1));
    num = zeros(1,n); den = zeros(1,n); nGood = zeros(1,n);
    for r = 1:nRun
        if ~isequal(size(D.SfullRaw{r}),size(D.SscotRaw{r})) || ...
           size(D.SfullRaw{r},2) ~= size(G,1) || ...
           ~isequal(size(S.Yfull{r}),size(S.Yscot{r})) || ...
           size(S.Yfull{r},1) ~= size(D.SfullRaw{r},1) || ...
           size(S.Yfull{r},2) ~= n
            error('fitFixedKByEcc:dimensionMismatch', ...
                  'Dimensions disagree for subject %d, run %d.',s,r);
        end
        Sf = double(D.SfullRaw{r}); Ss = double(D.SscotRaw{r});
        dPix = dPix | any(abs(Sf-Ss) > 1e-8,1);
        % Project onto each pRF, then convolve with that voxel's HRF.
        F{r} = centre(convByHemisphere(Sf*G,hrf,D.TR,hemIdx,r));
        P{r} = centre(convByHemisphere(Ss*G,hrf,D.TR,hemIdx,r));
        Yf = centre(double(S.Yfull{r}));
        Ys = centre(double(S.Yscot{r}));
        good0{r} = isfinite(F{r}) & isfinite(P{r}) & ...
                   isfinite(Yf) & isfinite(Ys);
        Fr = F{r}; Yr = Yf;
        Fr(~good0{r}) = 0; Yr(~good0{r}) = 0;
        num = num+sum(Fr.*Yr,1);
        den = den+sum(Fr.^2,1);
        nGood = nGood+sum(good0{r},1);
    end
    beta = num./den;
    beta(~isfinite(beta) | beta <= 0 | den <= eps | nGood < 20) = NaN;

    a = zeros(1,n); c = zeros(1,n); zz = zeros(1,n); ng = zeros(1,n);
    a2 = zeros(1,n); c2 = zeros(1,n); zz2 = zeros(1,n); ng2 = zeros(1,n);
    for r = 1:nRun
        Ys = centre(double(S.Yscot{r}));
        Z = Ys-P{r}.*beta;
        W = (F{r}-P{r}).*beta;
        W2 = F{r}.*beta;
        good = good0{r} & isfinite(Z) & isfinite(W);
        good2 = good0{r} & isfinite(Z) & isfinite(W2);
        Zk = Z; Wk = W; Zk(~good) = 0; Wk(~good) = 0;
        Z2 = Z; Wf = W2; Z2(~good2) = 0; Wf(~good2) = 0;
        a = a+sum(Wk.^2,1); c = c+sum(Wk.*Zk,1); zz = zz+sum(Zk.^2,1);
        a2 = a2+sum(Wf.^2,1); c2 = c2+sum(Wf.*Z2,1); zz2 = zz2+sum(Z2.^2,1);
        ng = ng+sum(good,1); ng2 = ng2+sum(good2,1);
    end
    ok = isfinite(beta) & a > eps & a2 > eps & ng >= 20 & ng2 >= 20 & ...
         isfinite(ecc(:).') & ecc(:).' >= opts.minEcc;
    mass = sum(G(dPix,:),1);
    bin = discretize(ecc,eccEdges);
    for b = 1:nBin
        q = ok(:) & bin(:) == b;
        if ~any(q), continue, end
        Ak(s,b) = sum(a(q)); ck(s,b) = sum(c(q)); zzk(s,b) = sum(zz(q));
        Ak2(s,b) = sum(a2(q)); ck2(s,b) = sum(c2(q)); zzk2(s,b) = sum(zz2(q));
        dfk(s,b) = sum(max(ng(q)-nRun,1));
        dfk2(s,b) = sum(max(ng2(q)-nRun,1));
        nVox(s,b) = nnz(q);
        massByBin{b} = [massByBin{b};mass(q).']; %#ok<AGROW>
    end
    fprintf('fixed-pRF fit: subject %d/%d\n',s,nSub);
end

ecc = ((eccEdges(1:end-1)+eccEdges(2:end))/2).';
[k,kLo,kHi,nTotal] = poolScalar(Ak,ck,zzk,dfk,nVox,nBoot,opts.alpha,opts.minVox);
[k2,k2Lo,k2Hi] = poolScalar(Ak2,ck2,zzk2,dfk2,nVox,nBoot,opts.alpha,opts.minVox);
massMedian = nan(nBin,1); massRMS = nan(nBin,1);
for b = 1:nBin
    if ~isempty(massByBin{b})
        massMedian(b) = median(massByBin{b});
        massRMS(b) = sqrt(mean(massByBin{b}.^2));
    end
end
nSubjects = sum(nVox > 0,1).';
T = table(ecc,massMedian,massRMS,nTotal,nSubjects,k,kLo,kHi,k2,k2Lo,k2Hi, ...
    'VariableNames',{'ecc','massIn_median','massIn_rms','nVox_fixedPRF', ...
                     'nSubjects','k','k_lo','k_hi','k2','k2_lo','k2_hi'});
end

function [b,lo,hi,nTotal] = poolScalar(A,c,zz,dof,nVox,nBoot,alpha,minVox)
nBin = size(A,2); nSub = size(A,1);
[b,lo,hi] = deal(nan(nBin,1));
nTotal = sum(nVox,1).';
for j = 1:nBin
    aa = sum(A(:,j)); cc = sum(c(:,j));
    if aa <= eps || nTotal(j) < minVox, continue, end
    b(j) = cc/aa;
    if nBoot > 0 && nSub >= 3
        idx = randi(nSub,nSub,nBoot);
        ab = sum(reshape(A(idx,j),nSub,nBoot),1);
        cb = sum(reshape(c(idx,j),nSub,nBoot),1);
        bs = cb./ab;
        lo(j) = percentile(bs,alpha/2);
        hi(j) = percentile(bs,1-alpha/2);
    else
        sse = max(sum(zz(:,j))-b(j)*cc,0);
        se = sqrt((sse/max(sum(dof(:,j))-1,1))/aa);
        z = sqrt(2)*erfinv(1-alpha);
        lo(j) = b(j)-z*se; hi(j) = b(j)+z*se;
    end
end
end

function v = percentile(z,p)
z = sort(z(isfinite(z)));
if isempty(z), v = NaN; return, end
i = max(1,min(numel(z),ceil(p*numel(z))));
v = z(i);
end

function Z = centre(Z)
Z = Z-mean(Z,1,'omitnan');
end
