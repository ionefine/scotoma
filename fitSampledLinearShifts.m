function out = fitSampledLinearShifts(subjectData,stimData,x,y,voxelTable,opts)
% fitSampledLinearShifts  Direct nonlinear pRF shift and filling-in fits.
%
% Independent x/y locations and full-field-fitted sigma/beta values come from
% fitSampledPRFCSS. Within each subject and eccentricity bin, this function
% directly minimizes scotoma-run SSE for three models using the same voxels:
%
%   fixed pRF:       k
%   shifted pRF:     delta-r, delta-theta, delta-sigma; k fixed at zero
%   shifted pRF + k: delta-r, delta-theta, delta-sigma, k
%
% This is a full nonlinear refit. There is no derivative approximation or
% iterative displacement step. Positive delta-r is radially outward; negative
% delta-r is toward the central scotoma. Positive delta-sigma enlarges the pRF.
% Beta remains fixed at its value estimated from the full-field runs.

if nargin < 6 || isempty(opts), opts = struct(); end
opts = defaults(opts);
if exist('fminsearchcon','file') ~= 2
    error('fitSampledLinearShifts:missingOptimizer','fminsearchcon must be on the MATLAB path.');
end
if ~iscell(subjectData), subjectData = {subjectData}; end
if ~iscell(stimData) || numel(stimData) ~= numel(subjectData)
    error('fitSampledLinearShifts:subjectCount','stimData must contain one entry per subject.');
end
required = {'subNum','voxelIndex','sourceVoxelIndex','originalEcc', ...
            'independentX','independentY','linearSigma','linearBeta','validLinearPRF'};
if ~istable(voxelTable) || ~all(ismember(required,voxelTable.Properties.VariableNames))
    error('fitSampledLinearShifts:badVoxelTable', ...
          'voxelTable must be the table produced by the revised fitSampledPRFCSS.');
end
validateattributes(x,{'numeric'},{'real','finite','nonempty'},mfilename,'x');
validateattributes(y,{'numeric'},{'real','finite','size',size(x)},mfilename,'y');
nSub = numel(subjectData); nBin = numel(opts.eccEdges)-1;
names = {'drNoFill','dthetaNoFill','dsigmaNoFill','drWithK','dthetaWithK', ...
         'dsigmaWithK','kFixed','kJoint','massInMedian','sseFixed','sseNoFill', ...
         'sseWithK','exitNoFill','exitWithK'};
R = struct();
for i = 1:numel(names), R.(names{i}) = nan(nSub,nBin); end
R.nVox = zeros(nSub,nBin); R.fixedAtBoundary = false(nSub,nBin);
R.noFillAtBoundary = false(nSub,nBin); R.withKAtBoundary = false(nSub,nBin);
R.subjectID = nan(nSub,1);
% Subjects are independent, so each one is fitted into its own row struct and
% the results are stitched together afterwards. parfor cannot slice fields of
% a struct, hence the cell of per-subject results.
rowOut = cell(nSub,1);
useParfor = opts.useParallel && license('test','Distrib_Computing_Toolbox');
if useParfor
    parfor s = 1:nSub
        rowOut{s} = fitOneSubject(subjectData{s},stimData{s},x,y,voxelTable,opts,nBin,s);
    end
else
    for s = 1:nSub
        rowOut{s} = fitOneSubject(subjectData{s},stimData{s},x,y,voxelTable,opts,nBin,s);
    end
end
for s = 1:nSub
    if isempty(rowOut{s}), continue, end
    R.subjectID(s) = rowOut{s}.subjectID;
    for i = 1:numel(names), R.(names{i})(s,:) = rowOut{s}.(names{i}); end
    R.nVox(s,:) = rowOut{s}.nVox;
    R.fixedAtBoundary(s,:) = rowOut{s}.fixedAtBoundary;
    R.noFillAtBoundary(s,:) = rowOut{s}.noFillAtBoundary;
    R.withKAtBoundary(s,:) = rowOut{s}.withKAtBoundary;
end
T = summarize(R,opts);
out = struct('options',opts,'summaryTable',T,'subjectResults',R);
if ismember('ROI',voxelTable.Properties.VariableNames)
    values = unique(voxelTable.ROI);
    if isscalar(values), out.ROI = values; else, out.ROI = NaN; end
else
    out.ROI = NaN;
end
end

function row = fitOneSubject(S,D,x,y,voxelTable,opts,nBin,s)
names = {'drNoFill','dthetaNoFill','dsigmaNoFill','drWithK','dthetaWithK', ...
         'dsigmaWithK','kFixed','kJoint','massInMedian','sseFixed','sseNoFill', ...
         'sseWithK','exitNoFill','exitWithK'};
row = struct();
for i = 1:numel(names), row.(names{i}) = nan(1,nBin); end
row.nVox = zeros(1,nBin); row.fixedAtBoundary = false(1,nBin);
row.noFillAtBoundary = false(1,nBin); row.withKAtBoundary = false(1,nBin);
    validateSubject(S,D,x,y,s);
    subNum = s; if isfield(S,'subNum'), subNum = S.subNum; end
    row.subjectID = subNum;
    rows = find(voxelTable.subNum == subNum & voxelTable.validLinearPRF & ...
        voxelTable.originalEcc >= opts.minShiftEcc & ...
        isfinite(voxelTable.linearSigma) & voxelTable.linearSigma > 0 & ...
        isfinite(voxelTable.linearBeta) & voxelTable.linearBeta > 0);
    if isempty(rows), return, end
    idx = double(voxelTable.voxelIndex(rows));
    if any(idx ~= round(idx) | idx < 1 | idx > size(S.prfXY,1)) || ...
            numel(unique(idx)) ~= numel(idx)
        error('fitSampledLinearShifts:badVoxelIndex','Invalid voxel indices for subject %g.',subNum);
    end
    if ~isfield(S,'sourceVoxelIndex') || any(double(S.sourceVoxelIndex(idx)) ~= ...
                                              double(voxelTable.sourceVoxelIndex(rows)))
        error('fitSampledLinearShifts:sourceVoxelMismatch', ...
            ['Voxel identity changed between the full-field and shift analyses for subject %g. ' ...
             'Use identical compiler settings.'],subNum);
    end
    xyTable = [voxelTable.independentX(rows),voxelTable.independentY(rows)];
    % Compare element-wise: xyTable is nVox-by-2, so flatten BOTH sides.
    % Flattening only the left side makes this an nVox*2 vs nVox-by-2
    % comparison, which errors for every subject with more than one voxel.
    if max(abs(xyTable(:)-reshape(double(S.prfXY(idx,:)),[],1))) > 1e-10
        error('fitSampledLinearShifts:locationMismatch', ...
              'Independent pRF locations no longer match subjectData for subject %g.',subNum);
    end
    P = struct('Yscot',{selectColumns(S.Yscot,idx)},'xy',xyTable, ...
        'sigma',double(voxelTable.linearSigma(rows)),'beta',double(voxelTable.linearBeta(rows)), ...
        'ecc',double(voxelTable.originalEcc(rows)),'hemIdx',double(S.hemIdx(idx)), ...
        'Sfull',{D.SfullRaw},'Sscot',{D.SscotRaw},'TR',D.TR, ...
        'hrf',{makeHRFs(D,numel(S.Yscot))},'x',double(x),'y',double(y));
    P.G0 = gaussianMatrix(P.xy,P.sigma,P.x,P.y);
    P.massIn = scotomaMass(P.G0,P.Sfull,P.Sscot);
    % Pixels the scotoma stimulus never drives contribute nothing to
    % Sscot*Gshift, so drop them from that product. The Gaussians are still
    % built and normalised on the FULL grid, so this is exact, not an
    % approximation - only the zero columns of Sscot are removed.
    P.pixKeep = false(1,size(P.G0,1));
    for r = 1:numel(P.Sscot)
        P.pixKeep = P.pixKeep | any(double(P.Sscot{r}) ~= 0,1);
    end
    P.SscotKeep = cell(1,numel(P.Sscot));
    for r = 1:numel(P.Sscot)
        P.SscotKeep{r} = double(P.Sscot{r}(:,P.pixKeep));
    end
    for b = 1:nBin
        q = discretize(P.ecc,opts.eccEdges) == b;
        if nnz(q) < opts.minVoxPerSubject, continue, end
        Q = subsetProblem(P,q);
        row.nVox(b) = nnz(q);
        [kFixed,sseFixed,fixedBoundary] = fitFixedK(Q,opts);
        if ~isfinite(kFixed)
            row.exitNoFill(b) = -2; row.exitWithK(b) = -2; continue
        end
        [pNo,sseNo,exitNo,noBoundary] = fitShift(Q,false,kFixed,opts);
        [pK,sseK,exitK,kBoundary] = fitShift(Q,true,kFixed,opts);
        row.exitNoFill(b) = exitNo; row.exitWithK(b) = exitK;
        row.fixedAtBoundary(b) = fixedBoundary;
        row.noFillAtBoundary(b) = noBoundary;
        row.withKAtBoundary(b) = kBoundary;
        if exitNo <= 0 || exitK <= 0 || any(~isfinite(pNo)) || any(~isfinite(pK)), continue, end
        row.drNoFill(b) = pNo(1); row.dthetaNoFill(b) = pNo(2);
        row.dsigmaNoFill(b) = pNo(3); row.drWithK(b) = pK(1);
        row.dthetaWithK(b) = pK(2); row.dsigmaWithK(b) = pK(3);
        row.kFixed(b) = kFixed; row.kJoint(b) = pK(4);
        row.massInMedian(b) = median(Q.massIn,'omitnan');
        row.sseFixed(b) = sseFixed; row.sseNoFill(b) = sseNo; row.sseWithK(b) = sseK;
    end
    if opts.verbose, fprintf('Voxel-shift fits: subject %d complete\n',s); end
end

function opts = defaults(opts)
if ~isfield(opts,'eccEdges'), opts.eccEdges = 0:0.25:5; end
if ~isfield(opts,'minVoxPerSubject'), opts.minVoxPerSubject = 5; end
if ~isfield(opts,'shiftBounds'), opts.shiftBounds = [-2 2;-2 2;-2 3]; end
if ~isfield(opts,'minSigma'), opts.minSigma = 0.05; end
if ~isfield(opts,'minShiftEcc'), opts.minShiftEcc = 1e-6; end
if ~isfield(opts,'kBounds'), opts.kBounds = [0 1]; end
if ~isfield(opts,'boundaryTol'), opts.boundaryTol = 0.01; end
if ~isfield(opts,'useParallel'), opts.useParallel = true; end
if ~isfield(opts,'verbose'), opts.verbose = true; end
if ~isfield(opts,'optim')
    opts.optim = optimset('Display','off','MaxIter',500,'MaxFunEvals',4000, ...
        'TolX',1e-4,'TolFun',1e-7);
end
opts.eccEdges = double(opts.eccEdges(:).'); opts.kBounds = double(opts.kBounds(:).');
if numel(opts.eccEdges) < 2 || any(~isfinite(opts.eccEdges)) || any(diff(opts.eccEdges) <= 0)
    error('fitSampledLinearShifts:badEdges','eccEdges must increase strictly.');
end
if ~isequal(size(opts.shiftBounds),[3 2]) || any(~isfinite(opts.shiftBounds(:))) || ...
        any(opts.shiftBounds(:,2) <= opts.shiftBounds(:,1))
    error('fitSampledLinearShifts:badShiftBounds','shiftBounds must be a valid 3-by-2 matrix.');
end
if ~isscalar(opts.minShiftEcc) || ~isfinite(opts.minShiftEcc) || opts.minShiftEcc < 0
    error('fitSampledLinearShifts:badMinEcc','minShiftEcc must be nonnegative.');
end
end

function validateSubject(S,D,x,y,s)
requiredS = {'Yscot','prfXY','sigma','hemIdx','sourceVoxelIndex'};
requiredD = {'SfullRaw','SscotRaw','hrfParams','tStim','TR'};
if ~isstruct(S) || ~all(isfield(S,requiredS)) || ~isstruct(D) || ~all(isfield(D,requiredD))
    error('fitSampledLinearShifts:badSubject','Subject %d is incomplete.',s);
end
nRun = numel(S.Yscot); nVox = size(S.prfXY,1);
if nRun == 0 || numel(D.SfullRaw) ~= nRun || numel(D.SscotRaw) ~= nRun || ...
        numel(D.tStim) ~= nRun || numel(S.hemIdx) ~= nVox
    error('fitSampledLinearShifts:dimensionMismatch','Subject %d dimensions disagree.',s);
end
for r = 1:nRun
    if size(S.Yscot{r},2) ~= nVox || ~isequal(size(D.SfullRaw{r}),size(D.SscotRaw{r})) || ...
            size(D.SfullRaw{r},1) ~= size(S.Yscot{r},1) || size(D.SfullRaw{r},2) ~= numel(x)
        error('fitSampledLinearShifts:runMismatch','Subject %d, run %d dimensions disagree.',s,r);
    end
end
if ~isequal(size(x),size(y)), error('fitSampledLinearShifts:gridMismatch','x and y differ.'); end
end

function [k,sse,atBoundary] = fitFixedK(P,opts)
objective = @(v) modelSSE([0 0 0 v],P,true,opts.minSigma);
grid = linspace(opts.kBounds(1),opts.kBounds(2),21);
values = arrayfun(objective,grid); [sse,ii] = min(values); k = grid(ii);
if ii > 1 && ii < numel(grid)
    [candidate,f] = fminbnd(objective,grid(ii-1),grid(ii+1),optimset('Display','off','TolX',1e-4));
    if f < sse, k = candidate; sse = f; end
end
if ~isfinite(sse), k = NaN; end
atBoundary = isfinite(k) && min(abs(k-opts.kBounds)) <= opts.boundaryTol;
end

function [p,sse,exitflag,atBoundary] = fitShift(P,includeK,kStart,opts)
lowerSigmaShift = max(opts.shiftBounds(3,1),opts.minSigma-min(P.sigma));
lb = [opts.shiftBounds(1,1),opts.shiftBounds(2,1),lowerSigmaShift];
ub = opts.shiftBounds(:,2).';
if lowerSigmaShift >= ub(3)
    p = nan(1,3+includeK); sse = Inf; exitflag = -3; atBoundary = false; return
end
starts = [0 0 0;0.1 0 0;-0.1 0 0;0 0 0.1];
if includeK
    lb = [lb opts.kBounds(1)]; ub = [ub opts.kBounds(2)];
    starts = [starts repmat(kStart,size(starts,1),1);0 0 0 0.5];
end
best = Inf; p = nan(1,numel(lb));
exitflag = -4;   % no start produced a finite SSE; must stay <= 0 so the
                 % failure counters in summarize() see it (NaN <= 0 is false)
for i = 1:size(starts,1)
    start = min(max(starts(i,:),lb+1e-6),ub-1e-6);
    [candidate,f,flag] = fminsearchcon(@(v) modelSSE(v,P,includeK,opts.minSigma), ...
        start,lb,ub,[],[],[],opts.optim);
    if isfinite(f) && f < best, best = f; p = candidate(:).'; exitflag = flag; end
end
sse = best;
atBoundary = any(abs(p-lb) <= opts.boundaryTol | abs(p-ub) <= opts.boundaryTol);
end

function sse = modelSSE(p,P,includeK,minSigma)
sigma = P.sigma+p(3);
if any(~isfinite(sigma) | sigma < minSigma), sse = Inf; return, end
theta = atan2(P.xy(:,2),P.xy(:,1));
xy = P.xy+[p(1)*cos(theta)-p(2)*sin(theta), ...
           p(1)*sin(theta)+p(2)*cos(theta)];
% Built on the full grid and normalised there, returned only for the pixels
% the scotoma stimulus can drive. Exactly equivalent, one fewer large copy.
Gshift = gaussianMatrix(xy,sigma,P.x,P.y,P.pixKeep);
if any(~isfinite(Gshift(:))), sse = Inf; return, end
if includeK, k = p(4); else, k = 0; end
sse = 0; nGood = 0;
for r = 1:numel(P.Yscot)
    % P.Kterm{r} is (Sfull-Sscot)*G0, which uses the UNSHIFTED pRF and so is
    % constant across the search. It is precomputed once per bin rather than
    % rebuilt on every objective evaluation.
    drive = P.SscotKeep{r}*Gshift+k*P.Kterm{r};
    pred = centre(convByHemisphere(drive,P.hrf,P.TR,P.hemIdx,r)).*P.beta(:).';
    residual = P.Yscot{r}-pred; good = isfinite(residual);
    sse = sse+sum(residual(good).^2); nGood = nGood+nnz(good);
end
if nGood == 0 || ~isfinite(sse), sse = Inf; end
end

function Q = subsetProblem(P,q)
Q = P; Q.xy = P.xy(q,:); Q.sigma = P.sigma(q); Q.beta = P.beta(q);
Q.ecc = P.ecc(q); Q.hemIdx = P.hemIdx(q); Q.G0 = P.G0(:,q); Q.massIn = P.massIn(q);
for r = 1:numel(P.Yscot), Q.Yscot{r} = P.Yscot{r}(:,q); end
% The missing-stimulus term uses the unshifted pRF, so it is constant for
% every objective evaluation in this bin. Only ~5% of pixels differ between
% the full and scotoma apertures, so this product is cheap to form once.
Q.Kterm = cell(1,numel(P.Yscot));
for r = 1:numel(P.Yscot)
    dS = double(P.Sfull{r})-double(P.Sscot{r});
    dPix = any(dS ~= 0,1);
    Q.Kterm{r} = dS(:,dPix)*Q.G0(dPix,:);
end
end

function mass = scotomaMass(G,Sfull,Sscot)
changed = false(1,size(G,1));
for r = 1:numel(Sfull)
    changed = changed | any(abs(double(Sfull{r})-double(Sscot{r})) > 1e-10,1);
end
mass = sum(G(changed,:),1).';
end

function T = summarize(R,opts)
nBin = numel(opts.eccEdges)-1; ecc = ((opts.eccEdges(1:end-1)+opts.eccEdges(2:end))/2).';
fields = {'drNoFill','dthetaNoFill','dsigmaNoFill','drWithK','dthetaWithK', ...
          'dsigmaWithK','kFixed','kJoint','massInMedian'};
M = struct(); SE = struct();
for i = 1:numel(fields)
    [M.(fields{i}),SE.(fields{i})] = columnMeanSE(R.(fields{i}));
end
valid = isfinite(R.kFixed) & isfinite(R.kJoint) & isfinite(R.drNoFill) & isfinite(R.drWithK);
nSubjects = sum(valid,1).'; nVox = sum(R.nVox.*valid,1).';
nFixedAtBoundary = sum(R.fixedAtBoundary & valid,1).';
nNoFillAtBoundary = sum(R.noFillAtBoundary & valid,1).';
nWithKAtBoundary = sum(R.withKAtBoundary & valid,1).';
nNoFillFailed = sum(R.exitNoFill <= 0 & R.nVox > 0,1).';
nWithKFailed = sum(R.exitWithK <= 0 & R.nVox > 0,1).';
T = table(ecc,nSubjects,nVox,M.massInMedian,SE.massInMedian, ...
    M.drNoFill,SE.drNoFill,M.dthetaNoFill,SE.dthetaNoFill,M.dsigmaNoFill,SE.dsigmaNoFill, ...
    M.drWithK,SE.drWithK,M.dthetaWithK,SE.dthetaWithK,M.dsigmaWithK,SE.dsigmaWithK, ...
    M.kFixed,SE.kFixed,M.kJoint,SE.kJoint,nFixedAtBoundary,nNoFillAtBoundary, ...
    nWithKAtBoundary,nNoFillFailed,nWithKFailed, ...
    'VariableNames',{'ecc','nSubjects','nVox','massIn_median','massIn_median_se', ...
    'dr_noFill','dr_noFill_se','dtheta_noFill','dtheta_noFill_se', ...
    'dsigma_noFill','dsigma_noFill_se','dr_withK','dr_withK_se', ...
    'dtheta_withK','dtheta_withK_se','dsigma_withK','dsigma_withK_se', ...
    'k_fixedPRF','k_fixedPRF_se','k_joint','k_joint_se','nFixedAtBoundary', ...
    'nNoFillAtBoundary','nWithKAtBoundary','nNoFillFailed','nWithKFailed'});
end

function [m,se] = columnMeanSE(X)
m = mean(X,1,'omitnan').'; n = sum(isfinite(X),1).';
se = std(X,0,1,'omitnan').'./sqrt(n); se(n < 2) = NaN;
end

function hrf = makeHRFs(D,nRun)
hrf = cell(numel(D.hrfParams),1);
for h = 1:numel(hrf)
    hrf{h} = cell(nRun,1);
    for r = 1:nRun
        hh = double(hrf_twogamma(D.hrfParams(h),D.tStim{r}));
        if isempty(hh) || any(~isfinite(hh)), error('fitSampledLinearShifts:badHRF','Invalid HRF.'); end
        hrf{h}{r} = hh(:);
    end
end
end

function Y = selectColumns(C,idx)
Y = cell(size(C)); for r = 1:numel(C), Y{r} = centre(double(C{r}(:,idx))); end
end

function G = gaussianMatrix(xy,sigma,x,y,keep)
% Unit-mass Gaussians on the stimulus grid, one column per voxel.
%
% keep (optional) is a pixel mask. The Gaussian is always evaluated and
% NORMALISED over the full grid, and only then restricted to the kept rows,
% so masking changes nothing about the returned values - it just avoids
% building and copying the full [nPix x nVox] matrix when the caller only
% needs the pixels the stimulus actually drives.
if nargin < 5 || isempty(keep), keep = true(numel(x),1); end
keep = keep(:);
G = nan(nnz(keep),size(xy,1));
for v = 1:size(xy,1)
    z = -((double(x)-xy(v,1)).^2+(double(y)-xy(v,2)).^2)/(2*sigma(v)^2);
    z = z-max(z(:)); g = exp(z(:)); mass = sum(g);
    if isfinite(mass) && mass > 0, G(:,v) = g(keep)/mass; end
end
end

function Z = centre(Z)
Z = Z-mean(Z,1,'omitnan');
end
