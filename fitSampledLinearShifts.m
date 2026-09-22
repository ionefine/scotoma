function out = fitSampledLinearShifts(subjectData,stimData,x,y,sampleTable,opts)
% fitSampledLinearShifts  Linear-pRF shifts with and without filling-in.
%
% This is deliberately separate from the CSS analysis. It reuses the
% voxel indices sampled by fitSampledPRFCSS, but retains voxels according
% to the linear full-field fit only. The full-field linear pRF estimates
% replace the original pRF parameters before the scotoma data are fitted.
%
% The fixed-pRF, shift-only, and shift-plus-k models are summarized from
% exactly the same successful voxels within each subject and bin:
%   no filling:  beta*[Sscot*G(p + delta)]
%   with k:      beta*[Sscot*G(p + delta) + k*Sdiff*G(p)]
%
% delta contains radial shift, tangential shift, and sigma change. beta is
% re-estimated from the full-field data and then fixed in both models.
% Gauss-Newton re-linearisation is handled by prfShiftIterate.

if nargin < 6 || isempty(opts), opts = struct(); end
opts = defaults(opts);
if ~iscell(subjectData), subjectData = {subjectData}; end
if ~iscell(stimData) || numel(stimData) ~= numel(subjectData)
    error('fitSampledLinearShifts:subjectCount', ...
          'stimData must contain one entry per subject.');
end
required = {'subNum','voxelIndex','linearX','linearY','linearSigma', ...
            'linearBeta','linearR2'};
if ~istable(sampleTable) || ~all(ismember(required,sampleTable.Properties.VariableNames))
    error('fitSampledLinearShifts:badSampleTable', ...
          'sampleTable must be the voxelTable produced by fitSampledPRFCSS.');
end
validateattributes(x,{'numeric'},{'real','finite','nonempty'},mfilename,'x');
validateattributes(y,{'numeric'},{'real','finite','size',size(x)},mfilename,'y');

[sampledData,sampleInfo] = buildSampledData(subjectData,sampleTable,x,y,opts);
if sum(sampleInfo.nVox) == 0
    error('fitSampledLinearShifts:noVoxels','No sampled voxels passed the linear-fit criteria.');
end

iter = struct('maxIter',opts.maxIter,'tol',opts.tol,'damp',opts.damp, ...
              'maxStep',opts.maxStep,'fitOpts',opts.fitOpts, ...
              'fillMeasure','k','useParallel',opts.useParallel);
% Column 4 is the fixed-pRF k model. It also supplies the common zero-offset
% cross-products used to start the two geometry models.
[fitFixed,~,~] = prfShiftIterate(sampledData,stimData,x,y,opts.eccEdges,4,iter);
iter.initialFits = fitFixed;
[fitNo,deltaNo,historyNo] = prfShiftIterate(sampledData,stimData,x,y, ...
                                             opts.eccEdges,1:3,iter);
[fitWithK,deltaWithK,historyWithK] = prfShiftIterate(sampledData,stimData,x,y, ...
                                                     opts.eccEdges,1:4,iter);

[T,subjectResults] = summarizeFits(fitFixed,fitNo,fitWithK,deltaNo,deltaWithK,opts);
out = struct();
if ismember('ROI',sampleTable.Properties.VariableNames)
    roiValues = unique(sampleTable.ROI);
    if isscalar(roiValues), out.ROI = roiValues; else, out.ROI = NaN; end
else
    out.ROI = NaN;
end
out.options = opts;
out.sampleInfo = sampleInfo;
out.summaryTable = T;
out.subjectResults = subjectResults;
out.fitFixed = fitFixed;
out.fitNoFill = fitNo;
out.fitWithK = fitWithK;
out.deltaNoFill = deltaNo;
out.deltaWithK = deltaWithK;
out.historyNoFill = historyNo;
out.historyWithK = historyWithK;
end

function opts = defaults(opts)
if ~isfield(opts,'eccEdges'), opts.eccEdges = 0:0.25:8; end
if ~isfield(opts,'minLinearR2'), opts.minLinearR2 = 0; end
if ~isfield(opts,'minVoxPerSubject'), opts.minVoxPerSubject = 5; end
if ~isfield(opts,'maxIter'), opts.maxIter = 3; end
if ~isfield(opts,'tol'), opts.tol = 0.01; end
if ~isfield(opts,'damp'), opts.damp = 1; end
if ~isfield(opts,'maxStep'), opts.maxStep = 0.5; end
if ~isfield(opts,'useParallel'), opts.useParallel = true; end
if ~isfield(opts,'fitOpts')
    opts.fitOpts = struct('minEcc',1e-6,'hPos',0.1,'hSig',0.1,'verbose',false);
end
opts.eccEdges = double(opts.eccEdges(:).');
if numel(opts.eccEdges) < 2 || any(~isfinite(opts.eccEdges)) || ...
   any(diff(opts.eccEdges) <= 0)
    error('fitSampledLinearShifts:badEdges','eccEdges must increase strictly.');
end
if ~isscalar(opts.minLinearR2) || ~isfinite(opts.minLinearR2) || ...
   ~isscalar(opts.minVoxPerSubject) || opts.minVoxPerSubject < 1 || ...
   opts.minVoxPerSubject ~= round(opts.minVoxPerSubject)
    error('fitSampledLinearShifts:badOptions','Invalid fit threshold or voxel-count option.');
end
end

function [sampled,info] = buildSampledData(subjectData,T,x,y,opts)
nSub = numel(subjectData);
sampled = cell(nSub,1);
nVox = zeros(nSub,1);
subNum = nan(nSub,1);
for s = 1:nSub
    S = subjectData{s};
    needed = {'Yfull','Yscot','prfXY','sigma','w_vox'};
    if ~isstruct(S) || ~all(isfield(S,needed))
        error('fitSampledLinearShifts:badSubject','Subject %d is incomplete.',s);
    end
    if isfield(S,'subNum'), subNum(s) = S.subNum; else, subNum(s) = s; end
    q = T.subNum == subNum(s) & isfinite(T.voxelIndex) & ...
        isfinite(T.linearX) & isfinite(T.linearY) & ...
        isfinite(T.linearSigma) & T.linearSigma > 0 & ...
        isfinite(T.linearBeta) & T.linearBeta > 0 & ...
        isfinite(T.linearR2) & T.linearR2 >= opts.minLinearR2;
    if ismember('validLinearPRF',T.Properties.VariableNames)
        q = q & T.validLinearPRF;
    end
    rows = find(q);
    idx = double(T.voxelIndex(rows));
    nAll = size(S.prfXY,1);
    if any(idx ~= round(idx) | idx < 1 | idx > nAll) || numel(unique(idx)) ~= numel(idx)
        error('fitSampledLinearShifts:badVoxelIndex', ...
              'Sampled voxel indices are invalid for subject %g.',subNum(s));
    end
    sampled{s} = S;
    for r = 1:numel(S.Yfull)
        sampled{s}.Yfull{r} = S.Yfull{r}(:,idx);
        sampled{s}.Yscot{r} = S.Yscot{r}(:,idx);
    end
    sampled{s}.prfXY = [T.linearX(rows),T.linearY(rows)];
    sampled{s}.sigma = T.linearSigma(rows);
    sampled{s}.w_vox = S.w_vox(idx);
    % hemisphere must follow the voxels; sampled{s} = S above copied the
    % full-length labels, which would misalign against the subset fields.
    if isfield(S,'hemisphere'), sampled{s}.hemisphere = S.hemisphere(idx); end
    if isfield(S,'hemIdx'),     sampled{s}.hemIdx     = S.hemIdx(idx);     end
    sampled{s}.Gprf = gaussianMatrix(sampled{s}.prfXY,sampled{s}.sigma,x,y);
    nVox(s) = numel(idx);
end
info = table(subNum,nVox);
end

function G = gaussianMatrix(xy,sigma,x,y)
nVox = size(xy,1);
G = zeros(numel(x),nVox);
for v = 1:nVox
    pRF.center = xy(v,:);
    pRF.sig = sigma(v);
    pRF.ar = 1;
    g = Gauss(pRF,x,y,1);
    z = sum(g);
    if isfinite(z) && z > 0, G(:,v) = g(:)/z; end
end
end

function [T,subjectResults] = summarizeFits(fitFixed,fitNo,fitWithK,deltaNo,deltaWithK,opts)
nBin = numel(opts.eccEdges)-1;
nSub = numel(fitNo);
ecc = ((opts.eccEdges(1:end-1)+opts.eccEdges(2:end))/2).';
z = nan(nBin,1);
[nSubjects,nVox,massInMedian] = deal(z);
[drNo,drNoSE,dtNo,dtNoSE,dsNo,dsNoSE] = deal(z);
[drK,drKSE,dtK,dtKSE,dsK,dsKSE] = deal(z);
[kJoint,kJointSE,kFixed,kFixedSE] = deal(z);
subjectResults = struct('drNoFill',nan(nSub,nBin),'dthetaNoFill',nan(nSub,nBin), ...
    'dsigmaNoFill',nan(nSub,nBin),'drWithK',nan(nSub,nBin), ...
    'dthetaWithK',nan(nSub,nBin),'dsigmaWithK',nan(nSub,nBin), ...
    'kFixed',nan(nSub,nBin),'kJoint',nan(nSub,nBin), ...
    'massInMedian',nan(nSub,nBin),'nVox',zeros(nSub,nBin));
for b = 1:nBin
    subjectNo = nan(nSub,3);
    subjectK = nan(nSub,4);
    subjectFixed = nan(nSub,1);
    subjectMass = nan(nSub,1);
    subjectNVox = zeros(nSub,1);
    for s = 1:nSub
        inBin = inEccentricityBin(fitNo{s}.ecc,opts.eccEdges,b);
        common = inBin & fitFixed{s}.ok & fitNo{s}.ok & fitWithK{s}.ok;
        if nnz(common) < opts.minVoxPerSubject, continue, end
        [bFixed,nFixed] = prfShiftPool(fitFixed{s},4,common,'k');
        [bNo,nNo] = prfShiftPool(fitNo{s},1:3,common,'k');
        [bK,nK] = prfShiftPool(fitWithK{s},1:4,common,'k');
        if min([nFixed,nNo,nK]) < opts.minVoxPerSubject || ...
                ~isfinite(bFixed) || any(~isfinite(bNo)) || any(~isfinite(bK))
            continue
        end
        subjectFixed(s) = bFixed;
        subjectNo(s,:) = bNo(:).'+deltaNo(b,:);
        subjectK(s,:) = [bK(1:3).'+deltaWithK(b,:),bK(4)];
        subjectMass(s) = median(fitFixed{s}.massIn(common),'omitnan');
        subjectNVox(s) = nNo;
    end
    subjectResults.drNoFill(:,b) = subjectNo(:,1);
    subjectResults.dthetaNoFill(:,b) = subjectNo(:,2);
    subjectResults.dsigmaNoFill(:,b) = subjectNo(:,3);
    subjectResults.drWithK(:,b) = subjectK(:,1);
    subjectResults.dthetaWithK(:,b) = subjectK(:,2);
    subjectResults.dsigmaWithK(:,b) = subjectK(:,3);
    subjectResults.kJoint(:,b) = subjectK(:,4);
    subjectResults.kFixed(:,b) = subjectFixed;
    subjectResults.massInMedian(:,b) = subjectMass;
    subjectResults.nVox(:,b) = subjectNVox;
    included = all(isfinite(subjectNo),2) & all(isfinite(subjectK),2) & ...
               isfinite(subjectFixed);
    nSubjects(b) = nnz(included);
    nVox(b) = sum(subjectNVox(included));
    if ~any(included), continue, end
    massInMedian(b) = mean(subjectMass(included),'omitnan');
    [drNo(b),drNoSE(b)] = meanAndSE(subjectNo(included,1));
    [dtNo(b),dtNoSE(b)] = meanAndSE(subjectNo(included,2));
    [dsNo(b),dsNoSE(b)] = meanAndSE(subjectNo(included,3));
    [drK(b),drKSE(b)] = meanAndSE(subjectK(included,1));
    [dtK(b),dtKSE(b)] = meanAndSE(subjectK(included,2));
    [dsK(b),dsKSE(b)] = meanAndSE(subjectK(included,3));
    [kJoint(b),kJointSE(b)] = meanAndSE(subjectK(included,4));
    [kFixed(b),kFixedSE(b)] = meanAndSE(subjectFixed(included));
end
subjectResults.ecc = ecc;
T = table(ecc,nSubjects,nVox,massInMedian, ...
    drNo,drNoSE,dtNo,dtNoSE,dsNo,dsNoSE, ...
    drK,drKSE,dtK,dtKSE,dsK,dsKSE,kFixed,kFixedSE,kJoint,kJointSE, ...
    'VariableNames',{'ecc','nSubjects','nVox','massIn_median', ...
    'dr_noFill','dr_noFill_se','dtheta_noFill','dtheta_noFill_se', ...
    'dsigma_noFill','dsigma_noFill_se','dr_withK','dr_withK_se', ...
    'dtheta_withK','dtheta_withK_se','dsigma_withK','dsigma_withK_se', ...
    'k_fixedPRF','k_fixedPRF_se','k_joint','k_joint_se'});
end

function q = inEccentricityBin(ecc,edges,b)
if b == numel(edges)-1
    q = ecc >= edges(b) & ecc <= edges(b+1);
else
    q = ecc >= edges(b) & ecc < edges(b+1);
end
end

function [m,se] = meanAndSE(x)
x = x(isfinite(x));
if isempty(x), m = NaN; se = NaN; return, end
m = mean(x);
if numel(x) < 2, se = NaN; else, se = std(x,0)/sqrt(numel(x)); end
end
