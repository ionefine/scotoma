function [subjectData,stimData,x,y] = compileStimAndSubData(compileOpts)
% compileStimAndSubData  Load matched BOLD, stimulus, and original pRF data.
%
% HRF parameters are read per hemisphere by loadHRFParams. The hrfParams
% stored inside each pRF MAT file is the right-hemisphere fit for every
% subject and is deliberately ignored; compileOpts.hrfParams overrides the
% table if supplied (one struct for both hemispheres, or a two-element
% struct/cell array ordered {L,R}).
%
% VOXEL ELIGIBILITY (applied once here and reused by every later analysis):
%   1. voxel belongs to the requested visual area;
%   2. x, y, sigma, and variance explained are finite;
%   3. the original pRF centre is inside radRange (endpoints included when
%      edgeSigma is zero; otherwise the requested Gaussian extent must fit);
%   4. original pRF variance explained is strictly greater than minvexpl;
%   5. original sigma is greater than or equal to minSigma.
% No voxel is rejected here because of its scotoma response.

if nargin < 1 || isempty(compileOpts), compileOpts = struct(); end
if ~isfield(compileOpts,'ROI'),       compileOpts.ROI = 1; end
if ~isfield(compileOpts,'radRange'),  compileOpts.radRange = [0,5]; end
if ~isfield(compileOpts,'minvexpl'),  compileOpts.minvexpl = 0.2; end
if ~isfield(compileOpts,'minSigma'),  compileOpts.minSigma = 0.1; end
if ~isfield(compileOpts,'subList'),   compileOpts.subList = 1:10; end
if ~isfield(compileOpts,'edgeSigma'), compileOpts.edgeSigma = 0; end
validateattributes(compileOpts.ROI,{'numeric'},{'scalar','integer','positive','finite'});
validateattributes(compileOpts.radRange,{'numeric'},{'vector','numel',2,'real','finite'});
if compileOpts.radRange(1) < 0 || compileOpts.radRange(2) <= compileOpts.radRange(1)
    error('compileStimAndSubData:badRadRange','radRange must be [nonnegativeMin largerMax].');
end
validateattributes(compileOpts.minvexpl,{'numeric'},{'scalar','real','finite','nonnegative'});
validateattributes(compileOpts.minSigma,{'numeric'},{'scalar','real','finite','positive'});
validateattributes(compileOpts.edgeSigma,{'numeric'},{'scalar','real','finite','nonnegative'});


% Compile data to generate the 'subjectData' and 'stimData' structures:
% subjectData is a compact cell array (one entry per requested subject).
% Each subject struct has fields:
%
%   .Gprf           [nPix x nVox]
%                   pRF weights for this subject
%
%   .prfXY          [nVox x 2]
%                   pRF centers [x y] in degrees
%
%   .w_vox          [nVox x 1]
%                   voxel weights, e.g. pRF fit quality in [0,1]
%
%   .Yfull{r}       [nT_r x nVox]
%                   measured full-stimulus timecourses for run r
%
%   .Yscot{r}       [nT_r x nVox]
%                   measured scotoma-stimulus timecourses for run r
%
% stimData is a cell array, one entry per subject, because each subject has
% their own fitted HRF. Each entry has:
%
%   .SfullRaw{r}    [nT_r x nPix]
%                   run-specific full-stimulus aperture, UNCONVOLVED
%
%   .SscotRaw{r}    [nT_r x nPix]
%                   run-specific scotoma-stimulus aperture, UNCONVOLVED
%
%   .hrfParams      1-by-2 struct array, ordered {'L','R'} as .hemiOrder
%   .hemiOrder      {'L','R'}
%   .tStim{r}       [1 x nT_r] stimulus time base, seconds
%   .TR             repetition time, seconds
%
% The designs are deliberately NOT convolved here. The HRF differs between
% hemispheres, so a single convolved design cannot be correct for every
% voxel. Callers project the raw design onto a voxel's pRF and then convolve
% that single time course with that voxel's own HRF, selected with
% subjectData.hemIdx (see convHRF and loadHRFParams). For the linear model
% this is exact, because convolution along time commutes with the spatial
% projection; CSS callers must project, apply the exponent, then convolve.

allSubs = [1,2,3,4,5,6,8,9,10,11];
subList = compileOpts.subList(:).';
if any(subList < 1 | subList > numel(allSubs) | subList ~= round(subList))
    error('compileStimAndSubData:badSubList','compileOpts.subList contains an invalid index.');
end
subjectData = cell(numel(subList),1);
stimData = cell(numel(subList),1);
x = [];
y = [];

TR = 1.2; % acquisition repetition time (seconds)

prfType = 'ses-study1_task-logbar';

for si = 1:numel(subList)
    sub = subList(si);
    subNum = allSubs(sub);


    % load prfs from 'study 1'

    prfName = sprintf('data/sub-%02g/sub-%02g_%s_prfs.mat',...
        subNum,subNum,prfType);

    prfFile = load(prfName);
    if ~isfield(prfFile,'prfs')
        error('compileStimAndSubData:missingPRFs','%s does not contain prfs.',prfName);
    end
    allPrfs = prfFile.prfs;

    % HRF parameters, one per hemisphere, ordered {'L','R'}.
    %
    % prfFile.hrfParams is deliberately ignored. It is the RIGHT-hemisphere
    % fit for every subject, because older2/CleanScotomaData.m assigned
    % `hrfParams = tmpprfs.hrfParams` after its {L,R} loop had closed, so the
    % last hemisphere loaded always won. See loadHRFParams for the details
    % and for the published table this reads instead.
    if isfield(compileOpts,'hrfParams')
        supplied = compileOpts.hrfParams;
        if iscell(supplied)
            if ~all(cellfun(@(v) isstruct(v) && isscalar(v),supplied))
                error('compileStimAndSubData:badHRFOverride', ...
                      'Every hrfParams cell must contain one scalar struct.');
            end
            supplied = [supplied{:}];
        end
        if ~isstruct(supplied)
            error('compileStimAndSubData:badHRFOverride', ...
                  'compileOpts.hrfParams must contain HRF parameter structs.');
        end
        if isscalar(supplied)
            hrfParamsThis = [supplied,supplied]; % caller forces one HRF on both
        elseif numel(supplied) == 2
            hrfParamsThis = supplied(:).';
        else
            error('compileStimAndSubData:badHRFOverride', ...
                  'compileOpts.hrfParams must hold one struct or two ({L,R}).');
        end
    else
        hrfParamsThis = loadHRFParams(subNum);
    end
    stimData{si} = struct();
    stimData{si}.hrfParams = hrfParamsThis;
    stimData{si}.hemiOrder = {'L','R'};
    stimData{si}.TR = TR;
    stimData{si}.tStim = cell(1,3);
    stimData{si}.SfullRaw = cell(1,3);
    stimData{si}.SscotRaw = cell(1,3);

    requiredPrfFields = {'varea','x0','y0','sigma','vexpl','hemisphere'};
    if ~all(isfield(allPrfs,requiredPrfFields))
        error('compileStimAndSubData:badPRFs','%s lacks one or more required pRF fields.',prfName);
    end
    x0 = double(allPrfs.x0(:));
    y0 = double(allPrfs.y0(:));
    sigma = double(allPrfs.sigma(:));
    vexpl = double(allPrfs.vexpl(:));
    varea = double(allPrfs.varea(:));
    hemAll = string(allPrfs.hemisphere(:));
    if numel(unique([numel(x0),numel(y0),numel(sigma),numel(vexpl),numel(varea), ...
                     numel(hemAll)])) ~= 1
        error('compileStimAndSubData:prfFieldSizes','pRF fields in %s have different lengths.',prfName);
    end
    if ~all(ismember(upper(strtrim(hemAll)),["L","R"]))
        error('compileStimAndSubData:badHemisphere', ...
              'prfs.hemisphere in %s must contain only L and R labels.',prfName);
    end
    allr = hypot(x0,y0);
    finitePrf = isfinite(x0) & isfinite(y0) & isfinite(sigma) & isfinite(vexpl);
    lowerExtent = allr-compileOpts.edgeSigma*sigma;
    upperExtent = allr+compileOpts.edgeSigma*sigma;
    eligible = varea == compileOpts.ROI & finitePrf & ...
        lowerExtent >= compileOpts.radRange(1) & ...
        upperExtent <= compileOpts.radRange(2) & ...
        vexpl > compileOpts.minvexpl & sigma >= compileOpts.minSigma;
    id = find(eligible);
    subjectData{si}.subNum = subNum;
    subjectData{si}.sourceVoxelIndex = id(:);
    subjectData{si}.selection = struct('ROI',compileOpts.ROI, ...
        'radRange',compileOpts.radRange,'minVarianceExplained',compileOpts.minvexpl, ...
        'minSigma',compileOpts.minSigma,'edgeSigma',compileOpts.edgeSigma, ...
        'nSourceVoxels',numel(x0),'nEligibleVoxels',numel(id));
    subjectData{si}.w_vox = vexpl(id);
    subjectData{si}.prfXY = [x0(id),y0(id)];
    subjectData{si}.sigma = sigma(id);
    % Hemisphere label per retained voxel, and its index into hrfParams.
    subjectData{si}.hemisphere = cellstr(upper(strtrim(hemAll(id))));
    subjectData{si}.hemIdx = uint8(1+strcmp(subjectData{si}.hemisphere,'R'));

    for runNum = 1:3
        fprintf('subject %d, run %d ',subNum,runNum)

        % load the two bold data sets

        boldName = sprintf('data/sub-%02g/sub-%02g_task-%s_run-%g_bold.mat', ...
            subNum,subNum,'scotoma',runNum);
        if ~isfile(boldName), error('compileStimAndSubData:missingBold','Missing %s.',boldName); end
        tmpBold = load(boldName);
        if ~isfield(tmpBold,'bold'), error('compileStimAndSubData:badBold','%s does not contain bold.',boldName); end
        allBold.scotoma = tmpBold.bold;
        boldName = sprintf('data/sub-%02g/sub-%02g_task-%s_run-%g_bold.mat', ...
            subNum,subNum,'logbar',runNum);
        if ~isfile(boldName), error('compileStimAndSubData:missingBold','Missing %s.',boldName); end
        tmpBold = load(boldName);
        if ~isfield(tmpBold,'bold'), error('compileStimAndSubData:badBold','%s does not contain bold.',boldName); end
        allBold.full = tmpBold.bold;

        % load in the full logbar stimulus
        % load 'funcOf' which contains the fields t,x,y

        [stim.full,funcOf] = loadScotomaStimuli(runNum,'logbar',TR);
        [stim.scotoma,funcOfScot] = loadScotomaStimuli(runNum,'scotoma',TR);
        sameGrid = isequal(size(funcOf.x),size(funcOfScot.x)) && ...
                   max(abs(funcOf.x(:)-funcOfScot.x(:))) < 1e-10 && ...
                   max(abs(funcOf.y(:)-funcOfScot.y(:))) < 1e-10;
        sameTime = isequal(size(funcOf.t),size(funcOfScot.t)) && ...
                   max(abs(funcOf.t(:)-funcOfScot.t(:))) < 1e-10;
        if ~sameGrid || ~sameTime || ~isequal(size(stim.full),size(stim.scotoma))
            error('compileStimAndSubData:stimulusMismatch', ...
                  'Full and scotoma stimuli differ in grid or timing for run %d.',runNum);
        end
        if any(stim.scotoma(:) > stim.full(:)+1e-10)
            error('compileStimAndSubData:notMaskedStimulus', ...
                  'The scotoma stimulus is not a masked version of the full stimulus in run %d.',runNum);
        end

        xr = funcOf.x;
        yr = funcOf.y;
        if si == 1 && runNum == 1
            x = xr; y = yr;
        elseif ~isequal(size(x),size(xr)) || ~isequal(size(y),size(yr)) || ...
               max(abs(x(:)-xr(:))) > 1e-10 || max(abs(y(:)-yr(:))) > 1e-10
            error('compileStimAndSubData:gridMismatch', ...
                  'Stimulus grid differs for subject %d, run %d.',subNum,runNum);
        end
        if runNum == 1
            stimData{si}.x = double(xr);
            stimData{si}.y = double(yr);
        end
        nx = size(funcOf.x,2);
        ny = size(funcOf.x,1);
        nt = length(funcOf.t);

        %%
        % select a subset of pRFs

        nVox = length(id);
        fprintf('%d voxels\n',nVox)
        [bold.full,prfs] = subData(allBold.full,allPrfs,id);
        bold.scotoma = subData(allBold.scotoma,allPrfs,id);

        subjectData{si}.Yfull{runNum} = bold.full.data;
        subjectData{si}.Yscot{runNum} = bold.scotoma.data;
        if size(subjectData{si}.Yfull{runNum},1) ~= nt || ...
                size(subjectData{si}.Yscot{runNum},1) ~= nt || ...
                size(subjectData{si}.Yfull{runNum},2) ~= nVox || ...
                size(subjectData{si}.Yscot{runNum},2) ~= nVox
            error('compileStimAndSubData:boldStimulusMismatch', ...
                  'BOLD and comparison-stimulus dimensions differ for subject %d, run %d.', ...
                  subNum,runNum);
        end

        %%
        % Generate the design matrix X

        if runNum == 1
            G = zeros(nx*ny,nVox);
            for i = 1:nVox
                pRF = struct('center',[prfs.x0(i),prfs.y0(i)], ...
                             'sig',prfs.sigma(i),'ar',1);
                G(:,i) = Gauss(pRF,x,y,1);
            end
            % Unit mass makes beta carry response amplitude independently of
            % pRF size and is required by every downstream decomposition.
            area = sum(G,1);
            if any(~isfinite(G(:))) || any(~isfinite(area) | area <= 0)
                error('compileStimAndSubData:badGaussian','A pRF has zero or invalid mass.');
            end
            subjectData{si}.Gprf = G./area;
        end


        % Reshape to the matrix S: columns of S are one pixel's time course,
        % rows are the pixel image at one time point.
        %
        % These are stored UNCONVOLVED. The HRF is hemisphere specific, so no
        % single convolved design is correct for every voxel; the fitting
        % code projects onto each voxel's pRF and convolves that one time
        % course with that voxel's own HRF. Convolution along time commutes
        % with the spatial projection, so the linear predictions are
        % unchanged; CSS predictions must apply the exponent in between.

        stimData{si}.SfullRaw{runNum} = reshape(double(stim.full),[nx*ny,nt])';
        stimData{si}.SscotRaw{runNum} = reshape(double(stim.scotoma),[nx*ny,nt])';
        stimData{si}.tStim{runNum} = double(funcOf.t(:).');

    end
    clear id
end
end
