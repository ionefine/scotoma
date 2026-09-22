function subjectDataSim = simulateSubjectDataWithFilling(subjectData, stimData, noiseSD, fillFrac, opts)
% simulateSubjectDataWithFilling
%
% Simulate either of the two filling-in parameterizations:
%
%   k:  Ssim = Sscot + fillFrac*(Sfull-Sscot)
%   k2: Ssim = Sscot + fillFrac*Sfull
%
% For k, fillFrac=0 is no filling-in and fillFrac=1 restores the full
% stimulus. For k2, fillFrac is the fraction of the full-field stimulus
% added to the scotoma stimulus; outside the scotoma it is a gain change.
%
% Returns a structure identical to subjectData in shape.
%
% INPUTS
%   subjectData   cell array of subject structs
%   stimData      cell array matching subjectData; each entry contains
%                 Sfull{r}, Sscot{r}. A single struct is accepted for one
%                 subject.
%   noiseSD       scalar noise SD
%   fillFrac      finite scalar (normally in [0,1])
%
% opts fields:
%   .useVoxelGain   default false
%   .gainRange      default [0.5 1.5]
%   .rngSeed        default []
%   .verbose        default true
%   .zeroMean       default true
%   .fillMeasure    'k' (default) or 'k2'
%
% OUTPUT
%   subjectDataSim  cell array matching subjectData, but with simulated
%                   Yfull and Yscot fields

    if nargin < 5 || isempty(opts), opts = struct(); end
    if ~isfield(opts,'useVoxelGain'), opts.useVoxelGain = false; end
    if ~isfield(opts,'gainRange'),    opts.gainRange = [0.5 1.5]; end
    if ~isfield(opts,'rngSeed'),      opts.rngSeed = []; end
    if ~isfield(opts,'verbose'),      opts.verbose = true; end
    if ~isfield(opts,'zeroMean'),     opts.zeroMean = true; end
    if ~isfield(opts,'fillMeasure'),  opts.fillMeasure = 'k'; end

    if ~iscell(subjectData), subjectData = {subjectData}; end
    if isempty(subjectData)
        error('simulateSubjectDataWithFilling:noSubjects','subjectData is empty.');
    end

    if ~isscalar(fillFrac) || ~isfinite(fillFrac)
        error('simulateSubjectDataWithFilling:badFillFrac', ...
              'fillFrac must be a finite scalar.');
    end
    if ~any(strcmpi(opts.fillMeasure,{'k','k2'}))
        error('simulateSubjectDataWithFilling:badFillMeasure', ...
              'opts.fillMeasure must be ''k'' or ''k2''.');
    end
    if ~isscalar(noiseSD) || ~isfinite(noiseSD) || noiseSD < 0
        error('simulateSubjectDataWithFilling:badNoise','noiseSD must be nonnegative.');
    end
    if numel(opts.gainRange) ~= 2 || any(~isfinite(opts.gainRange)) || ...
       opts.gainRange(2) < opts.gainRange(1)
        error('simulateSubjectDataWithFilling:badGainRange','opts.gainRange must be [min max].');
    end
    if ~iscell(stimData), stimData = {stimData}; end

    if ~isempty(opts.rngSeed)
        oldRng = rng;
        restoreRng = onCleanup(@() rng(oldRng)); %#ok<NASGU>
        rng(opts.rngSeed);
    end

    nSub = numel(subjectData);
    if numel(stimData) ~= nSub
        error('simulateSubjectDataWithFilling:stimulusCount', ...
              'stimData must contain one entry per subject.');
    end

    subjectDataSim = subjectData;

    for s = 1:nSub
        if opts.verbose
            fprintf('Simulating subject %d / %d, fillFrac = %.3f\n', s, nSub, fillFrac);
        end

        G = double(subjectData{s}.Gprf);
        gMass = sum(G,1);
        if any(~isfinite(G(:))) || any(~isfinite(gMass) | gMass <= 0)
            error('simulateSubjectDataWithFilling:badGaussian', ...
                  'Subject %d has nonfinite Gprf values or nonpositive pRF mass.',s);
        end
        G = G ./ gMass;
        [~, nVox] = size(G);
        nRun = numel(stimData{s}.SfullRaw);
        if nRun == 0 || numel(stimData{s}.SscotRaw) ~= nRun
            error('simulateSubjectDataWithFilling:runMismatch', ...
                  'SfullRaw and SscotRaw run counts disagree for subject %d.',s);
        end
        % One HRF per hemisphere; the stored designs are unconvolved.
        nHem = numel(stimData{s}.hrfParams);
        hemIdx = double(subjectData{s}.hemIdx(:).');
        if numel(hemIdx) ~= nVox || any(hemIdx < 1 | hemIdx > nHem)
            error('simulateSubjectDataWithFilling:badHemIdx', ...
                  'hemIdx for subject %d is missing or indexes outside hrfParams.',s);
        end
        hrf = cell(nHem,1);
        for h = 1:nHem
            hrf{h} = cell(nRun,1);
            for r = 1:nRun
                hh = double(hrf_twogamma(stimData{s}.hrfParams(h),stimData{s}.tStim{r}));
                hrf{h}{r} = hh(:);
            end
        end

        if opts.useVoxelGain
            gmin = opts.gainRange(1);
            gmax = opts.gainRange(2);
            voxelGain = gmin + (gmax - gmin) * rand(nVox,1);
        else
            voxelGain = ones(nVox,1);
        end

        for r = 1:nRun
            Sfull = double(stimData{s}.SfullRaw{r});
            Sscot = double(stimData{s}.SscotRaw{r});
            if ~isequal(size(Sfull),size(Sscot)) || size(Sfull,2) ~= size(G,1)
                error('simulateSubjectDataWithFilling:dimensionMismatch', ...
                      'Stimulus/G dimensions disagree for subject %d, run %d.',s,r);
            end

            switch lower(opts.fillMeasure)
                case 'k'
                    Ssim = Sscot + fillFrac*(Sfull-Sscot);
                case 'k2'
                    Ssim = Sscot + fillFrac*Sfull;
            end

            % Forward model: project onto each pRF, then convolve with that
            % voxel's own hemisphere HRF.
            Yfull = convByHemisphere(Sfull * G,hrf,stimData{s}.TR,hemIdx,r);
            Ysim  = convByHemisphere(Ssim  * G,hrf,stimData{s}.TR,hemIdx,r);

            % Optional voxel gain
            Yfull = Yfull .* voxelGain';
            Ysim  = Ysim  .* voxelGain';

            % Add independent Gaussian noise
            if noiseSD > 0
                Yfull = Yfull + noiseSD * randn(size(Yfull));
                Ysim  = Ysim  + noiseSD * randn(size(Ysim));
            end

            % Zero-mean each voxel timecourse within run
            if opts.zeroMean
                Yfull = Yfull - mean(Yfull, 1);
                Ysim  = Ysim  - mean(Ysim, 1);
            end

            subjectDataSim{s}.Yfull{r} = Yfull;
            subjectDataSim{s}.Yscot{r} = Ysim;
        end
    end
end
