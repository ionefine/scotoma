function subjectDataSim = simulateSubjectDataWithFilling(subjectData, stimData, noiseSD, fillFrac, opts)
% simulateSubjectDataWithFilling
%
% Simulate data with a filling-in fraction:
%
%   Ssim = (1-fillFrac)*Sscot + fillFrac*Sfull
%
% fillFrac = 0  -> no filling-in
% fillFrac = 1  -> full stimulus restored
%
% Returns a structure identical to subjectData in shape.
%
% INPUTS
%   subjectData   cell array of subject structs
%   stimData      struct with Sfull{r}, Sscot{r}
%   noiseSD       scalar noise SD
%   fillFrac      scalar in [0,1]
%
% opts fields:
%   .useVoxelGain   default false
%   .gainRange      default [0.5 1.5]
%   .rngSeed        default []
%   .verbose        default true
%   .zeroMean       default true
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

    if ~isempty(opts.rngSeed)
        rng(opts.rngSeed);
    end

    nSub = numel(subjectData);
    nRun = numel(stimData.Sfull);

    subjectDataSim = subjectData;

    for s = 1:nSub
        if opts.verbose
            fprintf('Simulating subject %d / %d, fillFrac = %.3f\n', s, nSub, fillFrac);
        end

        G = double(subjectData{s}.Gprf);
        [~, nVox] = size(G);

        if opts.useVoxelGain
            gmin = opts.gainRange(1);
            gmax = opts.gainRange(2);
            voxelGain = gmin + (gmax - gmin) * rand(nVox,1);
        else
            voxelGain = ones(nVox,1);
        end

        for r = 1:nRun
            Sfull = double(stimData.Sfull{r});
            Sscot = double(stimData.Sscot{r});

            % Mixture stimulus
            Ssim = (1 - fillFrac) * Sscot + fillFrac * Sfull;

            % Forward model
            Yfull = Sfull * G;
            Ysim  = Ssim  * G;

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