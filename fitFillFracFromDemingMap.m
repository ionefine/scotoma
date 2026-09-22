function fitOut = fitFillFracFromDemingMap(subjectReal, stimData, fillGrid, noiseSD, simOpts, demingOpts)
% fitFillFracFromDemingMap
%
% Fit a filling fraction by matching voxel-wise Deming slopes between
% real data and simulated data, without eccentricity binning.
%
% INPUTS
%   subjectReal   one subjectData struct
%   stimData      stimulus struct with SfullRaw / SscotRaw (unconvolved)
%   fillGrid      candidate missing-stimulus fractions k
%   noiseSD       scalar noise level for simulation
%   simOpts       options for simulateSubjectDataWithFilling
%   demingOpts    options for voxelDemingSlopeMap
%
% OUTPUT
%   fitOut struct with fields:
%       .fillGrid
%       .sse
%       .bestFillFrac
%       .bestIndex
%       .realSlope
%       .bestSimSlope
%       .simSlope
%       .nVoxUsed
%       .k, .kDeming           aliases for .bestFillFrac
%       .massIn                fraction of pRF mass inside the scotoma
%       .r                     slope - (1-m), the mass approximation
%       .rExact                slope - exact no-filling model slope
%       .kVoxel                r/m, reported only where m > 0.1
%       .commonVoxels          voxels used at every grid point

    if nargin < 4 || isempty(noiseSD), noiseSD = 0; end
    if nargin < 5 || isempty(simOpts), simOpts = struct(); end
    if nargin < 6 || isempty(demingOpts), demingOpts = struct(); end
    if isfield(simOpts,'fillMeasure') && ~strcmpi(simOpts.fillMeasure,'k')
        error('fitFillFracFromDemingMap:wrongFillMeasure', ...
              'This function estimates k; simOpts.fillMeasure must be ''k''.');
    end
    simOpts.fillMeasure = 'k';
    fillGrid = fillGrid(:);
    if isempty(fillGrid) || any(~isfinite(fillGrid))
        error('fitFillFracFromDemingMap:badGrid','fillGrid must contain finite values.');
    end
    if ~isscalar(noiseSD) || ~isfinite(noiseSD) || noiseSD < 0
        error('fitFillFracFromDemingMap:badNoise','noiseSD must be a nonnegative scalar.');
    end
    if noiseSD > 0 && (~isfield(simOpts,'rngSeed') || isempty(simOpts.rngSeed))
        error('fitFillFracFromDemingMap:stochasticObjective', ...
              'Set simOpts.rngSeed when noiseSD>0 so every grid point uses matched noise.');
    end


    % real voxel slopes
    realSlope = voxelDemingSlopeMap(subjectReal, demingOpts);

    nF = numel(fillGrid);
    nVox = numel(realSlope);
    simSlope = nan(nVox, nF);

    for k=1:nF
        f = fillGrid(k);

        simData = simulateSubjectDataWithFilling({subjectReal}, stimData, noiseSD, f, simOpts);
        tmpSlope = voxelDemingSlopeMap(simData{1}, demingOpts);
        
        simSlope(:,k) = tmpSlope;

        fprintf('%d of %d, f = %g\n',k,nF,f)

    end

    % Every grid point must be compared on the same voxel population.
    common = isfinite(realSlope(:)) & all(isfinite(simSlope),2);
    if nnz(common) < 3
        error('fitFillFracFromDemingMap:tooFewVoxels', ...
              'Only %d voxels have finite slopes across the entire grid.',nnz(common));
    end
    d = simSlope(common,:) - realSlope(common);
    sse = sum(d.^2,1).';
    nVoxUsed = repmat(nnz(common),nF,1);

    [~, idx] = min(sse);

    fitOut = struct();
    fitOut.fillGrid = fillGrid;
    fitOut.sse = sse;
    fitOut.bestFillFrac = fillGrid(idx);
    fitOut.k = fitOut.bestFillFrac;
    fitOut.kDeming = fitOut.bestFillFrac;
    fitOut.bestIndex = idx;

    fitOut.realSlope = realSlope;
    fitOut.bestSimSlope = simSlope(:,idx);
    fitOut.simSlope = simSlope;
    fitOut.nVoxUsed = nVoxUsed;
    fitOut.commonVoxels = common;

    % Geometry-only response excess r = slope-(1-m). This is stable near
    % m=0, but r=k*m is only an approximation for Deming slopes.
    if isempty(stimData.SfullRaw) || numel(stimData.SfullRaw) ~= numel(stimData.SscotRaw)
        error('fitFillFracFromDemingMap:runMismatch', ...
              'SfullRaw and SscotRaw must contain the same nonzero number of runs.');
    end
    dPix = false(1,size(stimData.SfullRaw{1},2));
    for run = 1:numel(stimData.SfullRaw)
        if ~isequal(size(stimData.SfullRaw{run}),size(stimData.SscotRaw{run})) || ...
           size(stimData.SfullRaw{run},2) ~= size(subjectReal.Gprf,1)
            error('fitFillFracFromDemingMap:dimensionMismatch', ...
                  'Stimulus/G dimensions disagree in run %d.',run);
        end
        dPix = dPix | any(abs(double(stimData.SfullRaw{run}) - ...
                 double(stimData.SscotRaw{run})) > 1e-8,1);
    end
    Gn = double(subjectReal.Gprf);
    gMass = sum(Gn,1);
    if any(~isfinite(Gn(:))) || any(~isfinite(gMass) | gMass <= 0)
        error('fitFillFracFromDemingMap:badGaussian', ...
              'Gprf contains nonfinite values or nonpositive pRF mass.');
    end
    Gn = Gn ./ gMass;
    fitOut.massIn = sum(Gn(dPix,:),1).';
    fitOut.r = realSlope(:) - (1-fitOut.massIn);
    fitOut.kVoxel = nan(size(fitOut.r));
    q = isfinite(fitOut.r) & fitOut.massIn > 0.1;
    fitOut.kVoxel(q) = fitOut.r(q)./fitOut.massIn(q);

    % Exact response excess relative to this forward model at k=0. Use a
    % noiseless prediction even if the legacy grid search included noise.
    [d0,idx0] = min(abs(fillGrid));
    if d0 <= 10*eps(max(1,max(abs(fillGrid)))) && noiseSD == 0
        fitOut.noFillSlope = simSlope(:,idx0);
    else
        sim0 = simulateSubjectDataWithFilling({subjectReal},stimData,0,0,simOpts);
        fitOut.noFillSlope = voxelDemingSlopeMap(sim0{1},demingOpts);
    end
    fitOut.rExact = realSlope(:) - fitOut.noFillSlope(:);
end
