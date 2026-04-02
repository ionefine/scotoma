function fitOut = fitFillFracFromVoxelDeming(subjectReal, stimData, fillGrid, noiseSD, simOpts, demingOpts)
% fitFillFracFromVoxelDeming
%
% Fit a filling fraction by matching voxel-wise Deming slopes between
% real data and simulated data, without eccentricity binning.
%
% INPUTS
%   subjectReal   one subjectData struct
%   stimData      stimulus struct with Sfull / Sscot
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

    if nargin < 4 || isempty(simOpts), simOpts = struct(); end
    if nargin < 5 || isempty(demingOpts), demingOpts = struct(); end


    % real voxel slopes
    realSlope = voxelDemingSlopeMap(subjectReal, demingOpts);

    nF = numel(fillGrid);
    nVox = numel(realSlope);
    simSlope = nan(nVox, nF);
    sse = nan(nF,1);
    nVoxUsed = nan(nF,1);

    for k=1:nF
        f = fillGrid(k);

        simData = simulateSubjectDataWithFilling({subjectReal}, stimData, noiseSD, f, simOpts);
        tmpSlope = voxelDemingSlopeMap(simData{1}, demingOpts);
        
        simSlope(:,k) = tmpSlope;

        ok = isfinite(realSlope) & isfinite(tmpSlope);
        d = realSlope(ok) - tmpSlope(ok);

        sse(k) = sum(d.^2);


        fprintf('%d of %d, f = %g, d = %0.2f\n',k,nF,f,sse(k))
        nVoxUsed(k) = nnz(ok);

    end

    [~, idx] = min(sse);

    fitOut = struct();
    fitOut.fillGrid = fillGrid(:);
    fitOut.sse = sse;
    fitOut.bestFillFrac = fillGrid(idx);
    fitOut.bestIndex = idx;

    fitOut.realSlope = realSlope;
    fitOut.bestSimSlope = simSlope(:,idx);
    fitOut.simSlope = simSlope;
    fitOut.nVoxUsed = nVoxUsed;
end