function fitOut = fitFillFracFromResidualDeming(subjectReal, stimData, fillGrid, noiseSD, simOpts, demingOpts)
% fitFillFracFromResidualDeming
%
% Fit a filling fraction by matching FEEDFORWARD-RESIDUAL voxel-wise
% Deming slopes between real and simulated data.
%
% Compared to fitFillFracFromDemingMap, the pure-feedforward prediction
% (Sscot*G, i.e. simulated Yscot at fillFrac = 0) is subtracted from
% both real and simulated Yscot before computing slopes against Yfull.
% Filling-in is therefore fit only to the part of the response that
% feedforward cannot explain: voxels whose pRF mostly extends outside
% the scotoma (and whose feedforward already predicts the data) have
% residual slope ~ 0 in both real and sim, so they no longer pull the
% fit toward fillFrac = 1.
%
% INPUTS
%   subjectReal   one subjectData struct
%   stimData      stimulus struct with Sfull / Sscot
%   fillGrid      vector of fillFrac values to search
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
%       .resRealSlope     deming slope of (Yscot - Sscot*G) vs Yfull
%       .resBestSimSlope  same for simulated data at bestFillFrac
%       .resSimSlope      [nVox x nF] all sim residual slopes
%       .realSlope        original (non-residual) deming slope, for plotting
%       .nVoxUsed

    if nargin < 5 || isempty(simOpts), simOpts = struct(); end
    if nargin < 6 || isempty(demingOpts), demingOpts = struct(); end

    nRun = numel(subjectReal.Yscot);

    % feedforward-only prediction: Yscot_ff = Sscot * G
    simFFOpts = simOpts;
    simFFOpts.verbose = false;
    simFF = simulateSubjectDataWithFilling({subjectReal}, stimData, 0, 0, simFFOpts);

    % residualised real subject struct
    subjectResReal = subjectReal;
    for r = 1:nRun
        subjectResReal.Yscot{r} = subjectReal.Yscot{r} - simFF{1}.Yscot{r};
    end

    resRealSlope = voxelDemingSlopeMap(subjectResReal, demingOpts);
    realSlope    = voxelDemingSlopeMap(subjectReal,    demingOpts);

    nF   = numel(fillGrid);
    nVox = numel(resRealSlope);
    resSimSlope = nan(nVox, nF);
    sse         = nan(nF, 1);
    nVoxUsed    = nan(nF, 1);

    simOptsLocal = simOpts;
    simOptsLocal.verbose = false;

    for k = 1:nF
        f = fillGrid(k);

        simData = simulateSubjectDataWithFilling({subjectReal}, stimData, noiseSD, f, simOptsLocal);

        subjectResSim = simData{1};
        for r = 1:nRun
            subjectResSim.Yscot{r} = simData{1}.Yscot{r} - simFF{1}.Yscot{r};
        end

        tmpSlope = voxelDemingSlopeMap(subjectResSim, demingOpts);
        resSimSlope(:,k) = tmpSlope;

        ok = isfinite(resRealSlope) & isfinite(tmpSlope);
        d  = resRealSlope(ok) - tmpSlope(ok);
        sse(k)      = sum(d.^2);
        nVoxUsed(k) = nnz(ok);

        fprintf('%d of %d, f = %g, sse = %0.2f\n', k, nF, f, sse(k));
    end

    [~, idx] = min(sse);

    fitOut = struct();
    fitOut.fillGrid        = fillGrid(:);
    fitOut.sse             = sse;
    fitOut.bestFillFrac    = fillGrid(idx);
    fitOut.bestIndex       = idx;
    fitOut.resRealSlope    = resRealSlope;
    fitOut.resBestSimSlope = resSimSlope(:,idx);
    fitOut.resSimSlope     = resSimSlope;
    fitOut.realSlope       = realSlope;
    fitOut.nVoxUsed        = nVoxUsed;
end
