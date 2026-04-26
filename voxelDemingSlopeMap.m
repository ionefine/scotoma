% DemingSlopeByVoxelStuff.m

function beta = voxelDemingSlopeMap(subjectData, opts)
% voxelDemingSlopeMap
%
% Compute Deming-regression slopes for every voxel comparing
% scotoma vs full fMRI timecourses.
%
% INPUT
%   subjectData.Yfull{r}   [nT_r x nVox]
%   subjectData.Yscot{r}   [nT_r x nVox]
%
% opts fields:
%   .lambda             default 1
%   .forceZeroIntercept default true
%   .centerData         default true
%   .combineRuns        'concat' or 'mean' (default 'concat')
%   .minDen             default 1e-12
%
% OUTPUT
%   beta [nVox x 1]
%
% NaN handling: samples where Yfull or Yscot is non-finite are excluded
% per voxel (per column) from the slope calculation. Centering uses
% finite samples only. Voxels with fewer than 2 valid samples produce
% NaN slopes; this is consistent across both 'concat' and 'mean'
% combineRuns modes.

    if nargin < 2 || isempty(opts), opts = struct(); end
    if ~isfield(opts, 'lambda'),             opts.lambda = 1; end
    if ~isfield(opts, 'forceZeroIntercept'), opts.forceZeroIntercept = true; end
    if ~isfield(opts, 'centerData'),         opts.centerData = true; end
    if ~isfield(opts, 'combineRuns'),        opts.combineRuns = 'concat'; end
    if ~isfield(opts, 'minDen'),             opts.minDen = 1e-12; end

    nRun = numel(subjectData.Yfull);
    if numel(subjectData.Yscot) ~= nRun
        error('Yfull and Yscot must have same number of runs.');
    end

    switch lower(opts.combineRuns)
        case 'concat'
            Xall = [];
            Yall = [];

            for r = 1:nRun
                X = double(subjectData.Yfull{r});
                Y = double(subjectData.Yscot{r});

                if opts.centerData
                    X = X - mean(X, 1, 'omitnan');
                    Y = Y - mean(Y, 1, 'omitnan');
                end

                Xall = [Xall; X]; %#ok<AGROW>
                Yall = [Yall; Y]; %#ok<AGROW>
            end

            beta = localDemingColumns(Xall, Yall, opts.lambda, ...
                opts.forceZeroIntercept, opts.minDen)';

        case 'mean'
            betaRuns = nan(nRun, size(subjectData.Yfull{1},2));

            for r = 1:nRun
                X = double(subjectData.Yfull{r});
                Y = double(subjectData.Yscot{r});

                if opts.centerData
                    X = X - mean(X, 1, 'omitnan');
                    Y = Y - mean(Y, 1, 'omitnan');
                end

                betaRuns(r,:) = localDemingColumns(X, Y, opts.lambda, ...
                    opts.forceZeroIntercept, opts.minDen);
            end

            beta = mean(betaRuns, 1, 'omitnan')';

        otherwise
            error('opts.combineRuns must be ''concat'' or ''mean''.');
    end
end


function beta = localDemingColumns(X, Y, lambda, forceZeroIntercept, minDen)
% Columnwise Deming slopes for paired matrices X and Y.
% NaN/Inf samples are excluded per column (per voxel).
% X, Y are [nT x nVox], beta is [1 x nVox]

    mask = isfinite(X) & isfinite(Y);
    X(~mask) = 0;
    Y(~mask) = 0;
    nValid = sum(mask, 1);
    den    = max(nValid, 1);  % avoid divide-by-zero; voxel filtered below

    if forceZeroIntercept
        sx2 = sum(X.^2,    1) ./ den;
        sy2 = sum(Y.^2,    1) ./ den;
        sxy = sum(X .* Y,  1) ./ den;
    else
        mx = sum(X, 1) ./ den;
        my = sum(Y, 1) ./ den;
        Xc = (X - mx) .* mask;   % invalid samples stay at 0
        Yc = (Y - my) .* mask;
        sx2 = sum(Xc.^2,    1) ./ den;
        sy2 = sum(Yc.^2,    1) ./ den;
        sxy = sum(Xc .* Yc, 1) ./ den;
    end

    disc = (sy2 - lambda .* sx2).^2 + 4 .* lambda .* (sxy.^2);

    beta = nan(1, size(X,2));
    ok = isfinite(sx2) & isfinite(sy2) & isfinite(sxy) & ...
         (abs(sxy) > minDen) & isfinite(disc) & (disc >= 0) & ...
         (nValid >= 2);

    beta(ok) = (sy2(ok) - lambda .* sx2(ok) + sqrt(disc(ok))) ./ ...
               (2 .* sxy(ok));

    beta(sy2 < minDen) = 0;  % no variability in y -> slope = 0
end
