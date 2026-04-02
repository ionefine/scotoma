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

                X(~isfinite(X)) = 0;
                Y(~isfinite(Y)) = 0;

                if opts.centerData
                    X = X - mean(X, 1);
                    Y = Y - mean(Y, 1);
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

                X(~isfinite(X)) = 0;
                Y(~isfinite(Y)) = 0;

                if opts.centerData
                    X = X - mean(X, 1);
                    Y = Y - mean(Y, 1);
                end

                betaRuns(r,:) = localDemingColumns(X, Y, opts.lambda, ...
                    opts.forceZeroIntercept, opts.minDen);
            end

            beta = mean(betaRuns, 1, 'omitnan')';
            beta(isnan(beta)) = 0; % deal with zero slope

        otherwise
            error('opts.combineRuns must be ''concat'' or ''mean''.');
    end
end


function beta = localDemingColumns(X, Y, lambda, forceZeroIntercept, minDen)
% Columnwise Deming slopes for paired matrices X and Y.
% X, Y are [nT x nVox], beta is [1 x nVox]

    if forceZeroIntercept
        sx2 = mean(X.^2, 1);
        sy2 = mean(Y.^2, 1);
        sxy = mean(X .* Y, 1);
    else
        mx = mean(X,1);
        my = mean(Y,1);
        Xc = X - mx;
        Yc = Y - my;
        sx2 = mean(Xc.^2, 1);
        sy2 = mean(Yc.^2, 1);
        sxy = mean(Xc .* Yc, 1);
    end

    disc = (sy2 - lambda .* sx2).^2 + 4 .* lambda .* (sxy.^2);

    beta = nan(1, size(X,2));
    ok = isfinite(sx2) & isfinite(sy2) & isfinite(sxy) & ...
         (abs(sxy) > minDen) & isfinite(disc) & (disc >= 0);

    beta(ok) = (sy2(ok) - lambda .* sx2(ok) + sqrt(disc(ok))) ./ ...
               (2 .* sxy(ok));

    beta(sy2<minDen) = 0;  % set slope to zero if no variabity in y.
end