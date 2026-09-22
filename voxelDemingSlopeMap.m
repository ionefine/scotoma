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
%   .lambda             error-variance ratio sigma_Y^2/sigma_X^2 (default 1)
%   .forceZeroIntercept default true
%   .centerData         default true
%   .combineRuns        'concat' or 'mean' (default 'concat')
%   .minDen             default 1e-12
%   .minPairs           minimum finite X/Y pairs (default 3)
%
% OUTPUT
%   beta [nVox x 1]

    if nargin < 2 || isempty(opts), opts = struct(); end
    if ~isstruct(opts) || ~isscalar(opts)
        error('voxelDemingSlopeMap:badOptions','opts must be a scalar struct.');
    end
    if ~isfield(opts, 'lambda'),             opts.lambda = 1; end
    if ~isfield(opts, 'forceZeroIntercept'), opts.forceZeroIntercept = true; end
    if ~isfield(opts, 'centerData'),         opts.centerData = true; end
    if ~isfield(opts, 'combineRuns'),        opts.combineRuns = 'concat'; end
    if ~isfield(opts, 'minDen'),             opts.minDen = 1e-12; end
    if ~isfield(opts, 'minPairs'),           opts.minPairs = 3; end
    if ~isnumeric(opts.lambda) || ~isreal(opts.lambda) || ...
       ~isscalar(opts.lambda) || ~isfinite(opts.lambda) || opts.lambda <= 0
        error('voxelDemingSlopeMap:badLambda','opts.lambda must be positive.');
    end
    if ~isnumeric(opts.minPairs) || ~isreal(opts.minPairs) || ...
       ~isscalar(opts.minPairs) || ~isfinite(opts.minPairs) || opts.minPairs < 2
        error('voxelDemingSlopeMap:badMinPairs','opts.minPairs must be at least 2.');
    end
    if opts.minPairs ~= round(opts.minPairs)
        error('voxelDemingSlopeMap:badMinPairs','opts.minPairs must be an integer.');
    end
    if ~isnumeric(opts.minDen) || ~isreal(opts.minDen) || ...
       ~isscalar(opts.minDen) || ~isfinite(opts.minDen) || opts.minDen < 0
        error('voxelDemingSlopeMap:badMinDen','opts.minDen must be nonnegative.');
    end
    if ~(islogical(opts.forceZeroIntercept) || isnumeric(opts.forceZeroIntercept)) || ...
       ~isscalar(opts.forceZeroIntercept) || ~ismember(opts.forceZeroIntercept,[0 1]) || ...
       ~(islogical(opts.centerData) || isnumeric(opts.centerData)) || ...
       ~isscalar(opts.centerData) || ~ismember(opts.centerData,[0 1])
        error('voxelDemingSlopeMap:badLogicalOption', ...
              'forceZeroIntercept and centerData must be scalar logical values.');
    end
    if ~(ischar(opts.combineRuns) || (isstring(opts.combineRuns) && isscalar(opts.combineRuns)))
        error('voxelDemingSlopeMap:badCombineRuns', ...
              'opts.combineRuns must be ''concat'' or ''mean''.');
    end
    combineRuns = lower(char(opts.combineRuns));

    if ~isstruct(subjectData) || ~isscalar(subjectData) || ...
       ~all(isfield(subjectData,{'Yfull','Yscot'})) || ...
       ~iscell(subjectData.Yfull) || ~iscell(subjectData.Yscot)
        error('voxelDemingSlopeMap:badSubjectData', ...
              'subjectData must contain cell arrays Yfull and Yscot.');
    end
    nRun = numel(subjectData.Yfull);
    if nRun == 0 || numel(subjectData.Yscot) ~= nRun
        error('Yfull and Yscot must have same number of runs.');
    end
    nVox = size(subjectData.Yfull{1},2);
    for r = 1:nRun
        if ~isnumeric(subjectData.Yfull{r}) || ~isreal(subjectData.Yfull{r}) || ...
           ~isnumeric(subjectData.Yscot{r}) || ~isreal(subjectData.Yscot{r}) || ...
           ~ismatrix(subjectData.Yfull{r}) || ~ismatrix(subjectData.Yscot{r}) || ...
           ~isequal(size(subjectData.Yfull{r}),size(subjectData.Yscot{r})) || ...
           size(subjectData.Yfull{r},2) ~= nVox
            error('voxelDemingSlopeMap:dimensionMismatch', ...
                  'Yfull and Yscot dimensions disagree in run %d.',r);
        end
    end

    switch combineRuns
        case 'concat'
            Xall = [];
            Yall = [];

            for r = 1:nRun
                X = double(subjectData.Yfull{r});
                Y = double(subjectData.Yscot{r});

                good = isfinite(X) & isfinite(Y);
                X(~good) = NaN;
                Y(~good) = NaN;

                if opts.centerData
                    X = X - mean(X,1,'omitnan');
                    Y = Y - mean(Y,1,'omitnan');
                end

                Xall = [Xall; X]; %#ok<AGROW>
                Yall = [Yall; Y]; %#ok<AGROW>
            end

            beta = localDemingColumns(Xall, Yall, opts.lambda, ...
                opts.forceZeroIntercept, opts.minDen, opts.minPairs)';

        case 'mean'
            betaRuns = nan(nRun,nVox);

            for r = 1:nRun
                X = double(subjectData.Yfull{r});
                Y = double(subjectData.Yscot{r});

                good = isfinite(X) & isfinite(Y);
                X(~good) = NaN;
                Y(~good) = NaN;

                if opts.centerData
                    X = X - mean(X,1,'omitnan');
                    Y = Y - mean(Y,1,'omitnan');
                end

                betaRuns(r,:) = localDemingColumns(X, Y, opts.lambda, ...
                    opts.forceZeroIntercept, opts.minDen, opts.minPairs);
            end

            beta = mean(betaRuns, 1, 'omitnan')';

        otherwise
            error('voxelDemingSlopeMap:badCombineRuns', ...
                  'opts.combineRuns must be ''concat'' or ''mean''.');
    end
end


function beta = localDemingColumns(X, Y, lambda, forceZeroIntercept, minDen, minPairs)
% Columnwise Deming slopes for paired matrices X and Y.
% X, Y are [nT x nVox], beta is [1 x nVox]

    if forceZeroIntercept
        sx2 = mean(X.^2,1,'omitnan');
        sy2 = mean(Y.^2,1,'omitnan');
        sxy = mean(X.*Y,1,'omitnan');
    else
        mx = mean(X,1,'omitnan');
        my = mean(Y,1,'omitnan');
        Xc = X - mx;
        Yc = Y - my;
        sx2 = mean(Xc.^2,1,'omitnan');
        sy2 = mean(Yc.^2,1,'omitnan');
        sxy = mean(Xc.*Yc,1,'omitnan');
    end

    a = sy2 - lambda .* sx2;
    disc = a.^2 + 4 .* lambda .* (sxy.^2);

    beta = nan(1, size(X,2));
    n = sum(isfinite(X) & isfinite(Y),1);
    ok = n >= minPairs & isfinite(sx2) & isfinite(sy2) & isfinite(sxy) & ...
         (abs(sxy) > minDen) & isfinite(disc) & (disc >= 0);
    root = sqrt(disc);
    direct = ok & a >= 0;
    stable = ok & a < 0;
    beta(direct) = (a(direct) + root(direct))./(2.*sxy(direct));
    beta(stable) = (2.*lambda.*sxy(stable))./(root(stable)-a(stable));

    % A flat Y with variable X has slope zero. If both are flat, the slope
    % is undefined and remains NaN.
    flatY = n >= minPairs & isfinite(sx2) & isfinite(sy2) & ...
            sx2 > minDen & sy2 <= minDen;
    beta(flatY) = 0;
end
