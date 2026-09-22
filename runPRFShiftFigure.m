% runPRFShiftFigure  Produce the pRF shift figure, start to finish.
%
%   prfShiftFit     one subject -> per-voxel cross-products   (geometry)
%   prfShiftPool    cross-products -> parameters + interval   (statistics)
%   prfShiftFigure  parameters -> figure                      (presentation)
% The returned table reports fixed-pRF k and k2 as the primary coefficients,
% plus k_joint and k2_joint from the shift/size models. All hold full-field
% beta fixed. Only fixed-pRF k2 reduces directly to g-1 outside the scotoma.
%
% Needs compileStimAndSubData and Gauss on the path.
% This is the slow, iterated GEOMETRY analysis. For the fixed-pRF k
% response decomposition, use runKResponseFigure instead.

clear subjectData                       % force a recompile if ROI changed

compileOpts.ROI       = 1;              % 1 = V1, 2 = V2, 3 = V3
compileOpts.radRange  = [0, 4];
compileOpts.minvexpl  = 0.2;
compileOpts.minSigma  = 0.1;
compileOpts.subList   = 1:10;
compileOpts.edgeSigma = 0;

[subjectData, stimData, x, y] = compileStimAndSubData(compileOpts);

% Iterated fit. maxIter = 0 gives the original single-pass estimate, which
% saturates dsigma near 0.45 deg; 3 committed updates remove that.
src = struct('subjectData', {subjectData}, 'stimData', {stimData}, ...
             'x', x, 'y', y);
iterOpts = struct('maxIter', 3, 'tol', 0.01, 'damp', 1.0, 'maxStep', 0.5, ...
                  'useParallel', true, ...
                  'fitOpts', struct('minEcc', 0.1, 'verbose', false));

% Start one pool and retain it for every re-linearisation pass. If a pool
% already exists, its current size and profile are left unchanged.
if iterOpts.useParallel && ~isempty(ver('parallel')) && ...
   license('test','Distrib_Computing_Toolbox') && isempty(gcp('nocreate'))
    try
        parpool('local');
    catch ME
        warning('runPRFShiftFigure:poolFailed', ...
                'Could not start a parallel pool (%s). Using serial fitting.',ME.message);
        iterOpts.useParallel = false;
    end
end

T = prfShiftFigure(src, [1.25 3], 0.25, 2000, iterOpts);
disp(T)

% Response-space view of the fixed-pRF k estimate. The band between the
% orange and blue curves is r=k*m; the gap from the blue curve to 1 is the
% reduction in BOLD response relative to the full-field condition.
[responseFigure,responseComponents] = plotKResponseDecomposition(T);
disp(responseComponents)
