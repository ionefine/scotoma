function [b,nVox,sse,dof] = prfShiftPool(fits,cols,mask,fillMeasure)
% prfShiftPool  Pool prfShiftFit cross-products and solve one model.
%
% Cross-products are summed over the selected voxels before the model is
% solved. This preserves the intended joint fit across voxels and avoids
% averaging unstable voxel-wise ratios.
%
% MODEL SELECTION
%   cols = 1:3   radial shift, tangential shift, and sigma change
%   cols = 4     fixed-pRF filling-in model
%   cols = 1:4   shifts and filling-in fitted jointly
%
% INPUTS
%   fits         prfShiftFit output, or a cell array with one per subject
%   cols         unique parameter columns selected from 1:4
%   mask         optional logical mask, or one mask per subject
%   fillMeasure  'k' (missing-stimulus component; default) or 'k2'
%
% OUTPUTS
%   b       fitted coefficients in the order specified by cols
%   nVox    number of contributing voxels
%   sse     residual sum of squares represented by the cross-products
%   dof     residual degrees of freedom represented by the cross-products
%
% This function intentionally performs no bootstrap or hypothesis test.
% Across-subject standard errors are calculated from subject-level fits in
% fitSampledLinearShifts.

if ~iscell(fits), fits = {fits}; end
if nargin < 3 || isempty(mask)
    mask = cellfun(@(F) true(size(F.ok)),fits,'UniformOutput',false);
elseif ~iscell(mask)
    mask = {mask};
end
if nargin < 4 || isempty(fillMeasure), fillMeasure = 'k'; end
if numel(mask) ~= numel(fits)
    error('prfShiftPool:maskCount','mask must contain one entry per fit.');
end
if isempty(cols) || any(~ismember(cols,1:4)) || numel(unique(cols)) ~= numel(cols)
    error('prfShiftPool:badColumns','cols must contain unique values from 1:4.');
end
switch lower(fillMeasure)
    case {'k','missing'}
        fn = struct('ok','ok','A','WtW','c','WtZ','zz','ZtZ','dof','dof');
    case {'k2','full'}
        fn = struct('ok','ok2','A','WtW2','c','WtZ2','zz','ZtZ2','dof','dof2');
    otherwise
        error('prfShiftPool:badFillMeasure', ...
              'fillMeasure must be ''k'' or ''k2'', not ''%s''.',fillMeasure);
end
nc = numel(cols);
A = zeros(nc); c = zeros(nc,1); zz = 0; dof = 0; nVox = 0;
for s = 1:numel(fits)
    F = fits{s};
    if numel(mask{s}) ~= numel(F.ok)
        error('prfShiftPool:maskSize','Mask %d has the wrong number of voxels.',s);
    end
    required = struct2cell(fn);
    if ~all(cellfun(@(name) isfield(F,name),required))
        error('prfShiftPool:missingFields', ...
              'Fit %d lacks fields required for fillMeasure=''%s''.',s,fillMeasure);
    end
    q = F.(fn.ok) & mask{s}(:);
    if ~any(q), continue, end
    allA = F.(fn.A); allC = F.(fn.c);
    allZZ = F.(fn.zz); allDof = F.(fn.dof);
    A = A+reshape(sum(allA(q,cols,cols),1),nc,nc);
    c = c+reshape(sum(allC(q,cols),1),nc,1);
    zz = zz+sum(allZZ(q));
    dof = dof+sum(allDof(q));
    nVox = nVox+nnz(q);
end
if nVox == 0 || ~all(isfinite(A(:))) || rcond(A) < 1e-12
    b = nan(nc,1); sse = NaN; return
end
b = A\c;
sse = max(zz-b.'*c,0);
end
