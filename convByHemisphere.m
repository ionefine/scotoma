function Y = convByHemisphere(X,hrf,TR,hemIdx,r)
% convByHemisphere  Convolve each column with its own hemisphere's HRF.
%
%   Y = convByHemisphere(X,hrf,TR,hemIdx,r) where
%     X       [nT x nVox] neural predictions, already projected onto each pRF
%     hrf     cell, hrf{h}{r} is the HRF for hemisphere h on run r
%     hemIdx  [nVox x 1] or scalar; which hemisphere each column belongs to
%     r       run index into hrf{h}
%
% Columns are grouped by hemisphere so each group is convolved once, which
% keeps this to one conv2 call per hemisphere rather than one per voxel.
%
% Callers must project (and, for CSS, apply the exponent) BEFORE calling
% this. See convHRF for why projecting first is exact for the linear model.

nVox = size(X,2);
if isscalar(hemIdx), hemIdx = repmat(hemIdx,1,nVox); end
hemIdx = double(hemIdx(:).');
if numel(hemIdx) ~= nVox
    error('convByHemisphere:badHemIdx', ...
          'hemIdx has %d entries but X has %d columns.',numel(hemIdx),nVox);
end
if any(~isfinite(hemIdx) | hemIdx ~= round(hemIdx) | ...
       hemIdx < 1 | hemIdx > numel(hrf))
    error('convByHemisphere:badHemIdx', ...
          'Every hemIdx entry must be an integer indexing hrf.');
end
if ~iscell(hrf) || isempty(hrf) || ~isscalar(r) || r ~= round(r) || r < 1
    error('convByHemisphere:badHRF','hrf must be a nonempty cell array and r a positive integer.');
end
Y = zeros(size(X));
for h = 1:numel(hrf)
    c = (hemIdx == h);
    if ~any(c), continue, end
    if ~iscell(hrf{h}) || numel(hrf{h}) < r
        error('convByHemisphere:missingRunHRF', ...
              'No HRF was supplied for hemisphere %d, run %d.',h,r);
    end
    Y(:,c) = convHRF(X(:,c),hrf{h}{r},TR);
end
end
