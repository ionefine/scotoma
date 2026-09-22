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
Y = zeros(size(X));
for h = 1:numel(hrf)
    c = (hemIdx == h);
    if ~any(c), continue, end
    Y(:,c) = convHRF(X(:,c),hrf{h}{r},TR);
end
end
