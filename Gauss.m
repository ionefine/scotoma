function out = Gauss(pRF,x,y,vectorize)
% Gauss  Evaluate an unnormalised elliptical Gaussian on an x/y grid.
%
% pRF.center = [x0 y0]
% pRF.sig    = standard deviation along the rotated y axis
% pRF.ar     = x-axis/y-axis standard-deviation ratio
% pRF.ang    = rotation in radians (default 0)
% vectorize  return out(:) when true (default false)
%
% The peak is one. Callers that require unit pRF mass must normalise the
% sampled Gaussian after this function returns.

if nargin < 4 || isempty(vectorize), vectorize = false; end
if ~isscalar(pRF) || ~all(isfield(pRF,{'center','sig','ar'}))
    error('Gauss:missingParameters','pRF needs center, sig, and ar.');
end
if ~isfield(pRF,'ang'), pRF.ang = 0; end
validateattributes(pRF.center,{'numeric'},{'real','finite','numel',2},mfilename,'pRF.center');
validateattributes(pRF.sig,{'numeric'},{'real','finite','scalar','positive'},mfilename,'pRF.sig');
validateattributes(pRF.ar,{'numeric'},{'real','finite','scalar','positive'},mfilename,'pRF.ar');
validateattributes(pRF.ang,{'numeric'},{'real','finite','scalar'},mfilename,'pRF.ang');
validateattributes(x,{'numeric'},{'real','finite','nonempty'},mfilename,'x');
validateattributes(y,{'numeric'},{'real','finite','nonempty'},mfilename,'y');
if ~isequal(size(x),size(y))
    error('Gauss:badGrid','x and y must be equal-sized finite grids.');
end
validateattributes(vectorize,{'numeric','logical'},{'scalar'},mfilename,'vectorize');
if ~ismember(double(vectorize),[0,1])
    error('Gauss:badVectorize','vectorize must be true or false.');
end

dx = x-pRF.center(1);
dy = y-pRF.center(2);
rampx = cos(pRF.ang)*dx+sin(pRF.ang)*dy;
rampy = -sin(pRF.ang)*dx+cos(pRF.ang)*dy;
out = exp(-0.5*((rampx/(pRF.ar*pRF.sig)).^2+(rampy/pRF.sig).^2));
if vectorize, out = out(:); end
end
