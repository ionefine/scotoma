function Y = convHRF(X,hrf,TR)
% convHRF  Convolve each column of X with hrf along time, scaled by TR.
%
%   Y = convHRF(X,hrf,TR) treats the rows of X as time points and returns a
%   matrix the same size as X, truncated to drop the convolution tail. This
%   is the single place the pipeline turns a neural prediction into a BOLD
%   prediction, so every caller applies the same TR scaling and truncation.
%
% Convolution along time commutes with the spatial projection onto a pRF:
% conv(S*G) equals conv(S)*G, because convolution acts independently down
% each pixel column and G mixes columns. Projecting first and convolving the
% resulting single time course is therefore exact for the linear model, and
% it avoids materialising one convolved [nT x nPix] design per hemisphere.
%
% It does NOT commute with the CSS exponent, so CSS callers must project,
% raise to the power n, and only then convolve.

Y = TR*conv2(X,hrf(:),'full');
Y = Y(1:size(X,1),:);
end
