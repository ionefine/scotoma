function [subjectDataSub,stimDataSub,keepIdx] = subsetSubjectDataByVoxel(subjectData,stimData,compileOpts)
% subsetSubjectDataByVoxel  Select voxels while preserving subject metadata.
%
% The radial filter includes both requested endpoints when edgeSigma is zero.

requiredOpts = {'radRange','minvexpl','minSigma','edgeSigma'};
if nargin < 3 || ~all(isfield(compileOpts,requiredOpts))
    error('subsetSubjectDataByVoxel:missingOptions','compileOpts lacks required fields.');
end
if numel(compileOpts.radRange) ~= 2 || any(isnan(compileOpts.radRange)) || ...
   compileOpts.radRange(2) <= compileOpts.radRange(1) || ...
   ~isscalar(compileOpts.edgeSigma) || ~isfinite(compileOpts.edgeSigma) || ...
   compileOpts.edgeSigma < 0
    error('subsetSubjectDataByVoxel:badOptions','Invalid radRange or edgeSigma.');
end
requiredData = {'Yfull','Yscot','prfXY','sigma','w_vox','Gprf'};
if ~all(isfield(subjectData,requiredData))
    error('subsetSubjectDataByVoxel:missingData','subjectData lacks required fields.');
end
nVox = size(subjectData.prfXY,1);
if size(subjectData.prfXY,2) ~= 2 || numel(subjectData.sigma) ~= nVox
    error('subsetSubjectDataByVoxel:badPRF','prfXY and sigma dimensions disagree.');
end
if numel(subjectData.w_vox) ~= nVox
    error('subsetSubjectDataByVoxel:badWeights','w_vox has the wrong length.');
end
w = subjectData.w_vox(:);

sigma = subjectData.sigma(:);
ecc = hypot(subjectData.prfXY(:,1),subjectData.prfXY(:,2));
lo = ecc-compileOpts.edgeSigma*sigma;
hi = ecc+compileOpts.edgeSigma*sigma;
keepIdx = lo >= compileOpts.radRange(1) & hi <= compileOpts.radRange(2) & ...
          w > compileOpts.minvexpl & sigma >= compileOpts.minSigma & ...
          isfinite(ecc) & isfinite(sigma) & isfinite(w);
fprintf('keeping %d of %d voxels\n',nnz(keepIdx),nVox)

subjectDataSub = subjectData;
for r = 1:numel(subjectData.Yfull)
    if size(subjectData.Yfull{r},2) ~= nVox
        error('subsetSubjectDataByVoxel:badYfull','Yfull{%d} has the wrong voxel count.',r);
    end
    subjectDataSub.Yfull{r} = subjectData.Yfull{r}(:,keepIdx);
end
for r = 1:numel(subjectData.Yscot)
    if size(subjectData.Yscot{r},2) ~= nVox
        error('subsetSubjectDataByVoxel:badYscot','Yscot{%d} has the wrong voxel count.',r);
    end
    subjectDataSub.Yscot{r} = subjectData.Yscot{r}(:,keepIdx);
end
subjectDataSub.prfXY = subjectData.prfXY(keepIdx,:);
subjectDataSub.sigma = subjectData.sigma(keepIdx);
subjectDataSub.w_vox = subjectData.w_vox(keepIdx);
% Hemisphere label and its index into stimData.hrfParams must follow the
% voxels, or a later fit will convolve with the wrong hemisphere's HRF.
if isfield(subjectData,'hemisphere')
    subjectDataSub.hemisphere = subjectData.hemisphere(keepIdx);
end
if isfield(subjectData,'hemIdx')
    subjectDataSub.hemIdx = subjectData.hemIdx(keepIdx);
end
if isfield(subjectData,'sourceVoxelIndex')
    subjectDataSub.sourceVoxelIndex = subjectData.sourceVoxelIndex(keepIdx);
end
if size(subjectData.Gprf,2) ~= nVox
    error('subsetSubjectDataByVoxel:badGprf','Gprf has the wrong voxel count.');
end
subjectDataSub.Gprf = subjectData.Gprf(:,keepIdx);
if isfield(subjectDataSub,'selection')
    subjectDataSub.selection.nEligibleVoxels = nnz(keepIdx);
end
stimDataSub = stimData;
end
