function [subjectDataSub, stimDataSub] = subsetSubjectDataByVoxel(subjectData, stimData, compileOpts)
% subsetSubjectDataByVoxel
%
% Subset a subjectData struct to only include selected voxels.
%
% INPUTS
%   subjectData   struct for one subject
%   stimData      stimulus struct (passed through unchanged)
%   compileOpts
%
% OUTPUTS
%   subjectDataSub   same structure but with voxel dimension reduced
%   stimDataSub      identical to input stimData
%
% The function subsets all fields that depend on voxels:
%   - Yfull{r}   [nT x nVox]
%   - Yscot{r}
%   - prfXY      [nVox x 2]
%   - sigma      [nVox x 1]
%   - w_vox      [nVox x 1] (optional)
%   - Gprf       [nPix x nVox] (optional)

allr = sqrt(subjectData.prfXY(:,1).^2+subjectData.prfXY(:,2).^2);  % pRF distance from fovea

keepIdx  = allr-compileOpts.edgeSigma*subjectData.sigma> compileOpts.radRange(1)  & ...
    allr+compileOpts.edgeSigma*subjectData.sigma < compileOpts.radRange(2) &...
    subjectData.w_vox(:) >compileOpts.minvexpl & subjectData.sigma(:)>compileOpts.minSigma;

nKeep = sum(keepIdx);

if nKeep == 0
    warning('No voxels selected (keepIdx all false).');
else
    fprintf('keeping %d of %d voxels\n',nKeep,length(allr));
end

subjectDataSub = subjectData;

% ---------------- Yfull / Yscot ----------------
if isfield(subjectData, 'Yfull')
    for r = 1:numel(subjectData.Yfull)
        Y = subjectData.Yfull{r};
        subjectDataSub.Yfull{r} = Y(:, keepIdx);
    end
end

if isfield(subjectData, 'Yscot')
    for r = 1:numel(subjectData.Yscot)
        Y = subjectData.Yscot{r};
        subjectDataSub.Yscot{r} = Y(:, keepIdx);
    end
end

% ---------------- pRF parameters ----------------
if isfield(subjectData, 'prfXY')
    subjectDataSub.prfXY = subjectData.prfXY(keepIdx, :);
end

if isfield(subjectData, 'sigma')
    subjectDataSub.sigma = subjectData.sigma(keepIdx);
end

if isfield(subjectData, 'w_vox')
    subjectDataSub.w_vox = subjectData.w_vox(keepIdx);
end

% ---------------- Gprf (if present) ----------------
if isfield(subjectData, 'Gprf')
    % Gprf is [nPix x nVox]
    subjectDataSub.Gprf = subjectData.Gprf(:, keepIdx);
end

% ---------------- stimData unchanged ----------------
stimDataSub = stimData;

end