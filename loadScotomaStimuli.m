function [stimImg,funcOf] = loadScotomaStimuli(runNum,stimulus,TR)
% loadScotomaStimuli  Load and resample one scotoma-experiment stimulus.
%
% stimImg is returned as [y x time], sampled at 0:TR:tEnd. The source MAT
% file must contain stimImg [time y x] and funcOf fields t, x, and y.

if ~isscalar(runNum) || ~ismember(runNum,1:3)
    error('loadScotomaStimuli:badRun','runNum must be 1, 2, or 3.');
end
stimulus = char(string(stimulus));
if ~any(strcmp(stimulus,{'logbar','scotoma'}))
    error('loadScotomaStimuli:badStimulus','stimulus must be ''logbar'' or ''scotoma''.');
end
if ~isscalar(TR) || ~isfinite(TR) || TR <= 0
    error('loadScotomaStimuli:badTR','TR must be a positive scalar.');
end

pattern = fullfile('data','stimuli',sprintf('task-%s_run-%d*.mat',stimulus,runNum));
files = dir(pattern);
if numel(files) ~= 1
    error('loadScotomaStimuli:fileCount', ...
          'Expected one file matching %s; found %d.',pattern,numel(files));
end
src = load(fullfile(files.folder,files.name));
if ~isfield(src,'stimImg') || ~isfield(src,'funcOf') || ...
   ~all(isfield(src.funcOf,{'t','x','y'}))
    error('loadScotomaStimuli:badFile','%s lacks stimImg or funcOf.t/x/y.',files.name);
end

raw = double(src.stimImg);
funcOf = src.funcOf;
tRaw = double(funcOf.t(:));
if numel(tRaw) < 2 || any(~isfinite(tRaw)) || any(diff(tRaw) <= 0) || ...
   abs(tRaw(1)) > 1e-9
    error('loadScotomaStimuli:badTime','funcOf.t must start at zero and increase strictly.');
end
if ndims(raw) ~= 3 || size(raw,1) ~= numel(tRaw) || ...
   ~isequal(size(raw,2),size(funcOf.x,1)) || ...
   ~isequal(size(raw,3),size(funcOf.x,2)) || ...
   ~isequal(size(funcOf.x),size(funcOf.y)) || any(~isfinite(raw(:))) || ...
   any(~isfinite(funcOf.x(:))) || any(~isfinite(funcOf.y(:)))
    error('loadScotomaStimuli:badDimensions','stimImg and funcOf dimensions disagree.');
end

t = 0:TR:tRaw(end);
flat = reshape(raw,numel(tRaw),[]);
flat = interp1(tRaw,flat,t,'linear');
if any(~isfinite(flat(:)))
    error('loadScotomaStimuli:interpolation','Stimulus interpolation produced invalid values.');
end
stimImg = permute(reshape(flat,[numel(t),size(raw,2),size(raw,3)]),[2,3,1]);
funcOf.t = t;
end
