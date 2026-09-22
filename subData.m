function [bold,prfs] = subData(fullBold,fullPrfs,id)
% [bold,prfs] = subData(fullBold,fullPrfs,id)  Subset matched bold and pRF data.
%
% The BOLD and pRF arrays must describe the same vertices in the same order.
% This check prevents a valid-looking index from attaching one vertex's BOLD
% time course to another vertex's pRF parameters.
requiredBold = {'vertex','varea','data','hemisphere'};
requiredPrf = {'vertex','varea','x0','y0','sigma','vexpl','hemisphere'};
if ~all(isfield(fullBold,requiredBold)) || ~all(isfield(fullPrfs,requiredPrf))
    error('subData:missingFields','BOLD or pRF structure lacks required fields.');
end
nBold = size(fullBold.data,2);
nPrf = numel(fullPrfs.vertex);
if numel(fullBold.vertex) ~= nBold || numel(fullBold.varea) ~= nBold || ...
        numel(fullBold.hemisphere) ~= nBold || nBold ~= nPrf || ...
        numel(fullPrfs.varea) ~= nPrf || numel(fullPrfs.hemisphere) ~= nPrf
    error('subData:unmatchedSizes','BOLD and pRF vertex arrays have different lengths.');
end
boldVertex = double(fullBold.vertex(:));
prfVertex = double(fullPrfs.vertex(:));
boldHem = upper(strtrim(string(fullBold.hemisphere(:))));
prfHem = upper(strtrim(string(fullPrfs.hemisphere(:))));
if any(boldVertex ~= prfVertex) || any(boldHem ~= prfHem) || ...
        any(double(fullBold.varea(:)) ~= double(fullPrfs.varea(:)))
    error('subData:vertexOrderMismatch', ...
        ['BOLD and pRF arrays are not in identical vertex/hemisphere order. ' ...
         'Align them by vertex identity before calling subData.']);
end
id = id(:);
if any(~isfinite(id) | id ~= round(id) | id < 1 | id > nPrf) || ...
        numel(unique(id)) ~= numel(id)
    error('subData:badIndex','id must contain unique valid vertex indices.');
end
bold.vertex = fullBold.vertex(id);
bold.varea = fullBold.varea(id);
bold.data = fullBold.data(:,id);
bold.hemisphere = fullBold.hemisphere(id);


prfs.vertex = fullPrfs.vertex(id);
prfs.varea = fullPrfs.varea(id);
prfs.x0 = fullPrfs.x0(id);
prfs.y0 = fullPrfs.y0(id);
prfs.sigma = fullPrfs.sigma(id);
prfs.vexpl = fullPrfs.vexpl(id);
prfs.hemisphere = fullPrfs.hemisphere(id);
