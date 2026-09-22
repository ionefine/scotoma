function hrfParams = loadHRFParams(subNum,tsvFile)
% loadHRFParams  Per-hemisphere HRF parameters for one subject.
%
%   hrfParams = loadHRFParams(subNum) returns a 1-by-2 struct array holding
%   the two-gamma parameters for that subject, ordered {'L','R'} to match
%   hemiOrder elsewhere in this pipeline. Each element carries the fields
%   hrf_twogamma requires: delta, c, a1, a2, b1, b2.
%
% WHY THIS EXISTS. The per-subject prfs MAT files store a single hrfParams
% struct, and that struct is the RIGHT-hemisphere fit for every subject:
% older2/CleanScotomaData.m assigns `hrfParams = tmpprfs.hrfParams` after
% its {L,R} loop has closed, so the last hemisphere loaded always wins.
% Roughly half of every subject's vertices are left-hemisphere and were
% therefore modelled with the wrong HRF. The pRF estimates in those same
% files were fitted per hemisphere with the correct HRF, so the stored
% value also contradicts the estimates it ships with.
%
% The authoritative values are published with the source dataset, so this
% function reads them directly rather than trusting the MAT files:
%   data/participants_hrf_parameters.tsv
% copied verbatim from OpenNeuro ds004698 (CC0),
% derivatives/prf-estimation/files/participants_hrf_parameters.tsv
% (Chang, Fine & Boynton 2025, J Vision 25(1):5).
%
% subNum is the BIDS participant number, i.e. 1..12 for sub-01..sub-12, not
% an index into any subject list.

if nargin < 2 || isempty(tsvFile)
    tsvFile = fullfile(fileparts(mfilename('fullpath')), ...
                       'data','participants_hrf_parameters.tsv');
end
validateattributes(subNum,{'numeric'}, ...
    {'real','finite','scalar','integer','positive'},mfilename,'subNum');
if ~isfile(tsvFile)
    error('loadHRFParams:missingTable', ...
          ['HRF parameter table not found:\n%s\n', ...
           'Copy it from ds004698 derivatives/prf-estimation/files/.'],tsvFile);
end

T = readtable(tsvFile,'FileType','text','Delimiter','\t', ...
              'ReadVariableNames',true);
required = {'participant_id','hemisphere','a1','a2','b1','b2','c','delta'};
if ~all(ismember(required,T.Properties.VariableNames))
    error('loadHRFParams:badTable','%s lacks one or more required columns.',tsvFile);
end

wanted = sprintf('sub-%02d',subNum);
rows = strcmp(string(T.participant_id),wanted);
if ~any(rows)
    error('loadHRFParams:missingSubject','%s has no row for %s.',tsvFile,wanted);
end

hemiOrder = {'L','R'};
hrfParams = repmat(struct('delta',[],'c',[],'a1',[],'a2',[],'b1',[],'b2',[]),1,2);
for h = 1:2
    q = rows & strcmpi(string(T.hemisphere),sprintf('hemi-%s',hemiOrder{h}));
    if nnz(q) ~= 1
        error('loadHRFParams:hemisphereCount', ...
              'Expected exactly one hemi-%s row for %s in %s; found %d.', ...
              hemiOrder{h},wanted,tsvFile,nnz(q));
    end
    hrfParams(h).delta = double(T.delta(q));
    hrfParams(h).c     = double(T.c(q));
    hrfParams(h).a1    = double(T.a1(q));
    hrfParams(h).a2    = double(T.a2(q));
    hrfParams(h).b1    = double(T.b1(q));
    hrfParams(h).b2    = double(T.b2(q));
end
end
