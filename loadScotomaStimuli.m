function [stimImg,funcOf] = loadScotomaStimuli(runNum,stimulus,TR)
% [stimImg,funcOf] = loadScotomaStimuli(runNum,stimulus,TR)
%
% Loads in stimulus for Kelly's scotoma experiment, shifts the dimensions
% and interpolates the stimulus to match the functional data using 'TR'
%
% Assumes files (e.g. task-logbar_run-1_seed-80_rec-down_stim.mat) are in the directory 'data/stimuli/'
%
% runNum: [1,2,3]  (three runs per session/subject/stimulus)  each run used
%                   a different stimulus
% stimulus: ['logbar','scotoma']
%
% load the stimulus

stimName = sprintf('data/stimuli/task-%s_run-%g*.mat',stimulus,runNum);

dr = dir(stimName);
stim = load(['data/stimuli/',dr.name]);

funcOf = stim.funcOf;

% resample stimulus to match TRs

t = 0:TR:max(stim.funcOf.t);
stimImg = zeros(size(stim.funcOf.x,1),size(stim.funcOf.x,2),length(t));

for i=1:size(stim.funcOf.x,1)
    for j=1:size(stim.funcOf.x,2)
        stimImg(i,j,:)= interp1(stim.funcOf.t,stim.stimImg(:,i,j),t);
    end
end

funcOf.t = t;


