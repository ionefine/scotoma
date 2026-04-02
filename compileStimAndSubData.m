function [subjectData,stimData,x,y] = compileStimAndSubData(compileOpts)


if ~exist('compileOpts','var')
    compileOpts.ROI = 1;
    compileOpts.radRange = [0,3];
    compileOpts.minvexpl = 0.1;
    compileOpts.minSigma = 0.1;
    compileOpts.subList = 1:10;
    compileOpts.edgeSigma = 0;
end


% Compile data to generate the 'subjectData' and 'stimData' structures:
% subjectData is a struct with fields:
%
%   .Gprf           [nPix x nVox]
%                   pRF weights for this subject
%
%   .prfXY          [nVox x 2]
%                   pRF centers [x y] in degrees
%
%   .w_vox          [nVox x 1]
%                   voxel weights, e.g. pRF fit quality in [0,1]
%
%   .Yfull{r}       [nT_r x nVox]
%                   measured full-stimulus timecourses for run r
%
%   .Yscot{r}       [nT_r x nVox]
%                   measured scotoma-stimulus timecourses for run r
%
% stimData is a struct with fields:
%
%   .Sfull{r}       [nT_r x nPix]
%                   run-specific full-stimulus predictor
%
%   .Sscot{r}       [nT_r x nPix]
%                   run-specific scotoma-stimulus predictor

allSubs = [1,2,3,4,5,6,8,9,10,11];

TR = 1.2;  % hard coded, boo

prfType = 'ses-study1_task-logbar';

for sub = compileOpts.subList
    subNum = allSubs(sub);


    % load prfs from 'study 1'

    prfName = sprintf('data/sub-%02g/sub-%02g_%s_prfs.mat',...
        subNum,subNum,prfType);

    load(prfName);
    allPrfs = prfs;
    clear prfs

    for runNum = 1:3
        fprintf('subject %d, run %d ',subNum,runNum)

        % load the two bold data sets

        boldName = sprintf('data/sub-%02g/sub-%02g_task-%s_run-%g_bold.mat', ...
            subNum,subNum,'scotoma',runNum);
        tmpBold = load(boldName);
        allBold.scotoma = tmpBold.bold;
        boldName = sprintf('data/sub-%02g/sub-%02g_task-%s_run-%g_bold.mat', ...
            subNum,subNum,'logbar',runNum);
        tmpBold = load(boldName);
        allBold.full = tmpBold.bold;

        % load in the full logbar stimulus
        % load 'funcOf' which contains the fields t,x,y

        [stim.full,funcOf] = loadScotomaStimuli(runNum,'logbar',TR);
        stim.scotoma = loadScotomaStimuli(runNum,'scotoma',TR);

        x = funcOf.x;
        y = funcOf.y;
        nx = size(funcOf.x,2);
        ny = size(funcOf.x,1);
        nt = length(funcOf.t);

        %%
        % select a subset of pRFs

        allr = sqrt(allPrfs.x0.^2 + allPrfs.y0.^2);  % pRF distance from fovea


        id = find(allPrfs.varea == compileOpts.ROI  &  allr-compileOpts.edgeSigma*allPrfs.sigma> compileOpts.radRange(1) & allr+ compileOpts.edgeSigma*allPrfs.sigma < compileOpts.radRange(2) &...
            allPrfs.vexpl >compileOpts.minvexpl & allPrfs.sigma>compileOpts.minSigma );

        nVox = length(id);
        fprintf('%d voxels\n',nVox)
        r = allr(id);

        [bold.full,prfs] = subData(allBold.full,allPrfs,id);
        bold.scotoma = subData(allBold.scotoma,allPrfs,id);

        subjectData{sub}.w_vox = allPrfs.vexpl(id);


        clear id % don't need this anymore
        subjectData{sub}.Yfull{runNum} = bold.full.data;
        subjectData{sub}.Yscot{runNum} = bold.scotoma.data;
        subjectData{sub}.prfXY =[prfs.x0(:),prfs.y0(:)];
        subjectData{sub}.sigma = prfs.sigma(:);

        %%
        % Generate the design matrix X

        G = zeros(nx*ny,nVox);


        pRF = [];
        for i=1:nVox
            pRF(i).center = [prfs.x0(i),prfs.y0(i)];
            pRF(i).sig = prfs.sigma(i);
            pRF(i).ar = 1;
            G(:,i) = Gauss(pRF(i),x,y,1);
        end



        % Normalize the pRFs to have equal area (not height).
        % This seems to matter
        G = G./repmat(sum(G),nx*ny,1);


        subjectData{sub}.Gprf = G;


        % Create an HDR based on Boynton et al. '96
        % p.tau = 1.5; %seconds
        % p.n = 3;
        % p.dt = TR;
        % p.delay = 2.25;  %seconds
        % th = 0:TR:30;

        hrf = shiftdim(hrf_twogamma(hrfParams,funcOf.t),-1);

        %hdr = shiftdim(gammapdf(p.n,p.tau,th-p.delay),-1);  %use shiftdim to make it a 1x1xn vector


        % convolve the time-course of each pixel in the stimulus with the HDR.
        convStim.full = TR*convn(stim.full,hrf);
        convStim.full = convStim.full(:,:,1:nt);  %truncate the extra padding after the convolution
        convStim.scotoma = TR*convn(stim.scotoma,hrf);
        convStim.scotoma = convStim.scotoma(:,:,1:nt);  %truncate the extra padding after the convolution



        % Reshape to the matrix S: Columns of S are the predicted time-course of each
        % convolved pixel image.  Rows are the convolved pixel image for each
        % time-point.

        stimData.Sfull{runNum} = reshape(convStim.full,[nx*ny,nt])';
        stimData.Sscot{runNum} = reshape(convStim.scotoma,[nx*ny,nt])';

    end
end
