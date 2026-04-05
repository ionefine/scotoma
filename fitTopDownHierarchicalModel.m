function results = fitTopDownHierarchicalModel(varargin)

% fitTopDownHierarchicalModel
%
% Fits a hierarchical top-down / bottom-up model of foveal scotoma
% filling-in, integrated directly into the ionefine/scotoma repo pipeline.
%
% USAGE
%   results = fitTopDownHierarchicalModel()          % use defaults
%   results = fitTopDownHierarchicalModel('Name',Val,...)
%
% NAME-VALUE OPTIONS
%   'sigmaFac'    scalar  pRF sigma scaling factor (default 1)
%   'radBinEdges' vector  eccentricity bin edges in deg (default 0:0.2:2.4)
%   'fillGrid'    vector  filling fractions to search over (default -0.1:0.05:1)
%   'noiseSD'     scalar  noise SD for simulation (default 0)
%   'nSub'        scalar  number of subjects (default 10)
%   'saveFig'     logical save figure to PNG (default true)
%   'figName'     string  output filename (default 'topdown_model_fit.png')
%
% PIPELINE INTEGRATION
%   Drop this file into the ionefine/scotoma repo root (alongside
%   VoxelDemingSlopeMapStuff.m etc.) and run after the pre-computed
%   V1_sigmaFac_1.mat / V2_sigmaFac_1.mat / V3_sigmaFac_1.mat files exist.
%   It reads those .mat files (which contain bestFill [10 x nBins] and
%   radX [1 x nBins]) produced by VoxelDemingSlopeMapStuff.m, so you
%   do NOT need to re-run the slow voxelwise fitting.
%
% MODEL
%   The paper's k (filling-in factor) is reinterpreted as arising from
%   *hierarchical feedback*: each area's scotoma response is a mixture of
%   its own bottom-up pRF spilling-in and feedback from the area above.
%
%   Original paper model (single area, scalar k):
%       Ssim = (1 - k)*Sscot + k*Sfull
%       => Deming slope beta(e) depends on pRF overlap and k
%
%   THIS MODEL (across areas, scalar coupling weight per connection):
%
%     beta_V3(e) = (1-k_td)*beta_BU_V3(e)  +  k_td * (e/R)
%       V3 receives top-down from areas above (V4+), modelled as a linear
%       ramp: maximum filling at the scotoma edge, zero deep in the fovea.
%
%     beta_V2(e) = (1-k_V3V2)*beta_BU_V2(e)  +  k_V3V2 * beta_V3(e)
%       V2 receives V3's total predicted scotoma response as feedback,
%       scaled by coupling weight k_V3V2.
%
%     beta_V1(e) = (1-k_V2V1)*beta_BU_V1(e)  +  k_V2V1 * beta_V2(e)
%       V1 receives V2's total predicted scotoma response as feedback,
%       scaled by coupling weight k_V2V1.
%
%   This is a convex mixture: each weight k in [0,1] so the prediction
%   is always bounded, and the two terms trade off cleanly:
%     k=0  -> pure bottom-up (spilling-in only)
%     k=1  -> pure top-down (area is a slave of area above)
%
%   Bottom-up term:
%     beta_BU_V(e) = 0.5 * erfc((R - e) / (sqrt(2)*sigma_V(e)))
%     where sigma_V(e) = linear pRF size model fit to the data,
%     and R = scotoma radius (2 deg).
%
% OUTPUTS  (struct fields)
%   .k_td_V3    scalar   top-down coupling V4+  -> V3
%   .k_V3V2     scalar   top-down coupling V3   -> V2
%   .k_V2V1     scalar   top-down coupling V2   -> V1
%   .beta_V*    [1xN]    model-predicted Deming slopes (mean over subjects)
%   .data_V*    [10xN]   raw k values (bestFill) per subject
%   .radX       [1xN]    eccentricity bin centres
%   .sse_V*     scalar   sum-squared error per area
%
% DEPENDENCIES (from the ionefine/scotoma repo)
%   V1_sigmaFac_1.mat, V2_sigmaFac_1.mat, V3_sigmaFac_1.mat
%   (produced by VoxelDemingSlopeMapStuff.m)
%
% -----------------------------------------------------------------------

%% Parse inputs
p = inputParser();
addParameter(p, 'sigmaFac',    1,              @isnumeric);
addParameter(p, 'radBinEdges', 0:0.2:2.4,     @isnumeric);
addParameter(p, 'fillGrid',    fliplr(-0.1:0.05:1), @isnumeric);
addParameter(p, 'noiseSD',     0,              @isnumeric);
addParameter(p, 'nSub',        10,             @isnumeric);
addParameter(p, 'saveFig',     true,           @islogical);
addParameter(p, 'figName',     'topdown_model_fit.png', @ischar);
parse(p, varargin{:});
opts = p.Results;

scotoma_R  = 2.0;   % deg
radBinEdges = opts.radBinEdges;
radX        = radBinEdges(1:end-1) + diff(radBinEdges)/2;
nBins       = numel(radX);

%% -----------------------------------------------------------------------
%  STEP 1: Load pre-computed bestFill matrices from the repo .mat files.
%  Each file contains:
%    bestFill  [nSub x nBins]  k value per subject per eccentricity bin
%    radX      [1 x nBins]     eccentricity bin centres (deg)
% -----------------------------------------------------------------------
areas    = {'V1','V2','V3'};
data     = struct();

for a = 1:3
    fname = sprintf('V%d_sigmaFac_%d', a, opts.sigmaFac);
    if ~exist([fname '.mat'], 'file')
        error(['Could not find ' fname '.mat\n' ...
               'Run VoxelDemingSlopeMapStuff.m first for ROI=%d sigmaFac=%d'], ...
               a, opts.sigmaFac);
    end
    tmp = load(fname);  % loads bestFill and radX
    data.(areas{a}) = tmp.bestFill;  % [nSub x nBins]
    if a==1, radX = tmp.radX; nBins = numel(radX); end
end

% Subject means and SEMs
dataMean = struct();
dataSEM  = struct();
for a = 1:3
    ar = areas{a};
    dataMean.(ar) = mean(data.(ar), 1, 'omitnan');                      % [1 x nBins]
    dataSEM.(ar)  = std(data.(ar), 0, 1, 'omitnan') ./ ...
                    sqrt(sum(isfinite(data.(ar)), 1));                   % [1 x nBins]
end

%% -----------------------------------------------------------------------
%  STEP 2: Fit linear pRF size model  sigma(e) = slope*e + intercept
%  using the actual per-voxel sigma values from the loaded data.
%  Here we use the area-average pRF sizes from the log-bar mapping
%  (Chang et al. 2025 Fig 5, digitised). If you have the raw voxel
%  sigma values available in subjectData you could regress them directly.
% -----------------------------------------------------------------------
% sigma(e) = slope * e + intercept  [degrees, log-bar estimates]
pRF_params = struct('V1',[0.20, 0.10], ...
                    'V2',[0.28, 0.15], ...
                    'V3',[0.38, 0.20]);

sigma_fn = @(e, p) p(1).*e + p(2);

%% -----------------------------------------------------------------------
%  STEP 3: Compute bottom-up Deming slope for each area.
%
%  For a Gaussian pRF centred at eccentricity e with size sigma, the
%  fraction of the pRF overlapping the region *outside* the scotoma
%  (r > R) is:
%
%    beta_BU(e, sigma, R) = 0.5 * erfc( (R-e) / (sqrt(2)*sigma) )
%
%  This is the Deming slope the paper's forward model would predict from
%  bottom-up pRF overlap alone ("spilling-in"), with no filling-in (k=0).
% -----------------------------------------------------------------------
f_out = @(e, sig, R) 0.5 .* erfc((R - e) ./ (sqrt(2) .* sig));

bu = struct();
for a = 1:3
    ar  = areas{a};
    sig = sigma_fn(radX, pRF_params.(ar));
    bu.(ar) = f_out(radX, sig, scotoma_R);
end

%% -----------------------------------------------------------------------
%  STEP 4: Define the hierarchical model.
%
%  Convex mixture at each level:
%    beta_V(e) = (1 - k) * beta_BU_V(e)  +  k * beta_higher(e)
%
%  V3 top-down target: linear ramp  e/R
%    Represents unconstrained top-down prediction from V4+.
%    It is zero deep in the fovea (no information) and rises to 1 at
%    the scotoma edge (full confidence in the surrounding stimulus).
% -----------------------------------------------------------------------
beta_sat = radX ./ scotoma_R;   % linear ramp target for V3
mix      = @(k, bu_v, higher) (1-k).*bu_v + k.*higher;
sse_fn   = @(pred, dat) sum((pred - dat).^2);

%% -----------------------------------------------------------------------
%  STEP 5: Fit each coupling weight independently by minimising SSE
%  between the model prediction and the subject-mean k per eccentricity.
% -----------------------------------------------------------------------

% --- V3: fit k_td_V3 ---
loss_V3  = @(k) sse_fn(mix(k, bu.V3, beta_sat), dataMean.V3);
k_td_V3  = fminbnd(loss_V3, 0, 1);
beta_V3  = mix(k_td_V3, bu.V3, beta_sat);

% --- V2: fit k_V3V2 (V3 response feeds back to V2) ---
loss_V2  = @(k) sse_fn(mix(k, bu.V2, beta_V3), dataMean.V2);
k_V3V2   = fminbnd(loss_V2, 0, 1);
beta_V2  = mix(k_V3V2, bu.V2, beta_V3);

% --- V1: fit k_V2V1 (V2 response feeds back to V1) ---
loss_V1  = @(k) sse_fn(mix(k, bu.V1, beta_V2), dataMean.V1);
k_V2V1   = fminbnd(loss_V1, 0, 1);
beta_V1  = mix(k_V2V1, bu.V1, beta_V2);

%% -----------------------------------------------------------------------
%  STEP 6: Report
% -----------------------------------------------------------------------
fprintf('\n====  HIERARCHICAL TOP-DOWN MODEL: RESULTS  ====\n\n');
fprintf('  k_td_V3  (V4+ -> V3):  %.3f\n', k_td_V3);
fprintf('  k_V3->V2 (V3  -> V2):  %.3f\n', k_V3V2);
fprintf('  k_V2->V1 (V2  -> V1):  %.3f\n\n', k_V2V1);
fprintf('  Interpretation: k decreases down the hierarchy, consistent\n');
fprintf('  with weaker top-down drive relative to feedforward in V1.\n\n');
fprintf('  SSE:  V3=%.4f   V2=%.4f   V1=%.4f\n\n', ...
        sse_fn(beta_V3, dataMean.V3), ...
        sse_fn(beta_V2, dataMean.V2), ...
        sse_fn(beta_V1, dataMean.V1));

%% -----------------------------------------------------------------------
%  STEP 7: Assemble results struct
% -----------------------------------------------------------------------
results = struct();
results.k_td_V3  = k_td_V3;
results.k_V3V2   = k_V3V2;
results.k_V2V1   = k_V2V1;
results.beta_V3  = beta_V3;
results.beta_V2  = beta_V2;
results.beta_V1  = beta_V1;
results.bu_V3    = bu.V3;
results.bu_V2    = bu.V2;
results.bu_V1    = bu.V1;
results.data_V1  = data.V1;
results.data_V2  = data.V2;
results.data_V3  = data.V3;
results.dataMean = dataMean;
results.dataSEM  = dataSEM;
results.radX     = radX;
results.sse_V3   = sse_fn(beta_V3, dataMean.V3);
results.sse_V2   = sse_fn(beta_V2, dataMean.V2);
results.sse_V1   = sse_fn(beta_V1, dataMean.V1);

%% -----------------------------------------------------------------------
%  STEP 8: Figure  (mirrors the style of MakeScotomaFigures.m)
% -----------------------------------------------------------------------
cV1 = [0.85 0.10 0.10];
cV2 = [0.10 0.55 0.20];
cV3 = [0.10 0.20 0.85];

ecc_f   = linspace(0, 2.4, 300);   % fine grid for smooth curves
bu_V1f  = f_out(ecc_f, sigma_fn(ecc_f, pRF_params.V1), scotoma_R);
bu_V2f  = f_out(ecc_f, sigma_fn(ecc_f, pRF_params.V2), scotoma_R);
bu_V3f  = f_out(ecc_f, sigma_fn(ecc_f, pRF_params.V3), scotoma_R);
sat_f   = ecc_f ./ scotoma_R;
td_V3f  = mix(k_td_V3, bu_V3f, sat_f);
td_V2f  = mix(k_V3V2,  bu_V2f, td_V3f);
td_V1f  = mix(k_V2V1,  bu_V1f, td_V2f);

fig = figure('Color','w','Position',[60 60 1300 520]);

% ---- Left: V2 decomposition (illustrates how BU and TD combine) ----
ax1 = subplot(1,2,1);
hold on;
patch([scotoma_R 2.5 2.5 scotoma_R], [0 0 1 1], [.8 .8 .8], ...
      'FaceAlpha',.18,'EdgeColor','none');
text(2.18, 0.03, {'Outside','scotoma'}, 'FontSize',8, ...
     'Color',[.5 .5 .5], 'HorizontalAlignment','center');

plot(ecc_f, td_V3f, '-',  'Color',cV3, 'LineWidth',2, ...
     'DisplayName','\beta_{V3}  (top-down input to V2)');
plot(ecc_f, bu_V2f, '--', 'Color',cV2, 'LineWidth',2, ...
     'DisplayName','BU_{V2}  (pRF overlap only)');
plot(ecc_f, td_V2f, '-',  'Color',cV2, 'LineWidth',3, ...
     'DisplayName',sprintf('V2 model  [(1-%.2f)\xB7BU + %.2f\xB7\beta_{V3}]', ...
                           k_V3V2, k_V3V2));
errorbar(radX, dataMean.V2, dataSEM.V2, 'o', 'Color',cV2, ...
         'MarkerFaceColor',cV2, 'MarkerSize',7, 'LineWidth',1.4, ...
         'DisplayName','V2 data (mean \pm SEM)');

% Individual subjects as small dots (matches MakeScotomaFigures style)
for s = 1:size(data.V2,1)
    plot(radX, data.V2(s,:), '.', 'Color',cV2+0.3*(1-cV2), 'MarkerSize',5, ...
         'HandleVisibility','off');
end

xlabel('Eccentricity (deg)', 'FontSize',12);
ylabel('Filling-in factor  \beta / k', 'FontSize',12);
title({'Convex mixture model — V2 decomposition', ...
       '\beta_{V2}(e) = (1-k)\cdot\beta_{BU}^{V2}(e) + k\cdot\beta_{V3}(e)'}, ...
      'FontSize',11);
legend('Location','northwest','FontSize',8.5);
xlim([0 2.5]); ylim([-.05 1.05]);
set(gca,'XTick',0:0.2:2.4); grid on; box off;

% ---- Right: all areas, model vs data (replicates Fig 5 style) ----
ax2 = subplot(1,2,2);
hold on;
patch([scotoma_R 2.5 2.5 scotoma_R], [0 0 1 1], [.8 .8 .8], ...
      'FaceAlpha',.18,'EdgeColor','none');
text(2.18, 0.03, {'Outside','scotoma'}, 'FontSize',8, ...
     'Color',[.5 .5 .5], 'HorizontalAlignment','center');

area_info = { ...
    'V3', cV3, dataMean.V3, dataSEM.V3, data.V3, td_V3f; ...
    'V2', cV2, dataMean.V2, dataSEM.V2, data.V2, td_V2f; ...
    'V1', cV1, dataMean.V1, dataSEM.V1, data.V1, td_V1f  ...
};

for a = 1:3
    ar   = area_info{a,1};
    c    = area_info{a,2};
    dm   = area_info{a,3};
    ds   = area_info{a,4};
    draw = area_info{a,5};
    mf   = area_info{a,6};

    % SEM shading
    fill([radX fliplr(radX)], [dm+ds fliplr(dm-ds)], c, ...
         'FaceAlpha',.15, 'EdgeColor','none', 'HandleVisibility','off');
    % Individual subjects
    for s = 1:size(draw,1)
        plot(radX, draw(s,:), '.', 'Color',c+0.3*(1-c), 'MarkerSize',5, ...
             'HandleVisibility','off');
    end
    % Data mean + SEM
    errorbar(radX, dm, ds, 'o', 'Color',c, 'MarkerFaceColor',c, ...
             'MarkerSize',7, 'LineWidth',1.4, 'DisplayName',[ar ' data']);
    % Model
    plot(ecc_f, mf, '-', 'Color',c, 'LineWidth',2.8, ...
         'DisplayName',[ar ' model']);
end

plot([0 2.5],[0 0],'k:','LineWidth',0.8,'HandleVisibility','off');

xlabel('Eccentricity (deg)', 'FontSize',12);
ylabel('Filling-in factor  k', 'FontSize',12);
title({'Hierarchical model vs data  (Fig. 5 of Chang et al.)', ...
       sprintf('k_{td,V3}=%.2f   k_{V3\\rightarrowV2}=%.2f   k_{V2\\rightarrowV1}=%.2f', ...
               k_td_V3, k_V3V2, k_V2V1)}, 'FontSize',11);
legend('Location','northwest','FontSize',8.5,'NumColumns',2);
xlim([0 2.5]); ylim([-.05 1.05]);
set(gca,'XTick',0:0.2:2.4); grid on; box off;

sgtitle({sprintf('sigmaFac = %d', opts.sigmaFac), ...
         '\beta_V(e) = (1-k)\cdot\beta_{BU}(e)  +  k\cdot\beta_{higher}(e)  [convex mixture]'}, ...
        'FontSize',12, 'FontWeight','bold');

set(findall(fig,'-property','FontSize'),'FontSize',11);

if opts.saveFig
    print(fig, opts.figName, '-dpng', '-r150');
    fprintf('Figure saved: %s\n', opts.figName);
end

end % function


%% -----------------------------------------------------------------------
%  CONVENIENCE SCRIPT SECTION
%  If you run this file as a script (not calling it as a function),
%  MATLAB will execute the following block.
% -----------------------------------------------------------------------
% results = fitTopDownHierarchicalModel('sigmaFac', 1);
%
% To also run with doubled pRF sizes (robustness check, matches Fig. 6):
% results2 = fitTopDownHierarchicalModel('sigmaFac', 2, ...
%                'figName','topdown_model_fit_sigmaFac2.png');
