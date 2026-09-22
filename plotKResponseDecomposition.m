function [fig,D] = plotKResponseDecomposition(T)
% plotKResponseDecomposition  Decompose normalized scotoma response using k.
%
% For each eccentricity bin,
%
%   bottom-up response       = 1-m
%   predictive component    = k*m
%   total scotoma response  = 1-m+k*m
%   reduction from full     = 1-total = m*(1-k)
%
% The full-field response is normalized to 1. The predictive component is
% drawn as a stacked band: its lower edge is the bottom-up response and its
% upper edge is the total response. This decomposition is for k, not k2.
%
% INPUT
%   T  table returned by prfShiftFigure
%
% OUTPUTS
%   fig  figure handle
%   D    table containing the four plotted response quantities

if ~istable(T)
    error('plotKResponseDecomposition:badInput','T must be a table.');
end
required = {'ecc','massIn_median','k','k_lo','k_hi'};
missing = required(~ismember(required,T.Properties.VariableNames));
if ~isempty(missing)
    error('plotKResponseDecomposition:missingVariables', ...
          'T lacks: %s.',strjoin(missing,', '));
end

[ecc,order] = sort(double(T.ecc(:)));
m = double(T.massIn_median(order));
k = double(T.k(order));
kLo = double(T.k_lo(order));
kHi = double(T.k_hi(order));
if any(isfinite(m) & (m < -1e-8 | m > 1+1e-8))
    error('plotKResponseDecomposition:badMass', ...
          'massIn_median must lie between 0 and 1.');
end

bottomUp = 1-m;
predictive = k.*m;
total = bottomUp+predictive;
reduction = 1-total;
totalA = bottomUp+kLo.*m;
totalB = bottomUp+kHi.*m;
totalLo = min(totalA,totalB);
totalHi = max(totalA,totalB);
D = table(ecc,m,k,bottomUp,predictive,total,reduction,totalLo,totalHi, ...
    'VariableNames',{'ecc','m','k','bottomUp','predictive','total', ...
                     'reduction','total_lo','total_hi'});

fig = figure('Color','w','Position',[100 100 720 440]);
ax = axes(fig); hold(ax,'on');
colBottom = [0.90 0.62 0.00];
colPredictive = [0.00 0.45 0.70];
colReduction = [0.82 0.82 0.82];

hReduction = addBand(ax,ecc,total,ones(size(total)),colReduction,0.45);
hBottom = addBand(ax,ecc,zeros(size(bottomUp)),bottomUp,colBottom,0.45);
hPredictive = addBand(ax,ecc,bottomUp,total,colPredictive,0.45);
h1 = yline(ax,1,'k-','LineWidth',1.2);
hBottomLine = plot(ax,ecc,bottomUp,'-','Color',colBottom,'LineWidth',2);
hTotal = plot(ax,ecc,total,'-o','Color',colPredictive, ...
              'MarkerFaceColor',colPredictive,'MarkerSize',5,'LineWidth',2);
plot(ax,ecc,totalLo,':','Color',colPredictive,'LineWidth',1);
plot(ax,ecc,totalHi,':','Color',colPredictive,'LineWidth',1);
xline(ax,2,'--','Color',[0.35 0.35 0.35],'LineWidth',1);

xlabel(ax,'pRF eccentricity (deg)');
ylabel(ax,'Response / full-field response');
title(ax,'Response decomposition implied by k and median pRF mass', ...
      'FontWeight','normal');
validY = [0;1;bottomUp(isfinite(bottomUp));totalLo(isfinite(totalLo)); ...
          totalHi(isfinite(totalHi))];
pad = max(0.05,0.06*(max(validY)-min(validY)));
ylim(ax,[min(validY)-pad,max(validY)+pad]);
if any(isfinite(ecc)), xlim(ax,[min(ecc),max(ecc)]); end
set(ax,'TickDir','out','Box','off','Layer','top');
grid(ax,'on');

h = [hBottom,hPredictive,hReduction,hTotal,h1];
labels = {'bottom-up response: 1-m','predictive addition: k m', ...
          'BOLD reduction: 1-total','total scotoma response', ...
          'full-field response'};
good = isgraphics(h);
legend(ax,h(good),labels(good),'Location','best','Box','off');
end

function h = addBand(ax,x,y1,y2,color,alpha)
% Draw each finite stretch separately so missing bins remain visible.
valid = isfinite(x) & isfinite(y1) & isfinite(y2);
starts = find(valid & [true;~valid(1:end-1)]);
stops = find(valid & [~valid(2:end);true]);
h = gobjects(1);
for j = 1:numel(starts)
    idx = starts(j):stops(j);
    if numel(idx) < 2, continue, end
    p = patch(ax,[x(idx);flipud(x(idx))],[y1(idx);flipud(y2(idx))],color, ...
              'EdgeColor','none','FaceAlpha',alpha);
    if ~isgraphics(h), h = p; end
end
end
