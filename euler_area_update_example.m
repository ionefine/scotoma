figure(1); clf;

% Time step and number of Euler iterations
dt = 0.1;
nIter = 100;

% Model parameters for each visual area
for a = 1:3
    area(a).slope = 0.08;
    area(a).tau_ff = 0.16;
    area(a).tau_td = 0.16;
end
area(1).intercept = 0.16;
area(2).intercept = 0.26;
area(3).intercept = 0.36;

x = linspace(-6, 6, 600);

% Inputs
stim = ones(size(x));
stim(x < 2) = 0;

nAreas = 3;
ff_stim = zeros(nAreas + 1, numel(x));
ff_stim(1, :) = stim;                    % external feedforward input to area 1

td_stim = zeros(nAreas + 1, numel(x));
td_stim(nAreas + 1, :) = 1;              % external top-down input to highest level

% State variables
resp = zeros(nAreas, numel(x));
resp_ff = zeros(nAreas, numel(x));
resp_td = zeros(nAreas, numel(x));

% Euler integration
for n = 1:nIter
    for a = 1:nAreas
        for l = 1:numel(x)
            sigma = area(a).slope * abs(x(l)) + area(a).intercept;
            rf = normpdf(x, x(l), sigma);
            rf_sum = sum(rf);

            if rf_sum > 0
                resp_ff(a, l) = (rf * ff_stim(a, :)') / rf_sum;
                resp_td(a, l) = (rf * td_stim(a, :)') / rf_sum;
            else
                resp_ff(a, l) = 0;
                resp_td(a, l) = 0;
            end
        end

        % Euler update: combine feedforward and top-down drives once
        dresp = (resp_ff(a, :) ./ area(a).tau_ff) + (resp_td(a, :) ./ area(a).tau_td);
        resp(a, :) = resp(a, :) + dt .* dresp;

        % Propagate updated activity to neighboring levels
        if a < nAreas
            ff_stim(a + 1, :) = resp(a, :);
        end
        if a > 1
            td_stim(a - 1, :) = resp(a, :);
        end
    end
end

% Plot responses for each area
for a = 1:nAreas
    subplot(nAreas, 1, a);
    plot(x, resp_ff(a, :), 'b', 'LineWidth', 1.2); hold on;
    plot(x, resp_td(a, :), 'r', 'LineWidth', 1.2);
    plot(x, resp(a, :), 'k--', 'LineWidth', 1.2);
    ylabel(sprintf('Area %d', a));
    if a == 1
        title('Feedforward (blue), Top-down (red), Total response (black dashed)');
    end
    if a == nAreas
        xlabel('x');
    end
    legend({'resp\_ff', 'resp\_td', 'resp'}, 'Location', 'best');
end
