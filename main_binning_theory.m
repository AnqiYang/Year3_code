% This script implements the theory that given light level decides optimal
% binning size.

maxI = @(p, f, L) L*p / 2 * abs(sinc(p*f) * p) + L / 2 * p^2;
minI = @(p, f, L) -L*p / 2 * abs(sinc(p*f) * p) + L / 2 * p^2;
mse = @(p, L) 4 + L / 2 * p^2;

% Given different scene light levels, varying from 2.5 to 1000 photoelectron per micron^2
L = 10 .^(linspace(log10(2.5), log10(1000), 10));
% Candidate pixel pitches vary from 0.5 micron to 3.5 micron
p = 3.5:-0.1:0.5;
% Visualize Contrast-to-noise curves at L0 light level
L0 = 50;
% Set the required contrast-to-noise threshold
cnrThreshold = 4;
fs = linspace(0, 2, 2e3);       % frequency cycles per microns
xs = linspace(0, 10e3, 2e3);    % spatial senosr coord in microns
color = colormap('lines');
cutoff = zeros(length(L), length(p));

cutoffFrequencies = {};
close all;

for lightLevel = 1: length(L)
    
    % create a blank figure
    if L(lightLevel) == L0
        figure, set(gcf, 'color', 'k');
        hold on;
        grid on; grid minor;
    end

    % Given lightLevel, compute contrast and noise for each pitch
    for ii = 1: length(p)

        % Plot contrast-to-noise curve
        contrast = maxI(p(ii), fs, L(lightLevel)) - ...
            minI(p(ii), fs, L(lightLevel));
        noiseStd = sqrt(mse(p(ii), L(lightLevel))) * ones(size(fs));

        % Find cutoff frequency
        cnr = contrast ./ noiseStd;
        ids = [];
        epsilon = 1e-9;
        while isempty(ids)
            epsilon = epsilon * 10;
            ids = find((cnr - cnrThreshold) < epsilon);
        end
        ids = min(ids);
        cutoff(lightLevel, ii) = fs(ids);

        % Visualize contrast-to-noise curves at L0 and every 0.5um pitch
        if mod(ii, 5) == 1
            if L(lightLevel) == L0
                plot(fs * 1e3, contrast, ...
                    'LineWidth', 1.5, 'Color', color(ii, :));
                plot(fs * 1e3, noiseStd, ...
                    'LineWidth', 1.5, ...
                    'LineStyle', '--', 'Color', color(ii, :));
                scatter(fs(ids) * 1e3, contrast(ids), 32, 'filled', ...
                    'LineWidth', 2, ...
                    'MarkerEdgeColor', 'w', ...
                    'MarkerFaceColor', color(ii, :));
            end
        end

    end
    % Visualization
    if L(lightLevel) == L0
        xlabel('cycles / mm');
        ylabel('electrons in one pixel');
        legend('3.5um', '', '',...
            '3.0um', '', '',...
            '2.5um', '', '',...
            '2.0um', '', '',...
            '1.5um', '', '', ...
            '1.0um', '', '', ...
            '0.5um', '', '');
        % legend('noise-free signal contrast', 'noise std');
        title(sprintf('Contrast-noise curves (%.1f e-/um^2)', L(lightLevel)));
        set(gca, 'Color', 'black');
        set(gca,'XColor',[1 1 1]);
        set(gca,'YColor',[1 1 1]);
        hold off;
    end
end

figure('Position', [10 10 1200 400]);
set(gcf, 'color', 'w');

hold on;
subplot(121);
grid on; grid minor;
plot(p, cutoff(1:7:end, :) * 1e3, 'LineWidth', 2);
set(gca, 'fontsize', 22);
xlabel('pixel pitch [um]');
ylabel('cutoff frequencey [cycles / mm]');
legend('2.5 e-/um^2', '4 e-/um^2', '10 e-/um^2', ...
    '20 e-/um^2', '50 e-/um^2', '100 e-/um^2', '1000 e-/um^2');

% find the best binning at different light level
pHat = p(1: 5: end);
cutoffHat = cutoff(:, 1: 5: end);
[~, bestBinning] = max(cutoffHat, [], 2);

% additional figure
[~, maxIds] = max(cutoff, [], 2);
bestPitch = zeros(size(maxIds));
for ii = 1: size(maxIds); bestPitch(ii) = p(maxIds(ii)); end;

subplot(122);
semilogx(L, bestPitch, '-', 'LineWidth', 2);
hold on;
set(gca, 'fontsize', 22);
% set(gcf, 'Color', 'black');
% set(gca, 'Color', 'black');
% set(gca,'XColor',[1 1 1]);
% set(gca,'YColor',[1 1 1]);
grid on; grid minor;
xlabel('Light level e-/um^2');
ylabel('Best pixel size');
yticks(0.5: 0.5: max(bestPitch)+0.5);
xticks([2, 4, 10, 100, 1000]);
ylim([1, 9] * 0.5);
hold off;

%% Build a continuous function between L and bestBinning
L = [0, L];
bestPitch = [bestPitch(1); bestPitch];

syms x;
bestBinningFunc(x) = piecewise(x >= 0, 0, x < 0, 0);
for i = 1: length(L) - 1
    L1 = L(i);
    L2 = L(i + 1);
    b1 = (bestPitch(i) / 0.5);
    b2 = (bestPitch(i + 1) / 0.5);

    slope = (b1 - b2) / (L1 - L2);
    intersect = b1 - slope * L1;
    pieceFunc(x) = piecewise(x < L1, 0, ...
        (x >= L1) & (x < L2), slope * x + intersect, ...
        x >= L2, 0);
    bestBinningFunc = bestBinningFunc + pieceFunc;
end
pieceFunc(x) = piecewise(x< L(end), 0, ...
    x >= L(end), 1);
bestBinningFunc = bestBinningFunc + pieceFunc;

save(sprintf('data/bestBinningFunc_by_threshold_%d.mat', cnrThreshold), ...
    'bestBinningFunc');
figure, plot(L, bestBinningFunc(L));