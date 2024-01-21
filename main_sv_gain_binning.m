% This script simulate noisy images from an HDR scene with these modes:
% 1, ground-truth noise free
% 2, no binning + const gain
% 3, no binning + spatially-varying (ROI) gain
% 4, no binning + per-pixel gain
% 5, additive binning + const gain
% 6, additive binning + spatially-varying gain
% 7, average binning + const gain
% 8, average binning + spatially-varying gain
% 9, ISP (digital) binning + const gain
% 10, ISP (digital) binning + spatially-varying gain


close all; clear; clc;
addpath utils;
addpath capture_modes;

no_binning_const_gain_snr = [];
additive_binning_const_gain_snr = [];
average_binning_const_gain_snr = [];
ISP_binning_const_gain_snr = [];

no_binning_sv_gain_snr = [];
additive_binning_sv_gain_snr = [];
average_binning_sv_gain_snr = [];
ISP_binning_sv_gain_snr = [];

no_binning_sv_gain_gain = [];
additive_binning_sv_gain_gain = [];
average_binning_sv_gain_gain = [];
ISP_binning_sv_gain_gain = [];

no_binning_sv_gain_binning = [];
additive_binning_sv_gain_binning = [];
average_binning_sv_gain_binning = [];
ISP_binning_sv_gain_binning = [];


% Define scene
imName = 'WallDrug';
scene = get_scene([imName, '.exr']);

% Define sensor
sensor.maxVolt = 1000;                          % electrons / voltage
sensor.ADC_std = 3.31;                          % post-amplifier readnoise
sensor.read_std = 0.33;                         % pre-amplifier readnoise
sensor.bit = 12;                                % digital image bit
sensor.blackLevel = round(0.02 * (2 ^ sensor.bit - 1)); % electrons / voltage
sensor.dx = 0.5;                                % pixel pitch
% binning and gain are two variables we want to control.
sensor.saturateVolt = (1 - 0.02) * sensor.maxVolt;
% sensor saturation occurs when #electrons are collected and amplified

sensor.gamma = 1 / 3.2;
sensor.offset = 0.000;
sensor.c = 10.;                                  % ISP tonemapping

load('data/bestBinningFunc_by_threshold_6.mat');
%% 1, noise free image (no noise, no quantization, no saturation)
[digital_gt, nl_gt] = capture_gt(scene, sensor, [], 1);
maxVal = max(digital_gt(:));

%% 2, capture with the smallest pixel pitch, 0.5um, constant gain
% Decide gain based on the max scene light level 'scene.L'
% g*q ~ N(g*L, g^2*L)
% Let u + 2mu be saturated value, so that theoretically only 2.2% pixel
% will be saturated.
% Thus we let g*L + 2*g*sqrt(L) = L_satuarte
% g = L_saturate / (L + 2 * sqrt(L))
% Note that g \in [1, 500]

gain = sensor.saturateVolt / (scene.L + 2 * sqrt(scene.L));
gain = min(500, max(1, gain));

[digital_1by1, nl_digital_1by1] = capture_with_no_binning( ...
    scene, sensor, 1, gain);

mse = norm(digital_1by1(:) - digital_gt(:), 2) / sqrt(length(digital_gt(:)));
psnr = 20 * log10(maxVal / mse);
[ssimVal, digital_1by1_ssimMap] = ssim(digital_1by1 / (2 ^ sensor.bit - 1), ...
                          digital_gt / (2 ^ sensor.bit - 1));

% Compute SNR wrt. varying patch brightness
patchSize = 128;
[patchDRs, patchSNRs] = compute_snr_by_patch(digital_1by1, digital_gt, ...
    patchSize);

% Visualize SNR curve
figure(3), hold on; 
set(gcf, 'color', 'w');
set(gca, 'FontSize', 24);
xlabel('Patch brightness log2(L/L_{max})');
colororder('default');
ylabel('SNR');
grid on; grid minor;
plot(patchDRs, patchSNRs, 'LineWidth', 2); 
SNR_legend = {'No binning + fix gain'}; hold off;

% Visualize captured image, binning, and gain map
figure(1), hold on;
subplot(4,3,1);
imshow(nl_digital_1by1);
title(sprintf('No binning + fix gain, psnr = %.2fdB, ssim = %.4f', psnr, ssimVal));

[h, w, ~] = size(scene.latentIm);
subf = subplot(4,3,2);
imagesc(ones(h, w)); axis off;
caxis(subf, [1, 10]);
title('No binning');

subf = subplot(4,3,3);
imagesc(ones(h, w) * gain); axis off;
title(sprintf('Fix gain = %.2f', gain));
caxis(subf, [1, 500]);
hold off;

%% 3, capture with the smallest pixel pitch (no binning), varying gain for each ROI
patchSize = 128;
for row = 1: patchSize: size(scene.latentIm, 1)
    for col = 1: patchSize: size(scene.latentIm, 2)

        latentImPatch = scene.latentIm(row: min(h, row + patchSize - 1), ...
            col: min(w, col + patchSize - 1), :);
        gtPatch = digital_gt(row: min(h, row + patchSize - 1), ...
            col: min(w, col + patchSize - 1), :);

        scenePatch.latentIm = latentImPatch;
        scenePatch.L = scene.L;

        % Decide gain based on max light level within each ROI
        maxL = scene.L * max(latentImPatch(:));
        gain = sensor.saturateVolt / (maxL + 2 * sqrt(maxL));
        gain = min(500, max(1, gain));
        [patch, nl_patch, saturated] = capture_with_no_binning( ...
            scenePatch, sensor, 1, gain);

        % Save binning and gain for visualization
        digital_1by1_varying(row: min(h, row + patchSize - 1), ...
            col: min(w, col + patchSize - 1), :) = patch;
        nl_digital_1by1_varying(row: min(h, row + patchSize - 1), ...
            col: min(w, col + patchSize - 1), :) = nl_patch;
        original_varying_gain(row: min(h, row + patchSize - 1), ...
            col: min(w, col + patchSize - 1), :) = gain;
        original_varying_binning(row: min(h, row + patchSize - 1), ...
            col: min(w, col + patchSize - 1), :) = 1;
        original_varying_saturated(row: min(h, row + patchSize - 1), ...
            col: min(w, col + patchSize - 1), :) = saturated;
    end
end
% Compute psnr of captured image
mse = norm(digital_1by1_varying(:) - digital_gt(:), 2) / sqrt(length(digital_gt(:)));
psnr = 20 * log10(maxVal / mse);
[ssimVal, digital_1by1_varying_ssimMap] = ssim(digital_1by1_varying / (2 ^ sensor.bit - 1), ...
                          digital_gt / (2 ^ sensor.bit - 1));

% Compute and visualize SNRs for varying patch brightness
[patchDRs, patchSNRs] = compute_snr_by_patch(digital_1by1_varying, ...
    digital_gt, patchSize);
figure(3), hold on; plot(patchDRs, patchSNRs, '*-', 'LineWidth', 2); 
SNR_legend = [SNR_legend, 'No binning + SV gain']; hold off;

% Visualize captured image, binning,and gain map
figure(2), hold on;
subplot(431);
imshow(nl_digital_1by1_varying);
title(sprintf('No binning + SV gain, psnr = %.2fdB, ssim = %.4f', psnr, ssimVal));

subf = subplot(432);
imagesc(original_varying_binning);
caxis(subf, [1, 10]);
title('No binning'); axis off;

subf = subplot(433);
imagesc(original_varying_gain);
caxis(subf, [1, 500]);
title('SV gain'); axis off;
hold off;

%% 4, capture with no binning, per-pixel adaptive gain
mode = 3;  % ignore
[h, w, ~] = size(scene.latentIm);
digital_1by1_pixelwise = zeros(h, w, 3);
nl_digital_1by1_pixelwise = zeros(h, w, 3);
original_pixelwise_gain = zeros(h, w, 3);
original_pixelwise_saturated = zeros(h, w, 3);

gain = 1;
for i = 1: size(scene.latentIm, 1)

    % sensor measure voltage, amplify, ADC and readout
    sceneLine.latentIm = scene.latentIm(i, :, :);
    sceneLine.L = scene.L;

    [digital, nl_digital, saturated] = capture_with_no_binning( ...
        sceneLine, sensor, 1, gain);
    digital_1by1_pixelwise(i, :, :) = digital;  % (only digital is observable)
    nl_digital_1by1_pixelwise(i, :, :) = nl_digital;
    original_pixelwise_gain(i, :, :) = gain;
    original_pixelwise_saturated(i, :, :) = saturated;

    % predict gain that scale voltage to maximum voltage
    % clip digital to avoid numerical issue
    gain = predict_gain(digital, saturated, sensor, mode);
end
% Compute PSNR of captured image
mse = norm(digital_1by1_pixelwise(:) - digital_gt(:), 2) / ...
    sqrt(length(digital_gt(:)));
psnr = 20 * log10(maxVal / mse);
[ssimVal, digital_1by1_pixelwise_ssimMap] = ssim(digital_1by1_pixelwise / (2 ^ sensor.bit - 1), ...
                          digital_gt / (2 ^ sensor.bit - 1));

% Compute and visualize SNRs for varying patch brightness
[patchDRs, patchSNRs] = compute_snr_by_patch(digital_1by1_pixelwise, ...
    digital_gt, patchSize);
figure(3), hold on; plot(patchDRs, patchSNRs, 'LineWidth', 2); 
SNR_legend = [SNR_legend, 'no binning + pixelwise gain ', num2str(mode)]; hold off;

% Visualize captured image, binning, and gain map
figure(4), hold on;
subplot(431);
imshow(nl_digital_1by1_pixelwise);
title(sprintf('No binning + pixelwise gain, psnr = %.2fdB, ssim = %.4f', psnr, ssimVal));

subf = subplot(432);
imagesc(ones(h, w));
caxis(subf, [1, 10]);
title('No binning'); axis off;

subf = subplot(433);
imagesc(original_pixelwise_gain(:,:,2));
caxis(subf, [1, 500]);
title('pixelwise gain'); axis off;
hold off;
%% 5, capture with spatially-varying additive binning, constant gain
patchSize = 128;
% Decide gain for the entire image
gain = sensor.saturateVolt / (scene.L + 2 * sqrt(scene.L));  % Fix gain!
gain = min(500, max(1, gain));

for row = 1: patchSize: size(scene.latentIm, 1)
    for col = 1: patchSize: size(scene.latentIm, 2)

        x0 = row;
        x1 = min(h, row + patchSize - 1);
        y0 = col;
        y1 = min(w, col + patchSize - 1);

        latentImPatch = scene.latentIm(x0: x1, y0: y1, :);
        gtPatch = digital_gt(x0: x1, y0: y1, :);

        scenePatch.latentIm = latentImPatch;
        scenePatch.L = scene.L;

        % Select optimal binning size based on theory
        maxL = scene.L * max(latentImPatch(:));
        max_binning = floor(double( ...
            sqrt(sensor.saturateVolt / maxL / gain)));
        binning = round(double(bestBinningFunc(maxL)));
        binning = min(max_binning, max(1, binning));

        [patch, nl_patch] = capture_with_additive_binning( ...
            scenePatch, sensor, binning, gain);

        % Save binning and gain for visualization
        digital_additive(x0: x1, y0: y1, :) = patch;
        nl_digital_additive(x0: x1, y0: y1, :) = nl_patch;
        additive_gain(x0: x1, y0: y1) = gain;
        additive_binning(x0: x1, y0: y1) = binning;
    end
end
% Compute psnr of the image
mse = norm(digital_additive(:) - digital_gt(:), 2) / ...
    sqrt(length(digital_gt(:)));
psnr = 20 * log10(maxVal / mse);
[ssimVal, digital_additive_ssimMap] = ssim(digital_additive / (2 ^ sensor.bit - 1), ...
                          digital_gt / (2 ^ sensor.bit - 1));

% Compute and visualize SNRs for varying patch brightness
[patchDRs, patchSNRs] = compute_snr_by_patch(digital_additive, ...
    digital_gt, patchSize);
figure(3), hold on; plot(patchDRs, patchSNRs, 'LineWidth', 2); 
SNR_legend = [SNR_legend, 'Addtive binning + fix gain']; hold off;

% Visualize captured image, binning, and gain map
figure(1), hold on;
subplot(434);
imshow(nl_digital_additive);
title(sprintf('Addtive binning + fix gain, psnr = %.2fdB, ssim = %.4f', psnr, ssimVal));

subf = subplot(435);
imagesc(additive_binning);
caxis(subf, [1, 10]);
title('Addtive binning'); axis off;

subf = subplot(436);
imagesc(additive_gain); axis off;
caxis(subf, [1, 500]);
title(sprintf('Fix gain = %.2f', gain));
hold off;

%% 6, capture with spatially-varying additive binning, spatially varying gain
patchSize = 128;

for row = 1: patchSize: size(scene.latentIm, 1)
    for col = 1: patchSize: size(scene.latentIm, 2)

        latentImPatch = scene.latentIm(row: min(h, row + patchSize - 1), ...
            col: min(w, col + patchSize - 1), :);
        gtPatch = digital_gt(row: min(h, row + patchSize - 1), ...
            col: min(w, col + patchSize - 1), :);

        scenePatch.latentIm = latentImPatch;
        scenePatch.L = scene.L;

        % Select optimal binning size
        maxL = scene.L * max(latentImPatch(:));
        max_binning = floor(double(sqrt( ...
            sensor.saturateVolt / maxL)));
        binning = round(double(bestBinningFunc(maxL)));
        binning = min(max_binning, max(1, binning));

        % Decide gain
        maxL_bin = maxL * (binning ^ 2);
        gain = sensor.saturateVolt / (maxL_bin + 2 * sqrt(maxL_bin));
        gain = min(500, max(1, gain));
        [patch, nl_patch] = capture_with_additive_binning( ...
            scenePatch, sensor, binning, gain);
        
        % Save binning and gain for visualization
        digital_additive_varying(row: min(h, row + patchSize - 1), ...
            col: min(w, col + patchSize - 1), :) = patch;
        nl_digital_additive_varying(row: min(h, row + patchSize - 1), ...
            col: min(w, col + patchSize - 1), :) = nl_patch;
        additive_varying_gain(row: min(h, row + patchSize - 1), ...
            col: min(w, col + patchSize - 1)) = gain;
        additive_varying_binning(row: min(h, row + patchSize - 1), ...
            col: min(w, col + patchSize - 1)) = binning;
    end
end
% compute pnsr of the captured image
mse = norm(digital_additive_varying(:) - digital_gt(:), 2) / ...
    sqrt(length(digital_gt(:)));
psnr = 20 * log10(maxVal / mse);
[ssimVal, digital_additive_varying_ssimMap] = ssim(digital_additive_varying / (2 ^ sensor.bit - 1), ...
                          digital_gt / (2 ^ sensor.bit - 1));

% compute and visualize SNRs for varying patch brightness
[patchDRs, patchSNRs] = compute_snr_by_patch(digital_additive_varying, ...
    digital_gt, patchSize);
figure(3), hold on; plot(patchDRs, patchSNRs, '>-', 'LineWidth', 3, 'MarkerSize', 10);
SNR_legend = [SNR_legend, 'Addtive binning + SV gain']; hold off;

% visualize captured, binning, and gain map
figure(2), hold on;
subplot(434);
imshow(nl_digital_additive_varying)
title(sprintf('Addtive binning + SV gain, psnr = %.2fdB, ssim = %.4f', psnr, ssimVal));

subf = subplot(435);
imagesc(additive_varying_binning);
caxis(subf, [1, 10]);
title('SV binning'); axis off;

subf = subplot(436);
imagesc(additive_varying_gain);
caxis(subf, [1, 500]);
title('SV gain'); axis off;
hold off;

%% 7, capture with spatially-varying average binning, constant gain
% Decide gain for the entire image
patchSize = 128;
gain = sensor.saturateVolt / (scene.L + 2 * sqrt(scene.L));
gain = min(500, max(1, gain));

for row = 1: patchSize: size(scene.latentIm, 1)
    for col = 1: patchSize: size(scene.latentIm, 2)

        latentImPatch = scene.latentIm(row: min(h, row + patchSize - 1), ...
            col: min(w, col + patchSize - 1), :);
        gtPatch = digital_gt(row: min(h, row + patchSize - 1), ...
            col: min(w, col + patchSize - 1), :);

        scenePatch.latentIm = latentImPatch;
        scenePatch.L = scene.L;

        % Select optimal binning size
        maxL = scene.L * max(latentImPatch(:));
        binning = round(double(bestBinningFunc([maxL])));
        binning = max(1, binning);
        [patch, nl_patch] = capture_with_average_binning( ...
            scenePatch, sensor, binning, gain);

        % Save binning and gain for visualization
        digital_average(row: min(h, row + patchSize - 1), ...
            col: min(w, col + patchSize - 1), :) = patch;
        nl_digital_average(row: min(h, row + patchSize - 1), ...
            col: min(w, col + patchSize - 1), :) = nl_patch;
        average_gain(row: min(h, row + patchSize - 1), ...
            col: min(w, col + patchSize - 1)) = gain;
        average_binning(row: min(h, row + patchSize - 1), ...
            col: min(w, col + patchSize - 1)) = binning;
    end
end
% Compute PNSR of the captured image
mse = norm(digital_average(:) - digital_gt(:), 2) / ...
    sqrt(length(digital_gt(:)));
psnr = 20 * log10(maxVal / mse);
[ssimVal, digital_average_ssimMap] = ssim(digital_average / (2 ^ sensor.bit - 1), ...
                          digital_gt / (2 ^ sensor.bit - 1));

% Compute and visualize SNRs for varying patch brightness
[patchDRs, patchSNRs] = compute_snr_by_patch(digital_average, ...
    digital_gt, patchSize);
figure(3), hold on; plot(patchDRs, patchSNRs, 'LineWidth', 2); 
SNR_legend = [SNR_legend, 'Average binning + fix gain']; hold off;

% Visualize captured image, binning, and gain map
figure(1), hold on;
subplot(437);
imshow(nl_digital_average)
title(sprintf('Average binning + fix gain, psnr = %.2fdB, ssim = %.4f', psnr, ssimVal));

subf = subplot(438);
imagesc(average_binning); axis off;
caxis(subf, [1, 10]);
title('SV binning');

subf = subplot(439);
imagesc(average_gain); axis off;
caxis(subf, [1, 500]);
title(sprintf('Fix gain = %.2f', gain));
hold off;

%% 8, capture with spatially-varying average binning, varying gain
patchSize = 128;

for row = 1: patchSize: size(scene.latentIm, 1)
    for col = 1: patchSize: size(scene.latentIm, 2)

        latentImPatch = scene.latentIm(row: min(h, row + patchSize - 1), ...
            col: min(w, col + patchSize - 1), :);
        gtPatch = digital_gt(row: min(h, row + patchSize - 1), ...
            col: min(w, col + patchSize - 1), :);

        scenePatch.latentIm = latentImPatch;
        scenePatch.L = scene.L;

        % Decide optimal binning size
        maxL = scene.L * max(latentImPatch(:));
        binning = round(double(bestBinningFunc([maxL])));
        binning = max(1, binning);

        % Decide gain
        gain = sensor.saturateVolt / (maxL + 2 * sqrt(maxL / binning ^ 2));
        gain = min(500, max(1, gain));
        [patch, nl_patch] = capture_with_average_binning( ...
            scenePatch, sensor, binning, gain);

        % Save binning and gain for visualizaiton
        digital_average_varying(row: min(h, row + patchSize - 1), ...
            col: min(w, col + patchSize - 1), :) = patch;
        nl_digital_average_varying(row: min(h, row + patchSize - 1), ...
            col: min(w, col + patchSize - 1), :) = nl_patch;
        average_varying_gain(row: min(h, row + patchSize - 1), ...
            col: min(w, col + patchSize - 1)) = gain;
        average_varying_binning(row: min(h, row + patchSize - 1), ...
            col: min(w, col + patchSize - 1)) = binning;
    end
end
% Compute psnr of the captured image
mse = norm(digital_average_varying(:) - digital_gt(:), 2) / ...
    sqrt(length(digital_gt(:)));
psnr = 20 * log10(maxVal / mse);
[ssimVal, digital_average_varying_ssimMap] = ssim(digital_average_varying / (2 ^ sensor.bit - 1), ...
                          digital_gt / (2 ^ sensor.bit - 1));

% Compute and visualize SNRs of varying patch brightness
[patchDRs, patchSNRs] = compute_snr_by_patch(digital_average_varying, ...
    digital_gt, patchSize);
figure(3), hold on; plot(patchDRs, patchSNRs, 'o-', 'LineWidth', 2); 
SNR_legend = [SNR_legend, 'Average binning + SV gain']; hold off;

% Visualize captured image, binning, and gain map
figure(2), hold on;
subplot(437);
imshow(nl_digital_average_varying)
title(sprintf('Average binning + SV gain, psnr = %.2fdB, ssim = %.4f', psnr, ssimVal));

subf = subplot(438);
imagesc(average_varying_binning);
caxis(subf, [1, 10]);
title('SV binning'); axis off;

subf = subplot(439);
imagesc(average_varying_gain);
caxis(subf, [1, 500]);
title('SV gain'); axis off;
hold off;

%% 9, capture with spatially-varying digital binning, constant gain
% Decide gain for the entire image
patchSize = 128;
gain = sensor.saturateVolt / (scene.L + 2 * sqrt(scene.L));
gain = min(500, max(1, gain));

for row = 1: patchSize: size(scene.latentIm, 1)
    for col = 1: patchSize: size(scene.latentIm, 2)

        latentImPatch = scene.latentIm(row: min(h, row + patchSize - 1), ...
            col: min(w, col + patchSize - 1), :);
        gtPatch = digital_gt(row: min(h, row + patchSize - 1), ...
            col: min(w, col + patchSize - 1), :);

        scenePatch.latentIm = latentImPatch;
        scenePatch.L = scene.L;

        % Select optimal binning size
        maxL = scene.L * max(latentImPatch(:));
        binning = round(double(bestBinningFunc(maxL)));
        binning = max(1, binning);
        [patch, nl_patch] = capture_with_ISP_binning( ...
            scenePatch, sensor, binning, gain);

        % Save binning and gain for visualization
        digital_ISP(row: min(h, row + patchSize - 1), ...
            col: min(w, col + patchSize - 1), :) = patch;
        nl_digital_ISP(row: min(h, row + patchSize - 1), ...
            col: min(w, col + patchSize - 1), :) = nl_patch;
        ISP_gain(row: min(h, row + patchSize - 1), ...
            col: min(w, col + patchSize - 1)) = gain;
        ISP_binning(row: min(h, row + patchSize - 1), ...
            col: min(w, col + patchSize - 1)) = binning;
    end
end
% Compute pnsr of the captured image
mse = norm(digital_ISP(:) - digital_gt(:), 2) / ...
    sqrt(length(digital_gt(:)));
psnr = 20 * log10(maxVal / mse);
[ssimVal, digital_ISP_ssimMap] = ssim(digital_ISP / (2 ^ sensor.bit - 1), ...
                          digital_gt / (2 ^ sensor.bit - 1));

% compute and visualize SNRs for varying patch brightness
[patchDRs, patchSNRs] = compute_snr_by_patch(digital_ISP, ...
    digital_gt, patchSize);
figure(3), hold on; plot(patchDRs, patchSNRs, 'LineWidth', 2); 
SNR_legend = [SNR_legend, 'ISP binning + fix gain']; hold off;

% visualize captured image, binning, and gain map
figure(1), hold on;
subplot(4,3,10);
imshow(nl_digital_ISP);
title(sprintf('ISP processing + fix gain, psnr = %.2fdB, ssim = %.4f', psnr, ssimVal));

subf = subplot(4,3,11);
imagesc(ISP_binning); axis off;
caxis(subf, [1, 10]);
title('ISP filter size');

subf = subplot(4,3,12);
imagesc(ISP_gain); axis off;
caxis(subf, [1, 500]);
title(sprintf('Fix gain = %.2f', gain));
hold off;

%% 10, capture with spatially-varying digital binning, varying gain
patchSize = 128;

for row = 1: patchSize: size(scene.latentIm, 1)
    for col = 1: patchSize: size(scene.latentIm, 2)

        latentImPatch = scene.latentIm(row: min(h, row + patchSize - 1), ...
            col: min(w, col + patchSize - 1), :);
        gtPatch = digital_gt(row: min(h, row + patchSize - 1), ...
            col: min(w, col + patchSize - 1), :);

        scenePatch.latentIm = latentImPatch;
        scenePatch.L = scene.L;

        % Decide optimal binning size
        maxL = scene.L * max(latentImPatch(:));
        binning = round(double(bestBinningFunc(maxL)));
        binning = max(1, binning);
        % Decide gain
        gain = sensor.saturateVolt / (maxL + 2 * sqrt(maxL));
        gain = min(500, max(1, gain));
        [patch, nl_patch] = capture_with_ISP_binning( ...
            scenePatch, sensor, binning, gain);
        
        % Save binning and gain for visualization
        digital_ISP_varying(row: min(h, row + patchSize - 1), ...
            col: min(w, col + patchSize - 1), :) = patch;
        nl_digital_ISP_varying(row: min(h, row + patchSize - 1), ...
            col: min(w, col + patchSize - 1), :) = nl_patch;
        ISP_varying_gain(row: min(h, row + patchSize - 1), ...
            col: min(w, col + patchSize - 1)) = gain;
        ISP_varying_binning(row: min(h, row + patchSize - 1), ...
            col: min(w, col + patchSize - 1)) = binning;
    end
end
% Compute PSNR of the captured image
mse = norm(digital_ISP_varying(:) - digital_gt(:), 2) / ...
    sqrt(length(digital_gt(:)));
psnr = 20 * log10(maxVal / mse);
[ssimVal, digital_ISP_varying_ssimMap] = ssim(digital_ISP_varying / (2 ^ sensor.bit - 1), ...
                          digital_gt / (2 ^ sensor.bit - 1));

% Compute and visualize SNRs for varying patch brightness
[patchDRs, patchSNRs] = compute_snr_by_patch(digital_ISP_varying, ...
    digital_gt, patchSize);
figure(3), hold on; plot(patchDRs, patchSNRs, '*--', 'LineWidth', 2);
SNR_legend = [SNR_legend, 'ISP binning + SV gain']; 
legend(SNR_legend);
hold off;

% Visualize the captured image, binning, and gain map
figure(2), hold on;
subplot(4,3,10);
imshow(nl_digital_ISP_varying)
title(sprintf('ISP post processing + SV gain, psnr = %.2fdB, ssim = %.4f', psnr, ssimVal));

subf = subplot(4,3,11);
imagesc(ISP_varying_binning);
caxis(subf, [1, 10]);
title('ISP filter gain'); axis off;

subf = subplot(4,3,12);
imagesc(ISP_varying_gain);
caxis(subf, [1, 500]);
title('SV gain'); axis off;
hold off;

%% save images
saveDir = ['figures/images/', imName];
imwrite(nl_gt, [saveDir, '_gt.png']);
imwrite(nl_digital_1by1, [saveDir, '_no_binning_const_gain.png']);
imwrite(nl_digital_1by1_varying, [saveDir, '_no_binning_sv_gain.png']);
imwrite(nl_digital_additive, [saveDir, '_additive_binning_const_gain.png']);
imwrite(nl_digital_additive_varying, [saveDir, '_additive_binning_sv_gain.png']);
imwrite(nl_digital_average, [saveDir, '_average_binning_const_gain.png']);
imwrite(nl_digital_average_varying, [saveDir, '_average_binning_sv_gain.png']);
imwrite(nl_digital_ISP, [saveDir, '_ISP_binning_const_gain.png']);
imwrite(nl_digital_ISP_varying, [saveDir, '_ISP_binning_sv_gain.png']);