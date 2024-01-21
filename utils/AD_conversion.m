function [digital, saturation_map] = AD_conversion(amplifiedVolt, maxVolt, ADC_std, bit, blackLevel)

% Add post ADC read noise
amplifiedVolt = amplifiedVolt + randn(size(amplifiedVolt)) * ADC_std;

% % ADC saturation
% amplifiedVolt = min(maxVolt, max(0, amplifiedVolt));

% Scale to bit-depth
digital = amplifiedVolt / maxVolt * (2 ^ bit - 1);

% quantize
digital = round(digital);

% add constant offset
digital = digital + blackLevel;

% fprintf('After applying gain, %.2f percent pixels are saturated.\n', ...
%     100 * mean(digital(:) > (2 ^ bit - 1)));
% fprintf('After applying gain, %.2f percent pixels are cropped to zero.\n', ...
%     100 * mean(digital(:) < 0));
% saturation and floor
saturation_map = (digital >= 2 ^ bit - 1);
digital = min(2 ^ bit - 1, max(0, digital));

end