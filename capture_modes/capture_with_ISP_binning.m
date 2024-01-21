function [digital_ISP, nl_digital] = capture_with_ISP_binning(scene, ...
                                                    sensor, binning, gain)


[h, w, ~] = size(scene.latentIm);
voltage = get_voltage(scene.latentIm, scene.L, ...
                    sensor.maxVolt, sensor.read_std);
amplifiedVolt = gain * voltage;
digital = AD_conversion(amplifiedVolt, sensor.maxVolt, ...
    sensor.ADC_std, sensor.bit, sensor.blackLevel);
digital = double(digital);

% normalize
digital = digital - sensor.blackLevel;
digital = digital / gain;

% ISP binning
for cc = 1: 3
    binnedDigital(:,:,cc) = binVoltage(digital(:,:,cc), binning, 'average');
end
for cc = 1: 3
    digital_ISP(:,:,cc) = resize(binnedDigital(:,:,cc), binning, [h, w]);
end
digital_ISP = double(digital_ISP);
% tonemapping
nl_digital = gammaCorrection(digital_ISP ./ (2 ^ sensor.bit - 1), ...
                sensor.gamma, sensor.offset, sensor.c);