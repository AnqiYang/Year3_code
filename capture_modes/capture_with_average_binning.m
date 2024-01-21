function [digital, nl_digital] = capture_with_average_binning(scene, ...
                                                    sensor, binning, gain)

[h, w, ~] = size(scene.latentIm);
voltage = get_voltage(scene.latentIm, scene.L, ...
                    sensor.maxVolt, sensor.read_std);
for cc = 1: 3
    binnedVoltage(:,:,cc) = binVoltage(voltage(:,:,cc), binning, 'average');
end
amplifiedVolt = gain * binnedVoltage;
digital = AD_conversion(amplifiedVolt, sensor.maxVolt, ...
    sensor.ADC_std, sensor.bit, sensor.blackLevel);
digital = double(digital);
for cc = 1: 3
    resizeDigital(:,:,cc) = resize(digital(:,:,cc), binning, [h, w]);
end
digital = resizeDigital;

% normalize
digital = digital - sensor.blackLevel;
digital = digital / gain;
digital = double(digital);

nl_digital = gammaCorrection(digital ./ (2 ^ sensor.bit - 1), ...
                sensor.gamma, sensor.offset, sensor.c);