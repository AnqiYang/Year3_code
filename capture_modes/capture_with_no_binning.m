function [digital, nl_digital, saturation_map] = capture_with_no_binning(scene, sensor, ...
                                                         binning, gain)

voltage = get_voltage(scene.latentIm, scene.L, ...
                      sensor.maxVolt, sensor.read_std);
amplifiedVolt = gain .* voltage;
[digital, saturation_map] = AD_conversion(amplifiedVolt, sensor.maxVolt, ...
                        sensor.ADC_std, sensor.bit, sensor.blackLevel);
digital = double(digital);

% normalize digital to equivlant to gain=1, binning = 1
digital = digital - sensor.blackLevel;
digital = digital / 1 ./ gain;
digital = double(digital);

nl_digital = gammaCorrection(digital ./ (2 ^ sensor.bit - 1), ...
                sensor.gamma, sensor.offset, sensor.c);