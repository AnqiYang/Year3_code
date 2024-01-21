function [digital, nl_digital] = capture_gt(scene, sensor, binning, gain)

gain = 1;
binning = 1;

voltage = scene.latentIm * scene.L;                 % get voltage
amplifiedVolt = gain * voltage;                     % amplifier gain
digital = amplifiedVolt / sensor.maxVolt * (2 ^ sensor.bit - 1);  % ADC conversion
digital = double(digital);

% normalize
digital = digital / 1 / gain;
% tonemapping
nl_digital = gammaCorrection(digital ./ max(digital(:)), ...
                sensor.gamma, sensor.offset, sensor.c);
