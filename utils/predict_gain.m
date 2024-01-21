function [next_gain] = predict_gain(curr_digital, curr_saturated, sensor, mode)

ADC_scale = (2 ^ sensor.bit - 1) / sensor.maxVolt;
estimatedL = max(2e-4, curr_digital / ADC_scale);
switch mode
    case 1
        next_gain = sensor.saturateVolt ./ (estimatedL + 3 * sqrt(estimatedL));
        % reset saturated pixel to smallest gain
        next_gain(curr_saturated) = 1;
    case 2
        next_gain = 0.6 * sensor.saturateVolt ./ estimatedL;
        % reset saturated pixel to smallest gain
        next_gain(curr_saturated) = 1;
    case 3
        next_gain = sensor.saturateVolt ./ (estimatedL + 4 * sqrt(estimatedL));
        % reset saturated pixel to smallest gain
        next_gain(curr_saturated) = 1;
    otherwise
        error('mode not implemented.\n');
end

% round up to bins in decibel
gain_db = 10*log10(next_gain);
gain_db = floor((gain_db * 10) / 2) * 2 / 10;
next_gain = 10 .^ (gain_db / 10);

% crop gain to reasonable value
next_gain = min(500, max(1, next_gain));