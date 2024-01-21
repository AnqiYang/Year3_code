function [voltage] = get_voltage(x, L, maxVolt, read_std, blackLevel)

% latent (clean) electron measurement
% x is normalized to have intensity 1
% voltage = poissrnd(x * L) + randn(size(x)) * read_std;
% voltage = voltage + blackLevel;
% voltage = min(maxVolt, max(0, voltage)); % full well capacity saturation

voltage = poissrnd(x * L);               % photon noise
voltage = min(maxVolt, voltage);         % full well saturation
voltage = voltage + randn(size(x)) * read_std; % readout noise

end