function [binnedVoltage] = binVoltage(voltage, binning, binningMode)

% Function to bin N x N pixels by summing or averaging them.
% Inputs:
%       voltage:       Analog voltage image [h, w] or [h, w, 3]
%       binning:       Binning size N. Bin N x N pixels together.
%       binningMode:   Additive or average.
% Outputs
%       binnedVoltage: immediately after binning, image size is [h/n, w/n];
%                      It is easier to compare results with the same size,
%                      so I kron and crop the final image back to [h, w].
%
%

% samplesPerPixel = round(p / dx);
sz = size(voltage);
new_sz = ceil(size(voltage) / binning) * binning;

% pad scene
padding = new_sz - sz;
prePadding = round(padding / 2);
postPadding = padding - round(padding / 2);
voltage = padarray(voltage, prePadding, 'replicate', 'pre');
voltage = padarray(voltage, postPadding, 'replicate', 'post');

% Binning n x n pixels 
numPixels = new_sz / binning;

sceneHat = mat2cell(voltage, binning*ones(1,numPixels(1)), ...
    binning*ones(1,numPixels(2)));
binnedVoltage = zeros(numPixels(1), numPixels(2));
for i = 1: numPixels(1)
    for j = 1: numPixels(2)
        binnedVoltage(i, j) = sum(sceneHat{i, j}, 'all');
    end
end
% binnedVoltage = kron(binnedVoltage, ones(binning));

% Implement binning mode (additive or average).
if strcmp('average', binningMode)
    binnedVoltage = binnedVoltage / (binning ^ 2);
elseif strcmp('additive', binningMode)
    binnedVoltage = binnedVoltage;
else
    error('invalid binning mode');
end

% % crop to orignal size
% binnedVoltage = binnedVoltage(prePadding(1) + 1: prePadding(1) + sz(1), ...
%                               prePadding(2) + 1: prePadding(2) + sz(2));
end