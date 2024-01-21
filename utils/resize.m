function [upVolt] = resize(dnVolt, binning, upSize)

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

% pad scene
sz = upSize;
new_sz = ceil(upSize / binning) * binning;

% pad scene
padding = new_sz - sz;
prePadding = round(padding / 2);


upVolt = kron(dnVolt, ones(binning));

% crop to orignal size
upVolt = upVolt(prePadding(1) + 1: prePadding(1) + upSize(1), ...
                prePadding(2) + 1: prePadding(2) + upSize(2));
end