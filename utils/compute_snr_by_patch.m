function [DRs, histSNRs] = compute_snr_by_patch(im, imRef, patchSize)

% Compute
[h, w, ~] = size(im);
numPatch = ceil(size(im, 1) / patchSize) * ceil(size(im, 2) / patchSize);
patchDRs = zeros(numPatch, 1);
patchMSEs = zeros(numPatch, 1);

maxVal = max(imRef(:));
count = 0;
for row = 1: patchSize: size(im, 1)
    for col = 1: patchSize: size(im, 2)

        count = count + 1;
        imPatch = im(row: min(h, row + patchSize - 1), ...
            col: min(w, col + patchSize - 1), :);
        refPatch = imRef(row: min(h, row + patchSize - 1), ...
            col: min(w, col + patchSize - 1), :);

        meanVal = mean(refPatch(:));
        % mse = norm(imPatch(:) - refPatch(:), 2) / sqrt(length(refPatch(:)));
        mse = sum((imPatch(:) - refPatch(:)) .^ 2) / (length(refPatch(:)));
        
        patchDRs(count) = log2(meanVal / maxVal);
        patchMSEs(count) = mse;
        patchSNRs(count) = meanVal / sqrt(mse);
    end
end

patchDRs = patchDRs(1: count);
patchMSEs = patchMSEs(1: count);
patchSNRs = patchSNRs(1: count);

binWidth = 1;
DRs = -15+binWidth:binWidth:0;
histMSEs = zeros(1, length(DRs));
histCounts = zeros(1, length(DRs));


for i = 1: count
    bin = round((patchDRs(i)+15+binWidth)/binWidth);
    if bin < 1 || bin > length(DRs); continue; end
    histMSEs(bin) = histMSEs(bin) + patchMSEs(i);
    histCounts(bin) = histCounts(bin) + 1;
end

DRs = DRs(histCounts > 0);
histMSEs = histMSEs(histCounts > 0);
histCounts = histCounts(histCounts > 0);

histMSEs = histMSEs ./ histCounts;
Ls = 2 .^ DRs * maxVal;
histSNRs = Ls ./ sqrt(histMSEs);