function [nl_rgb] = gammaCorrection(lin_rgb, gamma, offset, c)
nl_rgb = max(0, c*lin_rgb-offset) .^(gamma);
end

% function [nl_rgb] = gammaCorrection(lin_rgb, offset)
% nl_rgb = lin_rgb ./ (lin_rgb + offset);
% end