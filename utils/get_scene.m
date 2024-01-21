function scene = get_scene(dynamic_range)

if ~exist('dynamic_range', 'var')
    dynamic_range = 0;
end

if isa(dynamic_range, 'char')
    latentIm = exrread(['data/', dynamic_range]);
    latentIm = max(latentIm, 0);
    scene.latentIm = latentIm ./ max(latentIm(:));  % normalize to intensity 1
    scene.L = 800;
    [h, w, ~] = size(latentIm);

else
    scene.latentIm = zeros(1024, 1250, 3);
    scene.latentIm(:, 1:1024, :) = 2 ^ (dynamic_range) * 1;
    scene.latentIm(:, 1025:end, :) = 1;
    scene.L = 800;
    
end
