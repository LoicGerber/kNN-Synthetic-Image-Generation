function [di, partialTi, ki, sp] = createInputsMPS(lulcDir,mps,maskDir)

% CREATE MASK
% Mask with shape of image
mask = single(readgeoraster(maskDir));
sg = nan(size(mask));

% CREATE LATLON LAYERS
% Generate longitude (lon) and latitude (lat) maps
lon = zeros(size(sg));
lat = zeros(size(sg));
% Not real latlon, just info of where in image
lon_start = 0;  % Starting longitude
lat_start = 0;  % Starting latitude
% Populate lon and lat based on pixel positions
for la = 1:size(sg,1)
    for lo = 1:size(sg,2)
        lon(la,lo) = lon_start + lo;
        lat(la,lo) = lat_start + la;
    end
end
lat(mask==0) = nan;
lon(mask==0) = nan;
latlon = cat(3, lat, lon);

% LOAD LULC MAP
lulc = single(readgeoraster(lulcDir));
lulc(mask==0) = nan;

% CREATE PARTIAL TI
partialTi = cat(3, latlon, lulc);

% CREATE DI
di = cat(3, sg, partialTi);

% CREATE KI
ki = ones(mps.kernel_dims);
for w = 1:length(mps.kernel_weights)
    ki(:,:,w) = ki(:,:,w) .* mps.kernel_weights(w);
end

% CREATE SP
rng(mps.seed)
values = randperm(numel(mask));
sp = reshape(values, size(mask));
if mps.useMask == true
    sp(mask==0) = -inf;
end

end
