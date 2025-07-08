function spemOut = spem(x, y)
if ndims(x) == 3 && ndims(y) == 3
    % Handle 3D arrays: compute SPEM per layer
    numLayers = size(x, 3);
    spemOut = inf(numLayers,1);
    for n = 1:numLayers
        spemOut(n) = spem2D(x(:,:,n), y(:,:,n));
        if isnan(spemOut(n))
            spemOut(n) = 1 - hamming(x(:,:,n), y(:,:,n));
        end
    end
else
    % Handle 2D input
    spemOut = spem2D(x, y);
    if isnan(spemOut)
        spemOut = 1 - hamming(x, y);
    end
end
end

function spem = spem2D(x, y)
% Logical indexing for valid values
index = ~(isnan(x) | isnan(y));
validX = x(index);
validY = y(index);

% Alpha (Spearman Rank Correlation)
rs = corr(validY, validX, 'Type', 'Spearman');

% Beta (Coefficient of Variation Ratio)
meanY = mean(y(:), 'omitnan');
meanX = mean(x(:), 'omitnan');
stdY  = std(y(:), 'omitnan');
stdX  = std(x(:), 'omitnan');

if meanX == 0 || stdY == 0 || stdX == 0
    gamma = 0;
else
    gamma = (stdX / meanX) / (stdY / meanY);
end

% Gamma (Z-Score Difference)
zDiff = ((x - meanX) / stdX) - ((y - meanY) / stdY);
alpha = 1 - sqrt(mean(zDiff(:).^2, 'omitnan'));

% Bounded SPEM
spem = exp(-sqrt((rs - 1)^2 + (gamma - 1)^2 + (alpha - 1)^2));

end
