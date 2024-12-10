function spem = spem(x, y)
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

    % SPEM
    % spem = 1 - sqrt((rs - 1)^2 + (gamma - 1)^2 + (alpha - 1)^2);
    
    % Bounded SPEM
    spem = exp(-sqrt((rs - 1)^2 + (gamma - 1)^2 + (alpha - 1)^2));
end
