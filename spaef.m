function spaef = spaef(x,y)

index = ~(isnan(x) | isnan(y));
validX = x(index);
validY = y(index);

% Alpha
cc    = corrcoef(validY,validX); % Pearson correlation coefficient
alpha = cc(1,2);

% Beta
meanY = mean(y(:),'omitnan');
meanX = mean(x(:),'omitnan');
stdY  = std(y(:),'omitnan');
stdX  = std(x(:),'omitnan');
if meanX == 0 || stdY == 0 || stdX == 0
    beta = 0;
    % disp('WARNING: beta term in SPAEF was forced to be 0')
else
    beta  = (stdX/meanX) / (stdY/meanY); % cv ratio
end

% Gamma
% zScoreY      = zscore(y);
% zScoreX      = zscore(x);
zScoreY = (y - meanY) / stdY;
zScoreX = (x - meanX) / stdX;

bins = floor(sqrt(length(y)));

[yN,~] = histcounts(zScoreY,bins);
[xN,~] = histcounts(zScoreX,bins);

minOfHists     = min([yN; xN], [], 1);
overlappedHist = sum(minOfHists);
histogramMatch = overlappedHist/sum(yN);
gamma          = histogramMatch;

% SPAEF
spaef = 1 - sqrt((alpha-1)^2 + (beta-1)^2 + (gamma-1)^2);

end
