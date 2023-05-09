function spaef = spaef(synImage,refImage)

% Alpha
cc    = corrcoef(refImage,synImage); % Pearson correlation coefficient
alpha = cc(1,2);

% Beta
meanRef = mean(refImage(:));
meanSyn = mean(synImage(:));
stdRef  = std(refImage(:));
stdSyn  = std(synImage(:));
if meanSyn == 0 || stdRef == 0 || stdSyn == 0
    beta = 0;
    disp('WARNING: beta term in SPAEF was forced to be 0')
else
    beta  = (stdSyn/meanSyn) / (stdRef/meanRef); % cv ratio
end

% Gamma
zScoreRef      = zscore(refImage);
zScoreSyn      = zscore(synImage);

bins           = floor(sqrt(length(refImage)));

[refN,~]       = histcounts(zScoreRef,bins);
[synN,~]       = histcounts(zScoreSyn,bins);

minOfHists     = min([refN; synN], [], 1);
overlappedHist = sum(minOfHists);
histogramMatch = overlappedHist/sum(refN);
gamma          = histogramMatch;

% SPAEF
spaef = 1 - sqrt((alpha-1)^2 + (beta-1)^2 + (gamma-1)^2);

end
