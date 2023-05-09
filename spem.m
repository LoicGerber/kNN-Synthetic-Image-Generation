function spem = spem(synImage,refImage)

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
    disp('WARNING: beta term in SPEM was forced to be 0')
else
    beta  = (stdSyn/meanSyn) / (stdRef/meanRef); % cv ratio
end

% Gamma
zScoreSyn = zscore(synImage);
zScoreRef = zscore(refImage);
gamma     = 1 - sqrt(immse(zScoreSyn,zScoreRef));

% SPEM
spem = 1 - sqrt((1-alpha)^2 + (1-beta)^2 + (1-gamma)^2);

end
