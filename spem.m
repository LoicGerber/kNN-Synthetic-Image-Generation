function spem = spem(synImage,refImage)

% Alpha
index = ~(isnan(synImage) | isnan(refImage));
alpha = corr(refImage(index),synImage(index)); % Pearson correlation coefficient

% Beta
meanRef = mean(refImage(:),'omitnan');
meanSyn = mean(synImage(:),'omitnan');
stdRef  = std(refImage(:),'omitnan');
stdSyn  = std(synImage(:),'omitnan');
if meanSyn == 0 || stdRef == 0 || stdSyn == 0
    beta = 0;
    disp('WARNING: beta term in SPEM was forced to be 0')
else
    beta  = (stdSyn/meanSyn) / (stdRef/meanRef); % cv ratio
end

% Gamma
%zScoreSyn = zscore(synImage);
zScoreSyn = (synImage - meanSyn)/stdSyn;
%zScoreRef = zscore(refImage);
zScoreRef = (refImage - meanRef)/stdRef;
gamma     = 1 - sqrt(mean(mean((zScoreSyn-zScoreRef).^2,'omitnan'),'omitnan'));

% SPEM
spem = 1 - sqrt((1-alpha)^2 + (1-beta)^2 + (1-gamma)^2);

end
