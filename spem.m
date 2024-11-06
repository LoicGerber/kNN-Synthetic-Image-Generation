function spem = spem(synImage,refImage)

% Alpha
index = ~(isnan(synImage) | isnan(refImage));
rs = spear(refImage(index),synImage(index)); % Spearman Rank Correlation

% Beta
meanRef = mean(refImage(:),'omitnan');
meanSyn = mean(synImage(:),'omitnan');
stdRef  = std(refImage(:),'omitnan');
stdSyn  = std(synImage(:),'omitnan');
if meanSyn == 0 || stdRef == 0 || stdSyn == 0
    theta = 0;
%     fprintf('\n    WARNING: beta term in SPEM was forced to be 0')
else
    theta = (stdSyn/meanSyn) / (stdRef/meanRef); % cv ratio
end

% Gamma
%zScoreSyn = zscore(synImage);
zScoreSyn = (synImage - meanSyn)/stdSyn;
%zScoreRef = zscore(refImage);
zScoreRef = (refImage - meanRef)/stdRef;
alpha     = 1 - sqrt(mean(mean((zScoreSyn-zScoreRef).^2,'omitnan'),'omitnan'));

% SPEM
spem = 1 - sqrt((rs-1)^2 + (theta-1)^2 + (alpha-1)^2);

end
