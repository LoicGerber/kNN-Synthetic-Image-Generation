function [climateData,queryDates,learningDates,refValidation, ...
    Weights,sortedDates,synImages,validationMetric,sensitivityResults] = sensitivityAnalysis(rawData,nbImages_range,longWindow_range,inDir,outDir,targetVar,climateVars,normMethods,QdateStart,QdateEnd,LdateStart,LdateEnd,outputTime,targetDim, ...
    daysRange,metricKNN,useDOY,ensemble,generationType,parallelComputing,stochastic,stoSaveAll,metricV)

% Perform sensitivity analysis loop
iteration = 0;

% Initialize arrays to store sensitivity analysis results
sensitivityResults = cell(length(nbImages_range), length(longWindow_range));

for iRange = 1:length(nbImages_range)
    for jRange = 1:length(longWindow_range)
        iteration = iteration + 1;
        disp(['----- ITERATION ' num2str(iteration) ' -----'])
        % KNNDataGeneration
        longWindow  = longWindow_range(jRange);         % number of days to consider for the long climate window
        nbImages    = nbImages_range(iRange);         % K, number of days to consider for the generation of images

        % Validation switch
        if iRange == 1 && jRange == 1
            queryDates     = [];
            learningDates  = [];
            refValidation  = [];
            Weights        = [];
            climateData    = [];
            geoRef         = [];
        else
            sortedDates = [];
        end

        disp('--- 1. READING DATA ---')
        disp('Extracting climate informations...')
        climateData    = extractClimateData(climateVars,rawData,normMethods,QdateStart,QdateEnd,LdateStart,LdateEnd,longWindow,inDir,false);
        disp('Extracting Learning dates...')
        learningDates  = convertStructureToLearningDates(targetVar,LdateStart,LdateEnd,QdateStart,QdateEnd,rawData,climateData,targetDim,false,inDir,false);
        disp('Extracting Query dates...')
        [queryDates,learningDates,refValidation] = convertStructureToQueryDates(targetVar,targetDim,QdateStart,QdateEnd,learningDates,climateData,true,false,outputTime,inDir,useDOY);
        disp('Loading optimisedWeights.mat file...')
        Weights = createWeights(climateVars,metricKNN,useDOY,inDir);
        disp('--- 1. READING DATA DONE ---')

        disp('--- 2. KNN DATA SORTING ---')
        sortedDates = kNNDataSorting(climateVars,queryDates,learningDates,climateData,normMethods,0,longWindow,daysRange,Weights,nbImages,metricKNN,false,false,parallelComputing,inDir,false);
        disp('--- 2. KNN DATA SORTING DONE ---')

        disp('--- 3. SYNTHETIC IMAGES GENERATION ---')
        synImages = generateSynImages_ParamOptim(targetVar,targetDim,learningDates,sortedDates,nbImages,outDir,generationType,true,false,stochastic,stoSaveAll,ensemble);
        disp('--- 3. SYNTHETIC IMAGES GENERATION DONE ---')

        disp('--- 4. VALIDATION ---')
        validationMetric = validationMetrics(targetVar,targetDim,metricV,false,refValidation,synImages,stochastic,ensemble,outDir);
        disp('--- 4. VALIDATION DONE ---')

        % Store the sensitivity result in the array
        sensitivityResults(iRange, jRange) = {validationMetric.(lower(targetVar))(:,2)};
    end
end
save(fullfile(outDir,'sensitivityResults.mat'),'sensitivityResults');

% Visualise results
it     = 0;
nbIter = cumprod(size(sensitivityResults));
meanTS = nan(nbIter(2),2);
for i = 1:size(sensitivityResults, 1)
    for j = 1:size(sensitivityResults, 2)
        it           = it + 1;
        timeSeries   = sensitivityResults{i, j};
        meanTS(it,1) = mean(timeSeries);
        meanTS(it,2) = it;
    end
end

if metricV == 4 || metricV == 5
    sortedMeanTS = sortrows(meanTS,1,'descend');
else
    sortedMeanTS = sortrows(meanTS,1,'ascend');
end
best = sortedMeanTS(1,2);

metrics = {'MAE','RMSE','SPEM','SPAEF','KGE','NSE'};

it = 0;
bestDistImg = nan(size(sensitivityResults));
for i = 1:size(sensitivityResults, 1)
    for j = 1:size(sensitivityResults, 2)
        it = it + 1;
        if it == best
            disp(['Best parameter combination to optimise ' metrics{metricV} ':'])
            disp(['  k: ', num2str(nbImages_range(i))])
            disp(['  long: ', num2str(longWindow_range(j))])
            disp(['  Mean ' metrics{metricV} ': ', num2str(sortedMeanTS(1,1))])
            %break
        end
        bestDistImg(i,j) = meanTS(it);
    end
end

figure('WindowState', 'maximized');
hold on
imagesc(bestDistImg)
colormap(gca, turbo(256));
ylabel('\it k')
yticks(nbImages_range)
yticklabels({nbImages_range})
ylim([min(nbImages_range) max(nbImages_range)])
xlabel('Climate window length')
xticks(longWindow_range)
xticklabels({longWindow_range})
xlim([min(longWindow_range) max(longWindow_range)])
hcb=colorbar;
set(get(hcb,'label'),'string',['Mean ' metrics{metricV}],'Rotation',90);
set(gcf, 'color', 'white');
title('Best parameters combination')
saveas(gcf,strcat(outDir,'sensitivityAnalysis.png'))

end

