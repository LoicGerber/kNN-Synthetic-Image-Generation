function [geoRef,climateData,queryDates,learningDates,refValidation,additionalVars,Weights,sortedDates,synImages,validationMetric,optimisedWeights] = MAIN(...
    rawDir,inputDir,outputDir,var,vars,addVars,QdateStart,QdateEnd,LdateStart,LdateEnd,outputTime,precision, ...
    shortWindow,longWindow,nbImages,ensemble,GenerationType,OutputType,coordRefSysCode,parallelComputing, ...
    NetCDFtoInputs,createGenWeights,KNNsorting,generateImage,bootstrap,validationPrep,validation,metricViz, ...
    metric,optimPrep,saveOptimPrep,optimisation,nbOptiRuns)

%% Setup
close all
poolobj = gcp('nocreate');
delete(poolobj);
tStart = tic;

%% Reading the data needed for ranking learning dates using "KNNDataSorting" Function
disp('--- 1. READING DATA ---')

if NetCDFtoInputs == true || optimPrep == true || validationPrep == true
    disp('Formatting input data...')
    rawData        = ConvertNetCDFtoStructure(var,vars,addVars,precision,rawDir,inputDir);
    disp('Extracting georeference informations...')
    geoRef         = extractGeoInfo(var,coordRefSysCode,rawDir,inputDir);
    disp('Extracting climate informations...')
    climateData    = extractClimateData(vars,rawData,QdateStart,QdateEnd,LdateStart,LdateEnd,longWindow,inputDir);
    disp('Extracting Learning dates...')
    learningDates  = ConvertStructureToLearningDates(var,LdateStart,LdateEnd,QdateStart,QdateEnd,rawData,climateData,optimPrep,inputDir);
    disp('Extracting Query dates...')
    [queryDates,learningDates,refValidation] = ConvertStructureToQueryDates(var,QdateStart,QdateEnd,learningDates,climateData,longWindow,validationPrep,optimPrep,outputTime,inputDir,outputDir);
    disp('Extracting additional variables...')
    additionalVars = extractAdditionalVars(addVars,rawData,QdateStart,QdateEnd,LdateStart,LdateEnd,inputDir);
elseif NetCDFtoInputs == false && validationPrep == false
    disp('Loading QueryDates.mat file...')
    queryDates     = load(fullfile(inputDir,'queryDates.mat'));
    queryDates     = queryDates.queryDates;
    disp('Loading LearningDates.mat file...')
    learningDates  = load(fullfile(inputDir,'learningDates.mat'));
    learningDates  = learningDates.learningDates;
    disp('Loading climateData.mat file...')
    climateData    = load(fullfile(inputDir,'climateData.mat'));
    climateData    = climateData.climateData;
    disp('Loading additionalVars.mat file...')
    additionalVars = load(fullfile(inputDir,'additionalVars.mat'));
    additionalVars = additionalVars.additionalVars;
    disp('Loading GeoRef.mat file...')
    geoRef         = load(fullfile(inputDir,'GeoRef.mat'));
    geoRef         = geoRef.geoRef;
    if optimisation == true || validation == true || metricViz == true
        disp('Loading refValidation.mat file...')
        refValidation = load(fullfile(inputDir,'refValidation.mat'));
        refValidation = refValidation.refValidation;
    else
        refValidation = [];
    end
end
if createGenWeights == true || optimPrep == true || optimisation == true
    disp('Creating generic weights...')
    Weights = createWeights(var,vars,addVars,inputDir);
elseif createGenWeights == false
    disp('Loading optimisedWeights.mat file...')
    optimisedWeights = load(fullfile(inputDir,'optimisedWeights.mat'));
    Weights = optimisedWeights.optimisedWeights;
end

disp('--- 1. READING DATA DONE ---')

%% The function for Ranking the learning dates
disp('--- 2. KNN DATA SORTING ---')

% Generate ranked Learning Dates for each Query Date
if KNNsorting == true || validationPrep == true || optimPrep == true
    sortedDates = KNNDataSorting(var,vars,addVars,queryDates,learningDates,climateData,additionalVars,shortWindow,longWindow,Weights,nbImages,optimPrep,saveOptimPrep,parallelComputing,inputDir);
elseif KNNsorting == false && validationPrep == false && (optimPrep == false && optimisation == false)
    disp('Loading sortedDates.mat file...')
    sortedDates = load(fullfile(inputDir,'KNNSorting.mat'));
    sortedDates = sortedDates.sortedDates;
elseif optimPrep == false && optimisation == true
    disp('Loading KNNDistances.mat file...')
    sortedDates = load(fullfile(inputDir,'KNNDistances.mat'));
    sortedDates = sortedDates.sortedDates;
end

disp('--- 2. KNN DATA SORTING DONE ---')

%% Generation of Synthetic Images
disp('--- 3. SYNTHETIC IMAGES GENERATION ---')

if (generateImage == true || validation == true) && optimisation == false
    synImages = GenerateSynImages(var,learningDates,sortedDates,geoRef,outputDir,GenerationType,validation,optimisation,bootstrap,ensemble,OutputType);
elseif optimisation == true && validation == false
    disp('Optimisation run, synthetic image generation skipped...')
    synImages = [];
elseif metricViz == true
    disp('Loading synValidation.mat file...')
    synImages = load(fullfile(outputDir,'synValidation.mat'));
    synImages = synImages.synImages;
elseif generateImage == false && validation == false
    disp('Synthetic image generation skipped...')
    synImages = [];
end

disp('--- 3. SYNTHETIC IMAGES GENERATION DONE ---')

%% Validation
if (validation == true || metricViz == true) && optimisation == false
    disp('--- 4. VALIDATION ---')
    
    validationMetric = validationMetrics(var,metric,optimisation,refValidation,synImages,bootstrap,ensemble,outputDir);
    visualiseMetrics(var,refValidation,synImages,validationMetric,metric,LdateStart,LdateEnd,QdateStart,QdateEnd,bootstrap,outputDir);
    
    disp('--- 4. VALIDATION DONE ---')
else
    validationMetric = [];
end

%% Optimisation
if optimisation == true
    disp('--- 4. OPTIMISATION ---')
    
    if saveOptimPrep == true
        sortedDates = [];
    end
    
    % Get the table variable names
    variableNames = string(Weights.Properties.VariableNames);
    % Iterate over each variable
    for i = 1:numel(variableNames)
        bayesWeights(i) = optimizableVariable(variableNames{i}, [0, 1], 'Type', 'real');
    end
    % Set up the Bayesian optimization
    fun = @(x)computeObjectiveOptim(x.(1), x.(2), x.(3), x.(4), x.(5), x.(6), x.(7), x.(8), x.(9), ...
        var, addVars, learningDates, sortedDates, refValidation, saveOptimPrep, nbImages, ...
        geoRef, GenerationType, bootstrap, ensemble, metric, validation, optimisation, inputDir, outputDir);
    % Run the Bayesian optimization
    %if parallelComputing == true
    %    results = bayesopt(fun,bayesWeights,'Verbose',0,'AcquisitionFunctionName','expected-improvement-plus','MaxObjectiveEvaluations',nbOptiRuns,'UseParallel',true);
    %else
    results = bayesopt(fun,bayesWeights,'Verbose',0,'AcquisitionFunctionName','expected-improvement-plus','MaxObjectiveEvaluations',nbOptiRuns);
    %end
    % Retrieve the optimal weights
    disp('  Saving optimisedWeights.mat...')
    optimisedWeights = results.XAtMinObjective;
    optimisedWeights = array2table(table2array(optimisedWeights) ./ sum(table2array(optimisedWeights)),'VariableNames', results.XAtMinObjective.Properties.VariableNames)
    save(fullfile(inputDir,'optimisedWeights.mat'), 'optimisedWeights', '-v7.3','-nocompression');
    
    disp('--- 4. OPTIMISATION DONE ---')
else
    optimisedWeights = [];
end

tEnd = toc(tStart);

disp(['--- FINISHED IN ' num2str(tEnd) ' SECONDS ---'])

end

