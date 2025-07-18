function [geoRef,climateData,queryDates,learningDates,refValidation, ...
    Weights,sortedDates,synImages,validationMetric,sensitivityResults,optimisedWeights] = MAIN(...
    rawDir,outputDir,optiWeightsDir,maskDir,lulcDir,targetVar,climateVars,normMethods,QdateStart,QdateEnd,LdateStart,LdateEnd,outputTime,targetDim,saveMats, ...
    useDOY,shortWindow,longWindow,daysRange,nbImages,metricKNN,ensemble,generationType,mps,outputType,coordRefSysCode,parallelComputing, ...
    netCDFtoInputs,createGenWeights,kNNsorting,generateImage,saveNetCDF,stochastic,stoSaveAll,validationPrep,validation,pixelWise,createGIF, ...
    metricViz,metricV,nanValue,varLegend,varRange,errRange,sensiAnalysis,nbImages_range,longWindow_range,optimPrep,saveOptimPrep,optimisation,nbOptiRuns)

%% Setup
close all
poolobj = gcp('nocreate');
delete(poolobj);
tStart = tic;

if LdateStart > LdateEnd
    error('Learning date start > learning date end')
elseif QdateStart > QdateEnd
    error('Query date start > Query date end')
end

inDir = strcat(outputDir,'\inputData');
outDir = strcat(outputDir,'\output');
% Check if output directories exist, if not create them
if ~exist(inDir, 'dir')
    mkdir(inDir)
end
if ~exist(outDir, 'dir')
    mkdir(outDir)
end

%% Reading the data needed for ranking learning dates using "KNNDataSorting" Function
if sensiAnalysis == false
    disp('--- 1. READING DATA ---')
    
    if netCDFtoInputs == true || optimPrep == true || validationPrep == true
        disp('Formatting input data...')
        rawData        = convertRawDataToStructure(targetVar,targetDim,climateVars,rawDir);
        disp('Extracting georeference informations...')
        if targetDim ~= 1
            geoRef = extractGeoInfo(targetVar,coordRefSysCode,rawDir,inDir);
        else
            geoRef = [];
        end
        disp('Extracting climate informations...')
        climateData    = extractClimateData(climateVars,rawData,normMethods,QdateStart,QdateEnd,LdateStart,LdateEnd,longWindow,inDir,saveMats);
        disp('Extracting Learning dates...')
        learningDates  = convertStructureToLearningDates(targetVar,LdateStart,LdateEnd,QdateStart,QdateEnd,rawData,climateData,targetDim,optimPrep,inDir,saveMats);
        disp('Extracting Query dates...')
        [queryDates,learningDates,refValidation] = convertStructureToQueryDates(targetVar,targetDim,QdateStart,QdateEnd,learningDates,climateData,validationPrep,optimPrep,outputTime,inDir,saveMats);
    elseif netCDFtoInputs == false && validationPrep == false
        disp('Loading QueryDates.mat file...')
        queryDates     = load(fullfile(inDir,'queryDates.mat'));
        queryDates     = queryDates.queryDates;
        disp('Loading LearningDates.mat file...')
        learningDates  = load(fullfile(inDir,'learningDates.mat'));
        learningDates  = learningDates.learningDates;
        disp('Loading climateData.mat file...')
        climateData    = load(fullfile(inDir,'climateData.mat'));
        climateData    = climateData.climateData;
        disp('Loading GeoRef.mat file...')
        if targetDim ~= 1
            geoRef         = load(fullfile(inDir,'GeoRef.mat'));
            geoRef         = geoRef.geoRef;
        else
            geoRef = [];
        end
        if optimisation == true || validation == true || metricViz == true
            disp('Loading refValidation.mat file...')
            refValidation = load(fullfile(inDir,'refValidation.mat'));
            refValidation = refValidation.refValidation;
        else
            refValidation = [];
        end
    end
    if createGenWeights == true || optimPrep == true || optimisation == true
        disp('Creating generic weights...')
        Weights = createWeights(climateVars,metricKNN,useDOY,inDir);
    elseif createGenWeights == false
        disp('Loading optimisedWeights.mat file...')
        optimisedWeights = load(optiWeightsDir);
        Weights = optimisedWeights.optimisedWeights;
    end
    
    % TEMPORARY FIX FOR PIXELWISE SIZE MISMATCH
%     learningDatesMod = learningDates;
%     for row = 1:size(learningDatesMod, 1)
%         originalMatrix = cell2mat(learningDatesMod{row, 2});
%         topNaNRow = NaN(1, size(originalMatrix, 2));
%         matrixWithTopNaN = [topNaNRow; originalMatrix];
%         rightNaNColumn = NaN(size(matrixWithTopNaN, 1), 1);
%         matrixWithTopAndRightNaN = [matrixWithTopNaN, rightNaNColumn];
%         learningDatesMod{row, 2} = {matrixWithTopAndRightNaN};
%     end
%     learningDates = learningDatesMod;
%     refValidationMod = refValidation.e;
%     for row = 1:size(refValidationMod,3)
%         ogMatRef = refValidationMod(:,:,row);
%         topNaNRow = NaN(1, size(ogMatRef, 2));
%         matrixWithTopNaN = [topNaNRow; ogMatRef];
%         rightNaNColumn = NaN(size(matrixWithTopNaN, 1), 1);
%         matrixWithTopAndRightNaN = [matrixWithTopNaN, rightNaNColumn];
%         newRefValidationMod(:,:,row) = matrixWithTopAndRightNaN;
%     end
%     refValidation.e = newRefValidationMod;

    disp('--- 1. READING DATA DONE ---')
    
    %% The function for Ranking the learning dates
    disp('--- 2. KNN DATA SORTING ---')
    
    % Generate ranked Learning Dates for each Query Date
    if kNNsorting == true || validationPrep == true || optimPrep == true
        if pixelWise == false
            sortedDates = kNNDataSorting(climateVars,queryDates,learningDates,climateData,normMethods,shortWindow,longWindow,daysRange,Weights,nbImages,metricKNN,optimPrep,saveOptimPrep,parallelComputing,inDir,useDOY);
        else
            sortedDates = pixelWise_kNNDataSorting(maskDir,climateVars,queryDates,learningDates,climateData,longWindow,daysRange,nbImages,metricKNN,optimPrep,saveOptimPrep,parallelComputing,inDir);
        end
    elseif kNNsorting == false && validationPrep == false && (optimPrep == false && optimisation == false)
        disp('Loading sortedDates.mat file...')
        sortedDates = load(fullfile(inDir,'KNNSorting.mat'));
        sortedDates = sortedDates.sortedDates;
    elseif optimPrep == false && optimisation == true
        disp('Loading KNNDistances.mat file...')
        sortedDates = load(fullfile(inDir,'KNNDistances.mat'));
        sortedDates = sortedDates.sortedDates;
    end
    
    disp('--- 2. KNN DATA SORTING DONE ---')
    
    %% Generation of Synthetic Images
    disp('--- 3. SYNTHETIC IMAGES GENERATION ---')
    
    if (generateImage == true && validation == true) && optimisation == false
        if pixelWise == false
            synImages = generateSynImages(maskDir,targetVar,targetDim,learningDates,sortedDates,mps,lulcDir,geoRef,outDir,generationType,validation,saveNetCDF,stochastic,stoSaveAll,nbImages,ensemble,outputType);
        else
            synImages = pixelWise_generateSynImages(maskDir,targetVar,learningDates,sortedDates,geoRef,outDir,generationType,validation,optimisation,stochastic,stoSaveAll,nbImages,ensemble,outputType);
        end
    elseif (optimisation == true || optimPrep == true) && validation == false
        disp('Optimisation run, synthetic image generation skipped...')
        synImages = [];
    elseif metricViz == true
        disp('Loading synValidation.mat file...')
        synImages = load(fullfile(outDir,'synValidation.mat'));
        synImages = synImages.synImages;
    elseif generateImage == false && validation == false
        disp('Synthetic image generation skipped...')
        synImages = [];
    end
    
    disp('--- 3. SYNTHETIC IMAGES GENERATION DONE ---')
    
    %% Validation
    if (validation == true || metricViz == true) && optimisation == false
        disp('--- 4. VALIDATION ---')
        
        validationMetric = validationMetrics(targetVar,targetDim,metricV,optimisation,refValidation,synImages,stochastic,ensemble,outDir);
        visualiseMetrics(nbImages,pixelWise,targetVar,climateVars,targetDim,refValidation,synImages,validationMetric,sortedDates,climateData,metricV,nanValue,varLegend,varRange,errRange,metricKNN,LdateStart,LdateEnd,QdateStart,QdateEnd,daysRange,outputTime,stochastic,outDir,createGIF);

        disp('--- 4. VALIDATION DONE ---')
    else
        validationMetric = [];
    end
end
%% Sensitivity analysis
if sensiAnalysis == true
    disp('--- SENSITIVITY ANALYSIS')
    disp('Formatting input data...')
    rawData = convertRawDataToStructure(targetVar,targetDim,climateVars,rawDir);
    geoRef  = [];
    
    [climateData,queryDates,learningDates,refValidation, ...
    Weights,sortedDates,synImages,validationMetric,sensitivityResults] = sensitivityAnalysis(rawData,nbImages_range,longWindow_range,inDir,outDir,targetVar,climateVars,normMethods,QdateStart,QdateEnd,LdateStart,LdateEnd,outputTime,targetDim, ...
                                             daysRange,metricKNN,useDOY,ensemble,generationType,parallelComputing,stochastic,stoSaveAll,metricV);
    optimisedWeights = [];

    disp('--- SENSITIVITY ANALYSIS DONE')
else
    sensitivityResults = [];
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
        initialW.(variableNames{i}) = 1;
    end
    initialW = struct2table(initialW);
    % Set up the Bayesian optimization
    fun = @(x)computeObjectiveOptim(x.(1), x.(2), x.(3), x.(4), x.(5), x.(6), x.(7), x.(8), ...
        targetVar, targetDim, learningDates, sortedDates, refValidation, saveOptimPrep, metricKNN, nbImages, ...
        generationType, metricV, optimisation, useDOY, inDir, outDir);
    % Run the Bayesian optimization
    %if parallelComputing == true
    %    results = bayesopt(fun,bayesWeights,'Verbose',0,'AcquisitionFunctionName','expected-improvement-plus','MaxObjectiveEvaluations',nbOptiRuns,'UseParallel',true);
    %else
    results = bayesopt(fun,bayesWeights(1:8),'Verbose',0,'AcquisitionFunctionName','expected-improvement-plus','MaxObjectiveEvaluations',nbOptiRuns,'InitialX',initialW(:,1:8));
    %end
    % Retrieve the optimal weights
    disp('  Saving optimisedWeights.mat...')
    optimisedWeights = results.XAtMinObjective;
    % Normalize weights (variables and metrics separately)
    optimisedWeightsArray = table2array(optimisedWeights);
    variablesWeights = optimisedWeightsArray(1:7);
    metricsWeights   = optimisedWeightsArray(8);
    helWeight        = 1 - metricsWeights;
    metricsWeights   = [metricsWeights helWeight];
    
    variablesWeights = variablesWeights / sum(variablesWeights); % Normalize variables weights
    
    optimisedWeightsNormalized = [variablesWeights, metricsWeights]; % Combine normalized weights
    
    % Convert to table and save
    optimisedWeightsTable = array2table(optimisedWeightsNormalized, 'VariableNames', variableNames)
    save(fullfile(inDir, 'optimisedWeights.mat'), 'optimisedWeightsTable', '-v7.3', '-nocompression');

    disp('--- 4. OPTIMISATION DONE ---')
else
    optimisedWeights = [];
end

tEnd = toc(tStart);

disp(['--- FINISHED IN ' num2str(tEnd) ' SECONDS ---'])

end

