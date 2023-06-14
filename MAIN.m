%% SYNTHETIC IMAGE GENERATION USING A KNN ALGORITHM

% Created 04.2023 by Loïc Gerber

%
%
%
% REDO DOCUMENTATION
%
%
%

%% PARAMETERS DEFIINITION
% --- IMPORTANT ---
% Input files and their variable must be named as the variables defined below (pre.nc, 'pre', tmax.nc, 'tmax', etc.)
% --- IMPORTANT ---

% Setup
clc, clear, close all
poolobj = gcp('nocreate');
delete(poolobj);
tStart = tic;

% All directories
rawDir    = 'C:\Users\loger\OneDrive - Université de Lausanne\Documents\PhD\knn_image_generation\syntheticImageGeneration\voltaData\';           % Path to raw data
inputDir  = 'C:\Users\loger\OneDrive - Université de Lausanne\Documents\PhD\knn_image_generation\syntheticImageGeneration\test\inputData\';      % Path to saved input data
outputDir = 'C:\Users\loger\OneDrive - Université de Lausanne\Documents\PhD\knn_image_generation\syntheticImageGeneration\test\output\';         % Path to results

% ConvertStructureToInputs
var               = "Et";                          % Variable to be generated, with "example"
vars              = ["Tavg","Tmin","Tmax","Pre"];  % Input variables considered for the data generation, with ["example1","example2"]
addVars           = [];                            % Additional input variables, with ["example1","example2"]
QdateStart        = 20000601;                      % YYYYMMDD - Start of the Generation period
QdateEnd          = 20000610;                      % YYYYMMDD - End of the Generation period
LdateStart        = 20000101;                      % YYYYMMDD - Start of the Learning period
LdateEnd          = 20001231;                      % YYYYMMDD - End of the Learning period
outputTime        = 1;                             % Image generation timestep: 1 = DAILY, 2 = MONTHLY
precision         = 1;                             % Precision needed, 1 = single, 2 = double

% KNNDataGeneration
shortWindow       = 5;        % number of days to consider for the short climate window
longWindow        = 30;       % number of days to consider for the long climate window
nbImages          = 10;       % K, number of days to consider for the generation of images

% GenerateSynImages
ensemble          = 20;       % when using bootstrap, number of ensembles created
GenerationType    = 2;        % data generation type,  1 = BINARY,  2 = MEAN OF SELECTED IMAGES, 3 = MEDIAN OF SELECTED IMAGES
OutputType        = 1;        % output data file type, 1 = GeoTIFF, 2 = individual NetCDF files
coordRefSysCode   = 4326;     % Coordinate reference system code, WGS84 = 4326, https://epsg.org/home.html

% Functions switches
parallelComputing = false;    % true = parallel computing ON,  false = parallel computing OFF
NetCDFtoInputs    = false;    % true = create inputs,          false = load inputs
createOptiWeights = false;    % true = create generic weights, false = load optimised weights
KNNsorting        = false;    % true = create sorted data,     false = load sorted data
generateImage     = false;    % true = image generation ON,    false = image generation OFF
bootstrap         = false;    % true = bootstrap ON,           false = bootstrap OFF

% Validation switch
validationPrep    = false;    % true = validation preparation ON,    false = validation preparation OFF (!!! BYPASSES PREVIOUS SWITCHES !!!)
validation        = false;    % true = validation ON,    false = validation OFF (!!! BYPASSES PREVIOUS SWITCHES !!!)
metricViz         = false;    % true = visualisation ON, false = visualisation OFF
metric            = 1;        % 1 = RMSE, 2 = SPEM, 3 = SPAEF, 4 = Symmetric Phase-only Matched Filter-based Absolute Error Function (SPOMF)

% Bayesian optimisation switch
optimPrep         = true;    % true = optimisation preparation ON, false = optimisation preparation OFF (!!! BYPASSES PREVIOUS SWITCHES !!!)
saveOptimPrep     = false;
optimisation      = true;     % true = optimisation ON, false = optimisation OFF (!!! run AFTER optimisation preparation !!!)
nbOptiRuns        = 5;        % Number of runs for the Bayesian optimisation

%% Reading the data needed for ranking learning dates using "KNNDataSorting" Function
disp('--- 1. READING DATA ---')

if NetCDFtoInputs == true || optimPrep == true || validationPrep == true
    disp('Formatting input data for production/weights optimisation/validation run...')
    rawData                    = ConvertNetCDFtoStructure(var,vars,addVars,precision,rawDir,inputDir);
    GeoRef                     = extractGeoInfo(var,coordRefSysCode,rawDir,inputDir);
    climateData                = extractClimateData(vars,rawData,QdateStart,QdateEnd,LdateStart,LdateEnd,longWindow,inputDir);
    learningDates              = ConvertStructureToLearningDates(var,LdateStart,LdateEnd,QdateStart,QdateEnd,rawData,climateData,optimPrep,inputDir);
    [queryDates,learningDates] = ConvertStructureToQueryDates(var,QdateStart,QdateEnd,learningDates,climateData,longWindow,GeoRef,validationPrep,optimPrep,outputTime,inputDir,outputDir);
    additionalVars             = extractAdditionalVars(addVars,rawData,QdateStart,QdateEnd,LdateStart,LdateEnd,inputDir);
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
    GeoRef         = load(fullfile(inputDir,'GeoRef.mat'));
    GeoRef         = GeoRef.GeoRef;
end
if createOptiWeights == true || optimPrep == true || optimisation == true
    disp('Creating generic Weights.mat file...')
    Weights = createWeights(var,vars,addVars,inputDir);
elseif createOptiWeights == false
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
elseif KNNsorting == false && validationPrep == false && optimPrep == false
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

if generateImage == true && validation == false && optimisation == false
    GenerateSynImages(var,learningDates,sortedDates,GeoRef,outputDir,GenerationType,optimisation,bootstrap,ensemble,OutputType);
elseif validation == true
    OutputType = 1;
    GenerateSynImages(var,learningDates,sortedDates,GeoRef,outputDir,GenerationType,optimisation,bootstrap,ensemble,OutputType);
elseif optimisation == true
    disp('Optimisation run, synthetic image generation skipped...')
elseif generateImage == false && validation == false
    disp('Synthetic image generation skipped...')
end

disp('--- 3. SYNTHETIC IMAGES GENERATION DONE ---')

%% Validation
if (validation == true || metricViz == true) && optimisation == false
    disp('--- 4. VALIDATION ---')
    
    validationMetric = validationMetrics(metric,optimisation,outputDir);
    visualiseMetrics(validationMetric,metric,LdateStart,LdateEnd,outputDir);
    
    disp('--- 4. VALIDATION DONE ---')
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
    bayesWeights = [];
    for i = 1:numel(variableNames)
        %bayesWeights{i} = optimizableVariable(variableNames{i},[0,1],'Type','real');
        eval(sprintf('%s = optimizableVariable(variableNames{i}, [0,1], ''Type'', ''real'');', variableNames{i}));
        eval(sprintf('bayesWeights = [bayesWeights %s];',variableNames{i}));
    end
    % Set up the Bayesian optimization
    fun = @(x)computeObjectiveOptim(x.Et_W,x.Tavg_ShortW,x.Tmin_ShortW,x.Tmax_ShortW,x.Pre_ShortW,x.Tavg_LongW,x.Tmin_LongW,x.Tmax_LongW,x.Pre_LongW, ...
        var,addVars,learningDates,sortedDates,saveOptimPrep,nbImages,GeoRef,GenerationType,bootstrap,ensemble,metric,optimisation,inputDir,outputDir);
    % Run the Bayesian optimization
    if parallelComputing == true
        results = bayesopt(fun,[Et_W,Tavg_ShortW,Tmin_ShortW,Tmax_ShortW,Pre_ShortW,Tavg_LongW,Tmin_LongW,Tmax_LongW,Pre_LongW], ...
            'Verbose',1,'AcquisitionFunctionName','expected-improvement-plus','MaxObjectiveEvaluations',nbOptiRuns,'UseParallel',true);
    else
        results = bayesopt(fun,[Et_W,Tavg_ShortW,Tmin_ShortW,Tmax_ShortW,Pre_ShortW,Tavg_LongW,Tmin_LongW,Tmax_LongW,Pre_LongW], ...
            'Verbose',1,'AcquisitionFunctionName','expected-improvement-plus','MaxObjectiveEvaluations',nbOptiRuns);
    end
    % Retrieve the optimal weights
    disp('  Saving optimisedWeights.mat...')
    optimisedWeights = results.XAtMinObjective;
    save(fullfile(inputDir,'optimisedWeights.mat'), 'optimisedWeights', '-v7.3','-nocompression');
    
    disp('--- 4. OPTIMISATION DONE ---')
end

tEnd = toc(tStart);

disp(['--- FINISHED IN ' num2str(tEnd) ' SECONDS ---'])
