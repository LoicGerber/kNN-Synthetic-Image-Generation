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
rawDir    = 'X://LoicGerber\knn_image_generation\syntheticImageGeneration\voltaData\voltaClean\';           % Path to raw data
inputDir  = 'X://LoicGerber\knn_image_generation\syntheticImageGeneration\test\inputData\';     % Path to saved input data
outputDir = 'X://LoicGerber\knn_image_generation\syntheticImageGeneration\test\output\';        % Path to results

% ConvertStructureToInputs
var             = "Et";                          % Variable to be generated, with "example"
vars            = ["Tavg","Tmin","Tmax","Pre"];  % Input variables considered for the data generation, with ["example1","example2"]
addVars         = [];                            % Additional input variables, with ["example1","example2"]
QdateStart      = 20000101;                      % YYYYMMDD - Start of the Generation period
QdateEnd        = 20000105;                      % YYYYMMDD - End of the Generation period
LdateStart      = 19800101;                      % YYYYMMDD - Start of the Learning period
LdateEnd        = 19800131;                      % YYYYMMDD - End of the Learning period
outputTime      = 1;                             % Image generation timestep: 1 = DAILY, 2 = MONTHLY
longWindow      = 30;                            % number of days to consider for the long climate window

% KNNDataGeneration
shortWindow     = 5;    % number of days to consider for the short climate window
nbImages        = 5;    % number of days to consider for the generation of images

% GenerateSynImages
GenerationType  = 2;    % data generation type,  1 = BINARY,  2 = MEAN OF SELECTED IMAGES, 3 = MEDIAN OF SELECTED IMAGES
OutputType      = 1;    % output data file type, 1 = GeoTIFF, 2 = individual NetCDF files
coordRefSysCode = 4326; % Coordinate reference system code, WGS84 = 4326, https://epsg.org/home.html

% Functions switches
NetCDFtoInputs  = 1;    % 0 = create inputs,          1 = load inputs
loadOptiWeights = 1;    % 0 = create generic weights, 1 = load weights
KNNsorting      = 0;    % 0 = create sorted data,     1 = load sorted data
generateImage   = 1;    % 0 = image generation ON,    1 = image generation OFF

% Validation switch
validation      = 1;    % 0 = validation ON,    1 = validation OFF (!!! BYPASSES PREVIOUS SWITCHES !!!)
metricViz       = 1;    % 0 = visualisation ON, 1 = visualisation OFF
metric          = 1;    % 1 = RMSE, 2 = SPEM, 3 = SPAEF, 4 = Symmetric Phase-only Matched Filter-based Absolute Error Function (SPOMF)

% Bayesian optimisation switch
optimisation    = 0;    % 0 = optimisation ON, 1 = optimisation OFF

%% Reading the data needed for ranking learning dates using "KNNDataSorting" Function
disp('--- 1. READING DATA ---')

if NetCDFtoInputs == 0 && validation == 1
    disp('Formatting input data for production run...')
    rawData                    = ConvertNetCDFtoStructure(var,vars,addVars,rawDir,inputDir);
    GeoRef                     = extractGeoInfo(var,coordRefSysCode,rawDir,inputDir);
    climateData                = extractClimateData(vars,rawData,QdateStart,QdateEnd,LdateStart,LdateEnd,longWindow,validation,inputDir);
    learningDates              = ConvertStructureToLearningDates(var,LdateStart,LdateEnd,rawData,climateData,inputDir);
    [queryDates,learningDates] = ConvertStructureToQueryDates(var,QdateStart,QdateEnd,learningDates,climateData,longWindow,GeoRef,validation,outputTime,inputDir,outputDir);
    additionalVars             = extractAdditionalVars(addVars,rawData,QdateStart,QdateEnd,LdateStart,LdateEnd,inputDir);
elseif validation == 0
    disp('Formatting input data for validation run...')
    rawData                    = ConvertNetCDFtoStructure(var,vars,addVars,rawDir,inputDir);
    GeoRef                     = extractGeoInfo(var,coordRefSysCode,rawDir,inputDir);
    climateData                = extractClimateData(vars,rawData,QdateStart,QdateEnd,LdateStart,LdateEnd,longWindow,validation,inputDir);
    learningDates              = ConvertStructureToLearningDates(var,LdateStart,LdateEnd,rawData,climateData,inputDir);
    [queryDates,learningDates] = ConvertStructureToQueryDates(var,QdateStart,QdateEnd,learningDates,climateData,longWindow,GeoRef,validation,outputTime,inputDir,outputDir);
    additionalVars             = extractAdditionalVars(addVars,rawData,QdateStart,QdateEnd,LdateStart,LdateEnd,inputDir);
elseif (NetCDFtoInputs == 1 && validation == 1) || optimisation == 0
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
if loadOptiWeights == 0
    disp('Creating generic Weights.mat file...')
    Weights = createWeights(var,vars,addVars,inputDir);
elseif loadOptiWeights == 1
    disp('Loading Weights.mat file...')
    Weights = load(fullfile(inputDir,'Weights.mat'));
    Weights = Weights.Weights;
end

disp('--- 1. READING DATA DONE ---')

%% The function for Ranking the learning dates
disp('--- 2. KNN DATA SORTING ---')

% Generate ranked Learning Dates for each Query Date
if KNNsorting == 0 || validation == 0
    sortedDates = KNNDataSorting(var,vars,addVars,queryDates,learningDates,climateData,additionalVars,shortWindow,longWindow,Weights,nbImages,optimisation,inputDir);
elseif KNNsorting == 1 && validation == 1
    disp('Loading sortedDates.mat file...')
    sortedDates = load(fullfile(inputDir,'KNNSorting.mat'));
    sortedDates = sortedDates.sortedDates;
end

disp('--- 2. KNN DATA SORTING DONE ---')

%% Generation of Synthetic Images
disp('--- 3. SYNTHETIC IMAGES GENERATION ---')

if generateImage == 0 && validation == 1 && optimisation == 1
    GenerateSynImages(var,learningDates,sortedDates,GeoRef,outputDir,GenerationType,OutputType);
elseif validation == 0
    OutputType = 1;
    GenerateSynImages(var,learningDates,sortedDates,GeoRef,outputDir,GenerationType,OutputType);
elseif optimisation == 0
    disp('Optimisation run, synthetic image generation skipped...')
elseif generateImage == 1 && validation == 1
    disp('Synthetic image generation skipped...')
end

disp('--- 3. SYNTHETIC IMAGES GENERATION DONE ---')

%% Validation
if (validation == 0 || metricViz == 0) && optimisation == 1
    disp('--- 4. VALIDATION ---')

    validationMetric = validationMetrics(metric,outputDir);
    visualiseMetrics(validationMetric,metric,LdateStart,LdateEnd,outputDir);

    disp('--- 4. VALIDATION DONE ---')
end

%% Optimisation
if optimisation == 0
    disp('--- 4. OPTIMISATION ---')

    % --- Fati ---
    % fun = @(x)-(ErrorFunction(FindImageIndex(Day,Mean,x.Et_W,x.wCloseAggreTmin,x.wCloseAggreP,x.wTmax,x.wTmin,x.wP,x.wMODIS,x.wMYD,x.wShadow,x.wClosestLandsat),Workspace2));
    % results = bayesopt(fun,[wCloseAggreTmax,wCloseAggreTmin,wCloseAggreP,wTmax,wTmin,wP,wMODIS,wMYD,wShadow,wClosestLandsat],'Verbose',1,...
    %     'AcquisitionFunctionName','expected-improvement-plus','MaxObjectiveEvaluations',1);
    % --- Fati ---
    
    
    % Get the table variable names
    variableNames = Weights.Properties.VariableNames;
    % Iterate over each variable
    for i = 1:numel(variableNames)
        bayesWeights.variableNames{i} = optimizableVariable(variableNames{i},[0,1],'Type','real');
    end
    % Set up the Bayesian optimization
    fun = @(x)computeObjectiveOptim(Weights,var,addVars,learningDates,shortWindow,nbImages,GeoRef,GenerationType,metric,inputDir,outputDir);
    % Run the Bayesian optimization
    results = bayesopt(bOptim,Weights,'Verbose',1,'AcquisitionFunctionName','expected-improvement-plus','MaxObjectiveEvaluations',1);
    % Retrieve the optimal weights
    optimalWeights = results.XAtMinObjective;
    
    
    disp('--- 4. OPTIMISATION DONE ---')
end

tEnd = toc(tStart);

disp(['--- FINISHED IN ' num2str(tEnd) ' SECONDS ---'])
