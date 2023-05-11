%% PARAMETERS DEFIINITION

%
%
%
% REDO DOCUMENTATION
%
%
%

% --- IMPORTANT ---
% Input files and their variable must be named as the variables defined below (pre.nc, 'pre', tmax.nc, 'tmax', etc.)
% --- IMPORTANT ---

% Setup
clc, clear, close all
poolobj = gcp('nocreate');
delete(poolobj);
tStart = tic;

% All directories
rawDir    = 'X://LoicGerber\knn_image_generation\syntheticImageGeneration\voltaData\voltaClean\';
inputDir  = 'X://LoicGerber\knn_image_generation\syntheticImageGeneration\voltaResults1979\inputData\';
outputDir = 'X://LoicGerber\knn_image_generation\syntheticImageGeneration\voltaResults1979\output\';

% ConvertStructureToInputs
var             = "Et";                          % Variable to be generated, with "example"
vars            = ["Tavg","Tmin","Tmax","Pre"];  % Input variables considered for the data generation, with ["example1","example2"]
addVars         = [];                            % Additional input variables, with ["example1","example2"]
QdateStart      = 19790101;                      % YYYYMMDD - Start of the Generation period
QdateEnd        = 19791231;                      % YYYYMMDD - End of the Generation period
LdateStart      = 19800101;                      % YYYYMMDD - Start of the Learning period
LdateEnd        = 20201212;                      % YYYYMMDD - End of the Learning period
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
KNNsorting      = 1;    % 0 = create sorted data,     1 = load sorted data
generateImage   = 0;    % 0 = IMAGE GENERATION OFF,   1 = IMAGE GENERATION ON

% Validation switch
validation      = 0;    % 0 = VALIDATION OFF, 1 = VALIDATION ON (!!! BYPASSES PREVIOUS SWITCHES !!!)
metricEval      = 1;    % 0 = metrics evaluation off, 1 = metrics evaluation on
metric          = 0;    % 0 = RMSE, 1 = SPEM, 2 = SPAEF, 3 = Symmetric Phase-only Matched Filter-based Absolute Error Function (SPOMF)

%% Reading the data needed for ranking learning dates using "KNNDataSorting" Function
disp('--- 1. READING DATA ---')

if NetCDFtoInputs == 0 && validation == 0
    disp('Formatting input data for production run...')
    rawData = ConvertNetCDFtoStructure(rawDir,inputDir);
    climateData = extractClimateData(vars,rawData,QdateStart,QdateEnd,LdateStart,LdateEnd,longWindow,validation,inputDir);
    learningDates  = ConvertStructureToLearningDates(var,LdateStart,LdateEnd,rawData,climateData,inputDir);
    [queryDates,learningDates] = ConvertStructureToQueryDates(var,QdateStart,QdateEnd,learningDates,climateData,longWindow,validation,outputTime,inputDir,outputDir);
    additionalVars = extractAdditionalVars(addVars,rawData,QdateStart,QdateEnd,LdateStart,LdateEnd,inputDir);
    R = extractGeoInfo(var,coordRefSysCode,rawDir,inputDir);
    R = [];
elseif validation == 1
    disp('Formatting input data for validation run...')
    rawData = ConvertNetCDFtoStructure(rawDir,inputDir);
    climateData = extractClimateData(vars,rawData,QdateStart,QdateEnd,LdateStart,LdateEnd,longWindow,validation,inputDir);
    learningDates = ConvertStructureToLearningDates(var,LdateStart,LdateEnd,rawData,climateData,inputDir);
    [queryDates,learningDates] = ConvertStructureToQueryDates(var,QdateStart,QdateEnd,learningDates,climateData,longWindow,validation,outputTime,inputDir,outputDir);
    additionalVars = extractAdditionalVars(addVars,rawData,QdateStart,QdateEnd,LdateStart,LdateEnd,inputDir);
    R = extractGeoInfo(var,coordRefSysCode,rawDir,inputDir);
    R = [];
elseif NetCDFtoInputs == 1 && validation == 0
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
    disp('Loading R.mat file...')
    R              = load(fullfile(inputDir,'R.mat'));
    R              = R.R;
    R              = [];
end
disp('Loading Weights.mat file...')
Weights       = createWeights(var,vars,addVars,inputDir);
% Optional - calibrated weights
%Weights       = load('BayesResultXAtMinObjective.mat');
%Weights       = Weights.BayesResultXAtMinObjective;

disp('--- 1. READING DATA DONE ---')

%% The function for Ranking the learning dates
disp('--- 2. KNN DATA SORTING ---')

% Generate ranked Learning Dates for each Query Date
if KNNsorting == 0 || validation == 1
    sortedDates = KNNDataSorting(var,vars,addVars,queryDates,learningDates,climateData,additionalVars,shortWindow,longWindow,Weights,nbImages,inputDir);
elseif KNNsorting == 1 && validation == 0
    disp('Loading necessary datasets...')
    disp('KNNsorting flag == 1, load pre-existing data...')
    sortedDates = load(fullfile(inputDir,'KNNSorting.mat'));
    sortedDates = sortedDates.sortedDates;
end

disp('--- 2. KNN DATA SORTING DONE ---')

%% Generation of Synthetic Images
disp('--- 3. SYNTHETIC IMAGES GENERATION ---')

if generateImage == 1 && validation == 0
    GenerateSynImages(var,learningDates,sortedDates,R,outputDir,GenerationType,OutputType);
elseif validation == 1
    OutputType = 1;
    GenerateSynImages(var,learningDates,sortedDates,R,outputDir,GenerationType,OutputType);
elseif generateImage == 0 && validation == 0
    disp('generateImage flag == 0, no synthetic images generated...')
end

disp('--- 3. SYNTHETIC IMAGES GENERATION DONE ---')

%% Validation (optional)
if metricEval == 1 && validation == 1
    disp('--- 4. VALIDATION ---')

    validationMetric = validationMetrics(metric,outputDir);
    visualiseMetrics(validationMetric,metric,LdateStart,LdateEnd,outputDir);

    disp('--- 4. VALIDATION DONE ---')
else
    disp('--- VALIDATION SKIPPED ---')
end

tEnd = toc(tStart);

disp(['--- FINISHED IN ' num2str(tEnd) ' SECONDS ---'])
