%% SYNTHETIC IMAGE GENERATION USING A KNN ALGORITHM

% Created 04.2023 by Loïc Gerber

% PARAMETERS DEFIINITION
% --- IMPORTANT ---
% Input files and their variable must be named as the variables defined in the variable lists below (Pre.nc -> 'Pre', Tmax.nc -> 'Tmax', etc.)
% --- IMPORTANT ---

% Setup
clc, clear

% All directories
rawDir    = 'path/to/raw/data';         % Path to raw data
inputDir  = 'path/to/saved/inputs';     % Path to saved input data, will be automatically created
outputDir = 'path/to/outputs';          % Path to results, will be automatically created

% ConvertStructureToInputs
var               = ["Et"];                         % Variables to be generated, with ["example1","example2"]
vars              = ["Tavg","Tmin","Tmax","Pre"];   % Input variables considered for the data generation, with ["example1","example2"]
addVars           = [];                             % Additional input variables, with ["example1","example2"]
QdateStart        = 19500101;                       % YYYYMMDD - Start of the Generation period
QdateEnd          = 19791231;                       % YYYYMMDD - End of the Generation period
LdateStart        = 20000101;                       % YYYYMMDD - Start of the Learning period
LdateEnd          = 20201231;                       % YYYYMMDD - End of the Learning period
outputTime        = 1;                              % Image generation timestep: 1 = DAILY, 2 = MONTHLY
precision         = 1;                              % Precision needed, 1 = single, 2 = double

% KNNDataGeneration
shortWindow       = 5;          % number of days to consider for the short climate window
longWindow        = 30;         % number of days to consider for the long climate window
nbImages          = 10;         % K, number of days to consider for the generation of images

% GenerateSynImages
ensemble          = 10;         % when using bootstrap, number of ensembles created
generationType    = 2;          % generation type,  1 = BINARY,  2 = MEAN OF SELECTED IMAGES, 3 = MEDIAN OF SELECTED IMAGES
outputType        = 2;          % output data file type, 1 = GeoTIFF, 2 = NetCDF
coordRefSysCode   = 4326;       % Coordinate reference system code, WGS84 = 4326, https://epsg.org/home.html

% Functions switches
parallelComputing = false;      % true = parallel computing ON,  false = parallel computing OFF
netCDFtoInputs    = true;       % true = create inputs,          false = load inputs
createGenWeights  = true;       % true = create generic weights, false = load optimised weights
kNNsorting        = true;       % true = create sorted data,     false = load sorted data
generateImage     = true;       % true = image generation ON,    false = image generation OFF
bootstrap         = false;      % true = bootstrap ON,           false = bootstrap OFF

% Validation switch
validationPrep    = false;      % true = validation preparation ON,    false = validation preparation OFF (!!! BYPASSES PREVIOUS SWITCHES !!!)
validation        = false;      % true = validation ON,    false = validation OFF (!!! BYPASSES PREVIOUS SWITCHES !!!)
metricViz         = false;      % true = visualisation ON, false = visualisation OFF
metric            = 1;          % 1 = RMSE, 2 = SPEM, 3 = SPAEF, 4 = Symmetric Phase-only Matched Filter-based Absolute Error Function (SPOMF)

% Bayesian optimisation switch
optimPrep         = false;      % true = optimisation preparation ON, false = optimisation preparation OFF (!!! BYPASSES PREVIOUS SWITCHES !!!)
saveOptimPrep     = false;      % true = save optimisation inputs ON, false = save optimisation inputs OFF
optimisation      = false;      % true = optimisation ON, false = optimisation OFF (!!! run AFTER optimisation preparation !!!)
nbOptiRuns        = 50;         % Number of runs for the Bayesian optimisation

% Pass all arguments to MAIN function
[geoRef,climateData,queryDates,learningDates,refValidation,additionalVars,Weights,sortedDates,synImages,validationMetric,optimisedWeights] = MAIN(...
    rawDir,inputDir,outputDir,var,vars,addVars,QdateStart,QdateEnd,LdateStart,LdateEnd,outputTime,precision, ...
    shortWindow,longWindow,nbImages,ensemble,generationType,outputType,coordRefSysCode,parallelComputing, ...
    netCDFtoInputs,createGenWeights,kNNsorting,generateImage,bootstrap,validationPrep,validation,metricViz, ...
    metric,optimPrep,saveOptimPrep,optimisation,nbOptiRuns);

