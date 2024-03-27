%% SYNTHETIC IMAGE GENERATION USING A KNN ALGORITHM

% Created 04.2023 by Loïc Gerber

% PARAMETERS DEFIINITION
% --- IMPORTANT ---
% Input files and their variable must be named as the variables defined in the variable lists below (Pre.nc -> 'Pre', Tmax.nc -> 'Tmax', etc.)
% --- IMPORTANT ---

% Setup
clc, clear

% All directories
rawDir            = 'path/to/raw/data';                         % Path to raw data
inputDir          = 'path/to/saved/inputs';                     % Path to saved input data, will be automatically created
outputDir         = 'path/to/outputs';                          % Path to results, will be automatically created
optiWeightsDir    = 'path/to/optimised/weights/matrix.mat';     % Path to optimised weights matrix .mat file (if available, otherwise '')
maskDir           = 'path/to/binary/mask';                      % Path to binary mask delimiting the area to generate with the pixel-based approach (.tif image)

% ConvertStructureToInputs
targetVar         = ["Et"];                         % Variables to be generated, with ["example1","example2"]
climateVars       = ["Tavg","Tmin","Tmax","Pre"];   % Input variables considered for the data generation, with ["example1","example2"]
addVars           = ["sm","twsa"];                  % Additional input variables, with ["example1","example2"], if empty use []
normMethods       = [1,1,1,4];                      % Method of normalisation of the climate data. 1 = MinMax, 2 = Q10-Q90, 3 = log(x+0.01), 4 = log(x>0)
maxThreshold      = 30;                             % Days, max threshold for closest additional variables attribution
QdateStart        = 19500101;                       % YYYYMMDD - Start of the Generation period
QdateEnd          = 19791231;                       % YYYYMMDD - End of the Generation period
LdateStart        = 19800101;                       % YYYYMMDD - Start of the Learning period
LdateEnd          = 20201231;                       % YYYYMMDD - End of the Learning period
outputTime        = 1;                              % Image generation timestep: 1 = DAILY, 2 = MONTHLY, 3 = DEKADAL

% KNNDataGeneration
pixelWise         = true;       % Pixel-wise kNN (true), or domain-wise kNN (false)
shortWindow       = 5;          % Number of days to consider for the short climate window
longWindow        = 20;         % Number of days to consider for the long climate window (total days, including shortWindow)
nbImages          = 10;         % K, number of days to consider for the generation of images
metricKNN         = 2;          % 1 = RMSE, 2 = MAE, 3 = Manhattan distance, 4 = Euclidean distance
daysRange         = 90;         % Range of possible learning days (before and after) the query date's DOY

% GenerateSynImages
ensemble          = 10;         % When using bootstrap, number of ensembles created
generationType    = 2;          % Generation type,       1 = BINARY,  2 = MEAN OF SELECTED IMAGES, 3 = MEDIAN OF SELECTED IMAGES
outputType        = 2;          % Output data file type, 1 = GeoTIFF, 2 = NetCDF
coordRefSysCode   = 4326;       % Coordinate reference system code, WGS84 = 4326, https://epsg.org/home.html

% Functions switches
parallelComputing = false;      % true = parallel computing ON,         false = parallel computing OFF
netCDFtoInputs    = true;       % true = create inputs,                 false = load inputs
createGenWeights  = true;       % true = create generic weights,        false = load optimised weights
kNNsorting        = true;       % true = create sorted data,            false = load sorted data
generateImage     = true;       % true = image generation ON,           false = image generation OFF

% Bootstrap switches (NOT INFLUENCED BY VALIDATION SWITCHES)
bootstrap         = false;      % true = bootstrap ON,                  false = bootstrap OFF
bsSaveAll         = false;      % true = saves all bootstrap ensembles, false = saves only min, deterministic and max bootstrap ensembles as netCDF files

% Validation switches
validationPrep    = false;      % true = validation preparation ON, false = validation preparation OFF (!!! BYPASSES PREVIOUS SWITCHES !!!)
validation        = false;      % true = validation ON,             false = validation OFF (!!! BYPASSES PREVIOUS SWITCHES !!!)
metricViz         = false;      % true = visualisation ON,          false = visualisation OFF
metricV           = 1;          % 1 = RMSE, 2 = SPEM, 3 = SPAEF, 4 = Symmetric Phase-only Matched Filter-based Absolute Error Function (SPOMF)

% Bayesian optimisation switches
optimPrep         = false;      % true = optimisation preparation ON, false = optimisation preparation OFF (!!! BYPASSES PREVIOUS SWITCHES !!!)
saveOptimPrep     = false;      % true = save optimisation inputs ON, false = save optimisation inputs OFF
optimisation      = false;      % true = optimisation ON,             false = optimisation OFF (!!! run AFTER optimisation preparation !!!)
nbOptiRuns        = 50;         % Number of runs for the Bayesian optimisation algorithm

% Pass all arguments to MAIN function
[geoRef,climateData,queryDates,learningDates,refValidation,additionalVars, ...
    Weights,sortedDates,synImages,validationMetric,optimisedWeights] = ...
    MAIN(...
    rawDir,inputDir,outputDir,optiWeightsDir,maskDir,targetVar,climateVars,addVars,normMethods,QdateStart,QdateEnd,LdateStart,LdateEnd,outputTime, ...
    maxThreshold,shortWindow,longWindow,daysRange,nbImages,metricKNN,ensemble,generationType,outputType,coordRefSysCode,parallelComputing, ...
    netCDFtoInputs,createGenWeights,kNNsorting,generateImage,bootstrap,bsSaveAll,validationPrep,validation,pixelWise, ...
    metricViz,metricV,optimPrep,saveOptimPrep,optimisation,nbOptiRuns);

