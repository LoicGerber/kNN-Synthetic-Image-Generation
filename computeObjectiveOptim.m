function objective = computeObjectiveOptim(Pre_ShortW,Tmin_ShortW,Tmax_ShortW,Pre_LongW,Tmin_LongW,Tmax_LongW,doyW,spemW, ...
    targetVar,targetDim,learningDates,sortedDates,refValidation,saveOptimPrep,metricKNN,nbImages,generationType,metricV, ...
    optimisation,useDOY,inDir,outputDir)

%
%
%
% REDO DOCUMENTATION
%
%
%

% Tavg_ShortW = 1;
% Tmin_ShortW = 1;
% Tmax_ShortW = 1;
% Pre_ShortW  = 1;
% Tavg_LongW  = 1;
% Tmin_LongW  = 1;
% Tmax_LongW  = 1;
% Pre_LongW   = 1;
% doyW        = 1;
% spemW       = 1;
% helW        = 1;

totWeights  = Pre_ShortW + Tmin_ShortW + Tmax_ShortW + Pre_LongW + Tmin_LongW + Tmax_LongW + doyW;
% Tavg_ShortW = Tavg_ShortW/totWeights;
Tmin_ShortW = Tmin_ShortW/totWeights;
Tmax_ShortW = Tmax_ShortW/totWeights;
Pre_ShortW  = Pre_ShortW/totWeights;
% Tavg_LongW  = Tavg_LongW/totWeights;
Tmin_LongW  = Tmin_LongW/totWeights;
Tmax_LongW  = Tmax_LongW/totWeights;
Pre_LongW   = Pre_LongW/totWeights;
doyW        = doyW/totWeights;

helW        = 1 - spemW;

%sortedDates = KNNSortingOptim(var,addVars,shortWindow,bayesWeights,nbImages,inputDir);
sortedDatesOptim = kNNSortingOptim(sortedDates,Tmin_ShortW,Tmax_ShortW,Pre_ShortW,Tmin_LongW,Tmax_LongW,Pre_LongW,doyW,spemW,helW,metricKNN,nbImages,useDOY,saveOptimPrep,inDir);
synImagesOptim = generateSynImgsOptim(targetVar,targetDim,learningDates,sortedDatesOptim,nbImages,generationType);
% Compute the validation metric using the updated code
objective = validationMetrics(targetVar,targetDim,metricV,optimisation,refValidation,synImagesOptim,false,[],outputDir);
if metricV == 1
    disp(['    MAE: ' num2str(objective)])
elseif metricV == 2
    disp(['    RMSE: ' num2str(objective)])
elseif metricV == 3
    disp(['    SPEM: ' num2str(objective)])
elseif metricV == 4
    disp(['    SPAEF: ' num2str(objective)])
elseif metricV == 5
    disp(['    KGE: ' num2str(objective)])
elseif metricV == 6
    disp(['    NSE: ' num2str(objective)])
end
