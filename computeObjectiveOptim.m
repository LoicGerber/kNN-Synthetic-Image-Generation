function objective = computeObjectiveOptim(Et_W,Tavg_ShortW,Tmin_ShortW,Tmax_ShortW,Pre_ShortW,Tavg_LongW,Tmin_LongW,Tmax_LongW,Pre_LongW, ...
    targetVar,addVars,learningDates,sortedDates,refValidation,saveOptimPrep,nbImages,geoRef,generationType,bootstrap,ensemble,metricV,validation, ...
    optimisation,inputDir,outputDir)

%
%
%
% REDO DOCUMENTATION
%
%
%

totWeights  = Et_W + Tavg_ShortW + Tmin_ShortW + Tmax_ShortW  +Pre_ShortW + Tavg_LongW + Tmin_LongW + Tmax_LongW + Pre_LongW;
Et_W        = Et_W/totWeights;
Tavg_ShortW = Tavg_ShortW/totWeights;
Tmin_ShortW = Tmin_ShortW/totWeights;
Tmax_ShortW = Tmax_ShortW/totWeights;
Pre_ShortW  = Pre_ShortW/totWeights;
Tavg_LongW  = Tavg_LongW/totWeights;
Tmin_LongW  = Tmin_LongW/totWeights;
Tmax_LongW  = Tmax_LongW/totWeights;
Pre_LongW   = Pre_LongW/totWeights;

%sortedDates = KNNSortingOptim(var,addVars,shortWindow,bayesWeights,nbImages,inputDir);
sortedDatesOptim = kNNSortingOptim(sortedDates,addVars,Et_W,Tavg_ShortW,Tmin_ShortW,Tmax_ShortW,Pre_ShortW,Tavg_LongW,Tmin_LongW,Tmax_LongW,Pre_LongW,nbImages,saveOptimPrep,inputDir);
synImages = generateSynImages(targetVar,learningDates,sortedDatesOptim,geoRef,outputDir,generationType,validation,optimisation,bootstrap,false,ensemble,2);
% Compute the validation metric using the updated code
objective = validationMetrics(targetVar,metricV,optimisation,refValidation,synImages,bootstrap,ensemble,outputDir);
if metricV == 1
    % Calculate the RMSE
    disp(['RMSE: ' num2str(objective)])
elseif metricV == 2
    % Calculate the SPEM
    disp(['SPEM: ' num2str(objective)])
elseif metricV == 3
    % Calculate the SPAEF
    disp(['SPAEF: ' num2str(objective)])
elseif metricV == 4
    % Calculate the SPOMF absolute error
    disp(['SPOMF: ' num2str(objective)])
end
