function objective = computeObjectiveOptim(Et_W,Tavg_ShortW,Tmin_ShortW,Tmax_ShortW,Pre_ShortW,Tavg_LongW,Tmin_LongW,Tmax_LongW,Pre_LongW,var,addVars,learningDates,sortedDates,refValidation,saveOptimPrep,nbImages,GeoRef,GenerationType,bootstrap,ensemble,metric,validation,optimisation,inputDir,outputDir)

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
sortedDatesOptim = KNNSortingOptim(sortedDates,addVars,Et_W,Tavg_ShortW,Tmin_ShortW,Tmax_ShortW,Pre_ShortW,Tavg_LongW,Tmin_LongW,Tmax_LongW,Pre_LongW,nbImages,saveOptimPrep,inputDir);
OutputType  = 1;
synImages = GenerateSynImages(var,learningDates,sortedDatesOptim,GeoRef,outputDir,GenerationType,validation,optimisation,bootstrap,ensemble,OutputType);
% Compute the validation metric using the updated code
objective = validationMetrics(metric,optimisation,refValidation,synImages,outputDir);
if metric == 1
    % Calculate the RMSE
    disp(['RMSE: ' num2str(objective)])
elseif metric == 2
    % Calculate the SPEM
    disp(['SPEM: ' num2str(objective)])
elseif metric == 3
    % Calculate the SPAEF
    disp(['SPAEF: ' num2str(objective)])
elseif metric == 4
    % Calculate the SPOMF absolute error
    disp(['SPOMF: ' num2str(objective)])
end
