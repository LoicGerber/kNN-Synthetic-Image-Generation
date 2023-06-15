function objective = computeObjectiveOptim(Et_W,Tavg_ShortW,Tmin_ShortW,Tmax_ShortW,Pre_ShortW,Tavg_LongW,Tmin_LongW,Tmax_LongW,Pre_LongW,var,addVars,learningDates,sortedDates,refValidation,saveOptimPrep,nbImages,GeoRef,GenerationType,bootstrap,ensemble,metric,optimisation,inputDir,outputDir)

%
%
%
% REDO DOCUMENTATION
%
%
%

disp(['Et_W' num2str(Et_W)])
disp(['Tavg_ShortW' num2str(Tavg_ShortW)])
disp(['Tmin_ShortW' num2str(Tmin_ShortW)])
disp(['Tmax_ShortW' num2str(Tmax_ShortW)])
disp(['Pre_ShortW' num2str(Pre_ShortW)])
disp(['Tavg_LongW' num2str(Tavg_LongW)])
disp(['Tmin_LongW' num2str(Tmin_LongW)])
disp(['Tmax_LongW' num2str(Tmax_LongW)])
disp(['Pre_LongW' num2str(Pre_LongW)])

%sortedDates = KNNSortingOptim(var,addVars,shortWindow,bayesWeights,nbImages,inputDir);
sortedDatesOptim = KNNSortingOptim(sortedDates,addVars,Et_W,Tavg_ShortW,Tmin_ShortW,Tmax_ShortW,Pre_ShortW,Tavg_LongW,Tmin_LongW,Tmax_LongW,Pre_LongW,nbImages,saveOptimPrep,inputDir);
OutputType  = 1;
synImages = GenerateSynImages(var,learningDates,sortedDatesOptim,GeoRef,outputDir,GenerationType,optimisation,bootstrap,ensemble,OutputType);

% Compute the validation metric using the updated code
objective = validationMetrics(metric,optimisation,refValidation,synImages,outputDir);
disp(['Final result ' num2str(objective)])

end
