%function objective = computeObjectiveOptim(bayesWeights,var,addVars,learningDates,shortWindow,nbImages,GeoRef,GenerationType,metric,inputDir,outputDir)
function objective = computeObjectiveOptim(Et_W,Tavg_ShortW,Tmin_ShortW,Tmax_ShortW,Pre_ShortW,Tavg_LongW,Tmin_LongW,Tmax_LongW,Pre_LongW,var,addVars,learningDates,shortWindow,nbImages,GeoRef,GenerationType,metric,optimisation,inputDir,outputDir)

%
%
%
% REDO DOCUMENTATION
%
%
%

bayesWeights = [Et_W,Tavg_ShortW,Tmin_ShortW,Tmax_ShortW,Pre_ShortW,Tavg_LongW,Tmin_LongW,Tmax_LongW,Pre_LongW];

% ... Code to execute the relevant parts ...
%sortedDates = KNNSortingOptim(var,addVars,shortWindow,bayesWeights,nbImages,inputDir);
sortedDates = KNNSortingOptim(var,addVars,shortWindow,Et_W,Tavg_ShortW,Tmin_ShortW,Tmax_ShortW,Pre_ShortW,Tavg_LongW,Tmin_LongW,Tmax_LongW,Pre_LongW,nbImages,inputDir);
OutputType  = 1;
GenerateSynImages(var,learningDates,sortedDates,GeoRef,outputDir,GenerationType,OutputType);

% Compute the validation metric using the updated code
objective = validationMetrics(metric,optimisation,outputDir);

end
