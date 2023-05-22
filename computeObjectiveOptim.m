function objective = computeObjectiveOptim(Weights,var,addVars,learningDates,shortWindow,nbImages,GeoRef,GenerationType,metric,inputDir,outputDir)

%
%
%
% REDO DOCUMENTATION
%
%
%

% ... Code to execute the relevant parts ...
sortedDates = KNNSortingOptim(var,addVars,shortWindow,Weights,nbImages,inputDir);
OutputType = 1;
GenerateSynImages(var,learningDates,sortedDates,GeoRef,outputDir,GenerationType,OutputType);

% Compute the validation metric using the updated code
objective = validationMetrics(metric, outputDir);

end
