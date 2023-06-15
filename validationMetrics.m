function validationMetric = validationMetrics(metric,optimisation,refValidation,synImages,outputDir)

%
%
%
% REDO DOCUMENTATION
%
%
%

synValidation = synImages;

refImages = refValidation.maps;
refDates  = refValidation.date;
synImages = synValidation.maps;
synDates  = synValidation.date;

if size(refImages,3) ~= size(synImages,3)
    error('Numbers of reference and synthetic images do not match');
end

% Initialize an array to store the RMSE values
validationMetric = zeros(size(refImages,3), 2);

% Loop through each file in the reference directory
for i = 1:size(refImages,3)
    % Get the reference image filename and full path
    refImageDate = refDates(i);
    % Get the generated image filename and full path
    synImageDate = synDates(i);

    % Load the reference and generated images
    refImage = refImages(:,:,i);
    refImage(isnan(refImage)) = -999;
    synImage = synImages(:,:,i);
    synImage(isnan(synImage)) = -999;

    %currentDate = datetime(strrep(refImageDate,'.tif',''),'InputFormat','uuuuMMdd');
    if refImageDate == synImageDate
        currentDate = refImageDate;
    else
        error('Dates do not match')
    end

    if metric == 1
        % Calculate the RMSE
        validationMetric(i,1) = currentDate;
        validationMetric(i,2) = sqrt(immse(synImage,refImage));
        % Display the RMSE for this pair of images
        %fprintf('RMSE for %s: %.4f\n', currentDate, validationMetric(i,2));
    elseif metric == 2
        % Calculate the SPEM
        validationMetric(i,1) = currentDate;
        validationMetric(i,2) = spem(synImage,refImage);
        % Display the SPEM for this pair of images
        %fprintf('SPEM for %s: %.4f\n', currentDate, validationMetric(i,2));
    elseif metric == 3
        % Calculate the SPAEF
        validationMetric(i,1) = currentDate;
        validationMetric(i,2) = spaef(synImage,refImage);
        % Display the SPAEF for this pair of images
        %fprintf('SPAEF for %s: %.4f\n', currentDate, validationMetric(i,2));
    elseif metric == 4
        % Calculate the SPOMF absolute error
        validationMetric(i,1) = currentDate;
        validationMetric(i,2) = spae_metric(synImage,refImage);
        % Display the SPOMF absolute error for this pair of images
        %fprintf('SPOMF absolute error for %s: %.4f\n', currentDate, validationMetric(i,2));
    else
        error('Invalid metric flag...')
    end
end

if optimisation == false
    disp('Saving validationMetric.mat table...')
    validationSave = fullfile(outputDir,'validationMetric.mat');
    save(validationSave, 'validationMetric');
else % for optimisation run
    validationMetric = mean(validationMetric(:,2));
end

end
