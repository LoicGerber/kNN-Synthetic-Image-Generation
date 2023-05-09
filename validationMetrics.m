function validationMetric = validationMetrics(metric,outputDir)

%
%
%
% REDO DOCUMENTATION
%
%
%

tic

% Define the two directories
refDir = fullfile(outputDir,'referenceImages');
synDir = fullfile(outputDir,'syntheticImages');

% Get a list of all files in the reference directory
refImages = dir(fullfile(refDir, '*.tif'));
% Check that there are no extra generated images
synImages = dir(fullfile(synDir, '*.tif'));
if numel(refImages) ~= numel(synImages)
    error('Number of reference images and synthetic images does not match');
end

% Initialize an array to store the RMSE values
validationMetric = zeros(numel(refImages), 2);

% Loop through each file in the reference directory
for i = 1:numel(refImages)
    % Get the reference image filename and full path
    refImageDate = refImages(i).name;
    refImagePath = fullfile(refDir, refImageDate);

    % Get the generated image filename and full path
    synImageDate = refImages(i).name;
    synImagePath = fullfile(synDir, synImageDate);

    % Check that the generated image file exists and has the correct name
    if ~exist(synImagePath, 'file')
        error('Synthetic image %s does not exist or has the wrong name', synImageDate);
    end

    % Load the reference and generated images
    refImage = imread(refImagePath);
    refImage(isnan(refImage)) = -999;
    synImage = imread(synImagePath);
    synImage(isnan(synImage)) = -999;
    
    %currentDate = datetime(strrep(refImageDate,'.tif',''),'InputFormat','uuuuMMdd');
    currentDate = convertCharsToStrings(strrep(refImageDate,'.tif',''));

    if metric == 0
        % Calculate the RMSE
        validationMetric(i,1) = currentDate;
        validationMetric(i,2) = sqrt(immse(synImage,refImage));
        % Display the RMSE for this pair of images
        %fprintf('RMSE for %s: %.4f\n', currentDate, validationMetric(i,2));
    elseif metric == 1
        % Calculate the SPEM
        validationMetric(i,1) = currentDate;
        validationMetric(i,2) = spem(synImage,refImage);
        % Display the SPEM for this pair of images
        %fprintf('SPEM for %s: %.4f\n', currentDate, validationMetric(i,2));
    elseif metric == 2
        % Calculate the SPAEF
        validationMetric(i,1) = currentDate;
        validationMetric(i,2) = spaef(synImage,refImage);
        % Display the SPAEF for this pair of images
        %fprintf('SPAEF for %s: %.4f\n', currentDate, validationMetric(i,2));
    elseif metric == 3
        % Calculate the SPOMF absolute error
        validationMetric(i,1) = currentDate;
        validationMetric(i,2) = spae_metric(synImage,refImage);
        % Display the SPOMF absolute error for this pair of images
        %fprintf('SPOMF absolute error for %s: %.4f\n', currentDate, validationMetric(i,2));
    else
        error('Invalid metric flag...')
    end
end

% Calculate the average metric
avgMetricValues = mean(validationMetric(:,2));
% Display the average metric
fprintf('Average: %.4f\n', avgMetricValues);

disp('Saving validation metric table...')
validationSave = fullfile(outputDir,'validationMetric.mat');
save(validationSave, 'validationMetric');
disp('validationMetric table saved')

toc

end