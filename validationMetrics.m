function validationMetric = validationMetrics(metric,optimisation,refValidation,synImages,bootstrap,ensemble,outputDir)

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
synDates  = synValidation.date;

if bootstrap == false
    % Initialize an array to store the RMSE values
    validationMetric = zeros(size(refImages,3), 2);
    synImagesAll = synValidation.maps;
    if size(refImages,3) ~= size(synImagesAll,3)
        error('Numbers of reference and synthetic images do not match');
    end
else
    % Initialize an array to store the RMSE values
    validationMetric = cell(size(refImages,3), 3);
    synImagesAll = synValidation.bootstrap;
    maps = synValidation.maps;
    if size(refImages,3) ~= size(synImagesAll,1)
        error('Numbers of reference and synthetic images do not match');
    end
end



% Loop through each file in the reference directory
for i = 1:size(refImages,3)
    % Get the reference image filename and full path
    refImageDate = refDates(i);
    % Get the generated image filename and full path
    synImageDate = synDates(i);

    % Load the reference and generated images
    refImage = refImages(:,:,i);
    refImage(isnan(refImage)) = -999;
    if bootstrap == false
        synImage = synImagesAll(:,:,i);
        synImage(isnan(synImage)) = -999;
        %currentDate = datetime(strrep(refImageDate,'.tif',''),'InputFormat','uuuuMMdd');
        if refImageDate == synImageDate
            currentDate = refImageDate;
        else
            error('Dates do not match')
        end
        validationMetric(i,1) = currentDate;
        if metric == 1
            % Calculate the RMSE
            validationMetric(i,2) = sqrt(immse(synImage,refImage));
        elseif metric == 2
            % Calculate the SPEM
            validationMetric(i,2) = spem(synImage,refImage);
        elseif metric == 3
            % Calculate the SPAEF
            validationMetric(i,2) = spaef(synImage,refImage);
        elseif metric == 4
            % Calculate the SPOMF absolute error
            validationMetric(i,2) = spae_metric(synImage,refImage);
        else
            error('Invalid metric flag...')
        end
    else
        %currentDate = datetime(strrep(refImageDate,'.tif',''),'InputFormat','uuuuMMdd');
        if refImageDate == synImageDate
            currentDate = refImageDate;
        else
            error('Dates do not match')
        end
        validationMetric{i,1} = currentDate;
        % Bootstrap ensembles
        for j = 1:ensemble
            synImage = single(synImagesAll{i}(:,:,j));
            synImage(isnan(synImage)) = -999;
            if metric == 1
                % Calculate the RMSE
                validationMetric{i,2}(j) = sqrt(immse(synImage,refImage));
            elseif metric == 2
                % Calculate the SPEM
                validationMetric{i,2}(j) = spem(synImage,refImage);
            elseif metric == 3
                % Calculate the SPAEF
                validationMetric{i,2}(j) = spaef(synImage,refImage);
            elseif metric == 4
                % Calculate the SPOMF absolute error
                validationMetric{i,2}(j) = spae_metric(synImage,refImage);
            else
                error('Invalid metric flag...')
            end
        end
        % KNN result
        synImage = single(maps(:,:,i));
        synImage(isnan(synImage)) = -999;
        if metric == 1
            % Calculate the RMSE
            validationMetric{i,3} = sqrt(immse(synImage,refImage));
        elseif metric == 2
            % Calculate the SPEM
            validationMetric{i,3} = spem(synImage,refImage);
        elseif metric == 3
            % Calculate the SPAEF
            validationMetric{i,3} = spaef(synImage,refImage);
        elseif metric == 4
            % Calculate the SPOMF absolute error
            validationMetric{i,3} = spae_metric(synImage,refImage);
        else
            error('Invalid metric flag...')
        end
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
