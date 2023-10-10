function validationMetric = validationMetrics(targetVar,metricV,optimisation,refValidation,synImages,bootstrap,ensemble,outputDir)

%
%
%
% REDO DOCUMENTATION
%
%
%

synValidation = synImages;

validOptim = 0;

for j = 1:numel(targetVar)
    refImages = refValidation.(targetVar(j));
    refDates  = refValidation.date;
    synDates  = synValidation.date;

    if bootstrap == false
        % Initialize an array to store the RMSE values
        validationResult = zeros(size(refImages,3), 2);
        synImagesAll = synValidation.(targetVar(j));
        if size(refImages,3) ~= size(synImagesAll,3)
            error('Numbers of reference and synthetic images do not match');
        end
    else
        % Initialize an array to store the RMSE values
        validationResult = cell(size(refImages,3), 3);
        varBS = strcat(targetVar(j), "_Bootstrap");
        synImagesAll = synValidation.(varBS);
        maps = synValidation.(targetVar(j));
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
            validationResult(i,1) = currentDate;
            if metricV == 1
                % Calculate the RMSE
                validationResult(i,2) = sqrt(immse(synImage,refImage));
            elseif metricV == 2
                % Calculate the SPEM
                validationResult(i,2) = spem(synImage,refImage);
            elseif metricV == 3
                % Calculate the SPAEF
                validationResult(i,2) = spaef(synImage,refImage);
            elseif metricV == 4
                % Calculate the SPOMF absolute error
                validationResult(i,2) = spae_metric(synImage,refImage);
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
            validationResult{i,1} = currentDate;
            % Bootstrap ensembles
            for k = 1:ensemble
                synImage = single(synImagesAll{i}(:,:,k));
                synImage(isnan(synImage)) = -999;
                if metricV == 1
                    % Calculate the RMSE
                    validationResult{i,2}(k) = sqrt(immse(synImage,refImage));
                elseif metricV == 2
                    % Calculate the SPEM
                    validationResult{i,2}(k) = spem(synImage,refImage);
                elseif metricV == 3
                    % Calculate the SPAEF
                    validationResult{i,2}(k) = spaef(synImage,refImage);
                elseif metricV == 4
                    % Calculate the SPOMF absolute error
                    validationResult{i,2}(k) = spae_metric(synImage,refImage);
                else
                    error('Invalid metric flag...')
                end
            end
            % KNN result
            synImage = single(maps(:,:,i));
            synImage(isnan(synImage)) = -999;
            if metricV == 1
                % Calculate the RMSE
                validationResult{i,3} = sqrt(immse(synImage,refImage));
            elseif metricV == 2
                % Calculate the SPEM
                validationResult{i,3} = spem(synImage,refImage);
            elseif metricV == 3
                % Calculate the SPAEF
                validationResult{i,3} = spaef(synImage,refImage);
            elseif metricV == 4
                % Calculate the SPOMF absolute error
                validationResult{i,3} = spae_metric(synImage,refImage);
            else
                error('Invalid metric flag...')
            end
        end
    end
    if optimisation == false
        validationMetric.(targetVar(j)) = validationResult;
    else
        validOptim = validOptim + mean(validationResult(:,2));
    end
end

if optimisation == false
    disp('Saving validationMetric.mat table...')
    validationSave = fullfile(outputDir,'validationMetric.mat');
    save(validationSave, 'validationMetric');
else % for optimisation run
    validationMetric = validOptim;
end

end
