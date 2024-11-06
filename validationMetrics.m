function validationMetric = validationMetrics(targetVar,targetDim,nanValue,metricV,optimisation,refValidation,synImages,bootstrap,ensemble,outputDir)

%
%
%
% REDO DOCUMENTATION
%
%
%

synValidation = synImages;

validOptim = 0;

targetVarL = lower(targetVar);

for j = 1:numel(targetVar)
    refImages = refValidation.(targetVarL(j));
    refDates  = refValidation.date;
    synDates  = synValidation.date;

    if bootstrap == false
        % Initialize an array to store the RMSE values
        if targetDim ~= 1
            validationResult = zeros(size(refImages,3), 2);
            synImagesAll = synValidation.(targetVarL(j));
            if size(refImages,3) ~= size(synImagesAll,3)
                error('Numbers of reference and synthetic images do not match');
            end
        else
            validationResult = zeros(size(refImages,1),1);
            synImagesAll = synValidation.(targetVarL(j));
            if size(validationResult) ~= size(synImagesAll)
                error('Numbers of reference and synthetic images do not match');
            end
        end
    else
        % Initialize an array to store the RMSE values
        validationResult = cell(size(refImages,3), 3);
        varBS = strcat(targetVar(j), "_Bootstrap");
        synImagesAll = synValidation.(varBS);
        maps = synValidation.(targetVarL(j));
        if size(refImages,3) ~= size(synImagesAll,1)
            error('Numbers of reference and synthetic images do not match');
        end
    end
    
    % Initialize accumulation variables
    accumulatedRef = [];
    accumulatedSyn = [];

    % Loop through each file in the reference directory
    for i = 1:size(validationResult,1)
        % Get the reference image filename and full path
        refImageDate = refDates(i);
        % Get the generated image filename and full path
        synImageDate = synDates(i);

        % Load the reference and generated images
        if targetDim ~= 1
            refImage = double(refImages(:,:,i));
        else
            refImage = double(refImages(i));
        end
        refImage(isnan(refImage)) = -999;
        if bootstrap == false
            if targetDim ~= 1
                synImage = double(synImagesAll(:,:,i));
            else
                synImage = double(synImagesAll(i));
            end
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
            elseif metricV == 4  || metricV == 5
                if targetDim == 2
                    validationResult(i,2) = computeKGE(synImage,refImage,nanValue);
                else
                    % Accumulate data for KGE calculation
                    accumulatedRef = [accumulatedRef; refImage(:)];
                    accumulatedSyn = [accumulatedSyn; synImage(:)];
                end
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
                synImage = double(synImagesAll{i}(:,:,k));
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
                elseif metricV == 4 || metricV == 5
                    if targetDim == 2
                        validationResult{i,2}(k) = computeKGE(synImage,refImage,nanValue);
                    else
                        % Accumulate data for KGE calculation
                        accumulatedRef = [accumulatedRef; refImage(:)];
                        accumulatedSyn = [accumulatedSyn; synImage(:)];
                    end
                else
                    error('Invalid metric flag...')
                end
            end
            % KNN result
            synImage = double(maps(:,:,i));
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
            else
                error('Invalid metric flag...')
            end
        end
    end
    
    if metricV == 4 && targetDim == 1 % Compute KGE for the entire time series
        kgeValue = computeKGE(accumulatedSyn, accumulatedRef, nanValue);
        if bootstrap == false
            validationResult(:,2) = kgeValue;
        else
            validationResult{:,2} = kgeValue;
        end
    elseif metricV == 5 && targetDim == 1 % Compute NSE for the entire time series
        nseValue = computeNSE(accumulatedSyn, accumulatedRef);
        if bootstrap == false
            validationResult(:,2) = nseValue;
        else
            validationResult{:,2} = nseValue;
        end
    end

    if optimisation == false
        validationMetric.(targetVarL(j)) = validationResult;
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
