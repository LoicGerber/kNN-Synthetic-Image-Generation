function validationMetric = validationMetrics(targetVar,targetDim,metricV,optimisation,refValidation,synImages,stochastic,ensemble,outputDir)

%
%
%
% REDO DOCUMENTATION
%
%
%

synValidation = synImages;

targetVarL = lower(targetVar);

for j = 1:numel(targetVar)
    refImages = double(refValidation.(targetVarL(j)));
    refDates  = refValidation.date;
    synDates  = synValidation.date;

    if stochastic == false
        % Initialize an array to store the RMSE values
        if targetDim ~= 1
            validationResult = zeros(size(refImages,3), 2);
            synImagesAll = double(synValidation.(targetVarL(j)));
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
        varBS = strcat(targetVar(j), "_stochastic");
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
        %refImage(isnan(refImage)) = nanValue;
        if stochastic == false
            if targetDim ~= 1
                synImage = double(synImagesAll(:,:,i));
            else
                synImage = double(synImagesAll(i));
            end
            %synImage(isnan(synImage)) = nanValue;
            %currentDate = datetime(strrep(refImageDate,'.tif',''),'InputFormat','uuuuMMdd');
            if refImageDate == synImageDate
                currentDate = refImageDate;
            else
                error('Dates do not match')
            end
            validationResult(i,1) = currentDate;
            switch metricV
                case 1
                    % Calculate the MAE
                    validationResult(i,2) = mean(abs(synImage - refImage), 'all', 'omitnan');
                case 2
                    % Calculate the RMSE
                    validationResult(i,2) = sqrt(mean((synImage - refImage).^2, 'all', 'omitnan'));
                case 3
                    % Calculate the SPEM
                    validationResult(i,2) = spem(synImage,refImage);
                case 4
                    % Calculate the SPAEF
                    validationResult(i,2) = spaef(synImage,refImage);
                case 5 || 6
                    if targetDim == 2
                        validationResult(i,2) = computeKGE(synImage,refImage);
                    else
                        % Accumulate data for KGE calculation
                        accumulatedRef = [accumulatedRef; refImage(:)];
                        accumulatedSyn = [accumulatedSyn; synImage(:)];
                    end
            end
        else
            %currentDate = datetime(strrep(refImageDate,'.tif',''),'InputFormat','uuuuMMdd');
            if refImageDate == synImageDate
                currentDate = refImageDate;
            else
                error('Dates do not match')
            end
            validationResult{i,1} = currentDate;
            % stochastic ensembles
            for k = 1:ensemble
                synImage = double(synImagesAll{i}(:,:,k));
                %synImage(isnan(synImage)) = nanValue;
                switch metricV
                    case 1
                        % Calculate the MAE
                        validationResult{i,2}(k) = mean(abs(synImage - refImage), 'all', 'omitnan');
                    case 2
                        % Calculate the RMSE
                        validationResult{i,2}(k) = sqrt(mean((synImage - refImage).^2, 'all', 'omitnan'));
                    case 3
                        % Calculate the SPEM
                        validationResult{i,2}(k) = spem(synImage,refImage);
                    case 4
                        % Calculate the SPAEF
                        validationResult{i,2}(k) = spaef(synImage,refImage);
                    case 5 || 6
                        if targetDim == 2
                            validationResult{i,2}(k) = computeKGE(synImage,refImage);
                        else
                            % Accumulate data for KGE calculation
                            accumulatedRef = [accumulatedRef; refImage(:)];
                            accumulatedSyn = [accumulatedSyn; synImage(:)];
                        end
                end
            end
            % KNN result
            synImage = double(maps(:,:,i));
            %synImage(isnan(synImage)) = nanValue;
            switch metricV
                case 1
                    % Calculate the MAE
                    validationResult{i,3} = mean(abs(synImage - refImage), 'all', 'omitnan');
                case 2
                    % Calculate the RMSE
                    validationResult{i,3} = sqrt(mean((synImage - refImage).^2, 'all', 'omitnan'));
                case 3
                    % Calculate the SPEM
                    validationResult{i,3} = spem(synImage,refImage);
                case 4
                    % Calculate the SPAEF
                    validationResult{i,3} = spaef(synImage,refImage);
            end
        end
    end

    if metricV == 5 && targetDim == 1 % Compute KGE for the entire time series
        kgeValue = computeKGE(accumulatedSyn, accumulatedRef);
        if stochastic == false
            validationResult(:,2) = kgeValue;
        else
            validationResult{:,2} = kgeValue;
        end
    elseif metricV == 6 && targetDim == 1 % Compute NSE for the entire time series
        nseValue = computeNSE(accumulatedSyn, accumulatedRef);
        if stochastic == false
            validationResult(:,2) = nseValue;
        else
            validationResult{:,2} = nseValue;
        end
    end

    if optimisation == false
        validationMetric.(targetVarL(j)) = validationResult;
    else
        validOptim = mean(validationResult(:,2));
    end
end

switch metricV
    case 1, disp(['  Mean MAE: ' num2str(mean(validationResult(:,2)))])
    case 2, disp(['  Mean RMSE: ' num2str(mean(validationResult(:,2)))])
    case 3, disp(['  Mean SPEM: ' num2str(mean(validationResult(:,2)))])
    case 4, disp(['  Mean SPAEF: ' num2str(mean(validationResult(:,2)))])
    case 5, disp(['  Mean KGE: ' num2str(mean(validationResult(:,2)))])
    case 6, disp(['  Mean NSE: ' num2str(mean(validationResult(:,2)))])

        if optimisation == false
            disp('Saving validationMetric.mat table...')
            validationSave = fullfile(outputDir,'validationMetric.mat');
            save(validationSave, 'validationMetric');
        else % for optimisation run
            validationMetric = validOptim;
        end

end
