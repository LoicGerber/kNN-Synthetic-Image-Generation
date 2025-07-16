function sortedDates = pixelWise_kNNDataSorting(maskDir,climateVars,queryDates,learningDates,climateData,longWindow,daysRange,nbImages,metricKNN,optimPrep,saveOptimPrep,parallelComputing,inputDir)

%
%
%
% REDO DOCUMENTATION
%
%
%

maskData = readgeoraster(maskDir);

% checks that at least one learning and query dates are present
if any(size(learningDates)==0)
    error('At least one dimension of LearningDates is 0! Code exited...')
elseif any(size(queryDates)==0)
    error('At least one dimension of QueryDates is 0! Code exited...')
end

climateDates = table2array(climateData(:,'date'));
climateCells = table2array(removevars(climateData,'date'));
% Get sizes
[nTime, nVar]  = size(climateCells);
[xSize, ySize] = size(climateCells{1,1});  % Assuming all cells contain same-size [x y] matrices
% Flatten each [x y] matrix into a column vector
flatData       = cellfun(@(m) reshape(m, [], 1), climateCells, 'UniformOutput', false);
% Stack all vectors horizontally: size will be [x*y x (nTime*nVar)]
stacked        = cat(2, flatData{:});
% Reshape into 4D array: [x y time variable]
climateMaps    = reshape(stacked, xSize, ySize, nTime, nVar);

queryDatesDate    = queryDates;
learningDatesDate = table2array(learningDates(:,1));

if optimPrep == false
    ismem = ismember(learningDatesDate, queryDatesDate);
    learningDatesDate = learningDatesDate(~ismem);
end

totQDates = size(queryDatesDate,1);
totLDates = size(learningDatesDate,1);
distanceAll = cell(size(maskData,1),size(maskData,2),totQDates);

disp('Starting loop to sort learning dates for each query date...')

if parallelComputing == true
    parfor qd = 1:totQDates
        currentQDate = queryDatesDate(qd);
        dayOfYearQ = day(datetime(currentQDate,'ConvertFrom','yyyyMMdd'),'dayofyear');
        minRangeQ = dayOfYearQ - daysRange;
        if minRangeQ <= 0, minRangeQ = 365 + minRangeQ; end
        maxRangeQ = dayOfYearQ + daysRange;
        if maxRangeQ > 365, maxRangeQ = maxRangeQ - 365; end
        if minRangeQ < maxRangeQ
            rangeQ = minRangeQ:1:maxRangeQ;
        else
            rangeQmin = minRangeQ:1:365;
            rangeQmax = 1:1:maxRangeQ;
            rangeQ    = [rangeQmin rangeQmax];
        end

        fprintf(['\n  Processing day ' num2str(qd) '/' num2str(totQDates) ' (' num2str(currentQDate) ')'])
        
        distMapQd = cell(size(maskData));

        % Extract the longWindow climate for the current query date
        idx = find(climateDates == currentQDate);
        if idx > longWindow
            queryClimate = climateMaps(:, :, idx-(longWindow-1):idx, :);

            % Compute the distances between the query climate and the climate for each learning date
            climateDistance = cell(totLDates,2);
            % Display progress - only for serial computing
            for ld = 1:totLDates
                currentLDate    = learningDatesDate(ld);
                dayOfYearL      = day(datetime(currentLDate,'ConvertFrom','yyyyMMdd'),'dayofyear');
                if dayOfYearL == 366
                    dayOfYearL = 1;
                end
                idx = find(climateDates == currentLDate);
                if ismember(dayOfYearL,rangeQ) % if learning date is not within 3 months of the query date, it is skipped
                    if idx >= longWindow % skips learning dates that are in the longWindow
                        % Learning dates climate
                        learningClimate = climateMaps(:, :, idx-(longWindow-1):idx, :);
                        
                        % Compute absolute difference and mean across time (3rd dimension)
                        if metricKNN == 1 % RMSE
                            climateDistAll = sqrt(mean((learningClimate - queryClimate).^2, 3, 'omitnan'));
                        elseif metricKNN == 2 % MAE
                            climateDistAll = mean(abs(learningClimate - queryClimate), 3, 'omitnan');
                        else
                            error('Bad metricKNN parameter')
                        end
                        climateDistance{ld,1} = currentLDate;
                        climateDistance{ld,2} = sum(sum(climateDistAll, 3, 'omitnan'), 4, 'omitnan');
                    else
                        continue
                    end
                else
                    continue
                end
            end

            distance = climateDistance(~cellfun('isempty',climateDistance(:,1)),:);
            
            % Get matrix size from the first entry
            [X, Y] = size(distance{1, 2});
            
            % Pre-extract all dates and matrices
            allDates    = distance(:, 1);
            allMatrices = cat(3, distance{:, 2});
            
            % Loop through each pixel position
            for pixY = 1:X
                for pixX = 1:Y
                    if maskData(pixY, pixX) == 0
                        % Skip masked pixels
                        %distMapQd{pixY, pixX} = NaN;
                        continue;
                    end
                    % Extract the time series for pixel (r, c)
                    pixelValues = squeeze(allMatrices(pixY, pixX, :));
                    [sortedVals, idx] = sort(pixelValues, 'ascend');

                    % Sort dates accordingly
                    sortDates = cell2mat(allDates(idx));
            
                    % Keep only the nbImages first
                    sortedVals = sortedVals(1:min(nbImages, end));
                    sortDates  = sortDates(1:min(nbImages, end));
            
                    % Store as cell array of [date, value] pairs
                    distMapQd{pixY, pixX} = [sortDates, double(sortedVals)];
                end
            end
        else
            fprintf('\n')
            disp('    Not enough learning dates. Query date skipped...')
            fprintf('\b')
            continue
        end
        distanceAll(:,:,qd) = distMapQd;
    end
    fprintf('\n')
    % Shut down parallel pool
    poolobj = gcp('nocreate');
    delete(poolobj);
else % serial computing
    for qd = 1:totQDates
        currentQDate = queryDatesDate(qd);
        dayOfYearQ = day(datetime(currentQDate,'ConvertFrom','yyyyMMdd'),'dayofyear');
        minRangeQ = dayOfYearQ - daysRange;
        if minRangeQ <= 0, minRangeQ = 365 + minRangeQ; end
        maxRangeQ = dayOfYearQ + daysRange;
        if maxRangeQ > 365, maxRangeQ = maxRangeQ - 365; end
        if minRangeQ < maxRangeQ
            rangeQ = minRangeQ:1:maxRangeQ;
        else
            rangeQmin = minRangeQ:1:365;
            rangeQmax = 1:1:maxRangeQ;
            rangeQ    = [rangeQmin rangeQmax];
        end

        fprintf(['\n  Processing day ' num2str(qd) '/' num2str(totQDates) ' (' num2str(currentQDate) ')'])
        
        distMapQd = cell(size(maskData));

        % Extract the longWindow climate for the current query date
        idx = find(climateDates == currentQDate);
        if idx > longWindow
            queryClimate = climateMaps(:, :, idx-(longWindow-1):idx, :);

            % Compute the distances between the query climate and the climate for each learning date
            climateDistance = cell(totLDates,2);
            % Display progress - only for serial computing
            progress = 0;
            fprintf(1,'\n    Progress for current query date: %3.0f%%\n',progress);
            for ld = 1:totLDates
                currentLDate    = learningDatesDate(ld);
                dayOfYearL      = day(datetime(currentLDate,'ConvertFrom','yyyyMMdd'),'dayofyear');
                if dayOfYearL == 366
                    dayOfYearL = 1;
                end
                idx = find(climateDates == currentLDate);
                if ismember(dayOfYearL,rangeQ) % if learning date is not within 3 months of the query date, it is skipped
                    if idx >= longWindow % skips learning dates that are in the longWindow
                        % Learning dates climate
                        learningClimate = climateMaps(:, :, idx-(longWindow-1):idx, :);
                        
                        % Compute absolute difference and mean across time (3rd dimension)
                        if metricKNN == 1 % RMSE
                            climateDistAll = sqrt(mean((learningClimate - queryClimate).^2, 3, 'omitnan'));
                        elseif metricKNN == 2 % MAE
                            climateDistAll = mean(abs(learningClimate - queryClimate), 3, 'omitnan');
                        else
                            error('Bad metricKNN parameter')
                        end
                        climateDistance{ld,1} = currentLDate;
                        climateDistance{ld,2} = sum(sum(climateDistAll, 3, 'omitnan'), 4, 'omitnan');
                    else
                        continue
                    end
                else
                    continue
                end
                % Display computation progress - only for serial computing
                progress = (100*(ld/totLDates));
                fprintf(1,'\b\b\b\b%3.0f%%',progress);
            end

            distance = climateDistance(~cellfun('isempty',climateDistance(:,1)),:);
            
            % Get matrix size from the first entry
            [X, Y] = size(distance{1, 2});
            
            % Pre-extract all dates and matrices
            allDates    = distance(:, 1);
            allMatrices = cat(3, distance{:, 2});
            
            % Loop through each pixel position
            for pixY = 1:X
                for pixX = 1:Y
                    if maskData(pixY, pixX) == 0
                        % Skip masked pixels
                        %distMapQd{pixY, pixX} = NaN;
                        continue;
                    end
                    % Extract the time series for pixel (r, c)
                    pixelValues = squeeze(allMatrices(pixY, pixX, :));
                    [sortedVals, idx] = sort(pixelValues, 'ascend');

                    % Sort dates accordingly
                    sortDates = cell2mat(allDates(idx));
            
                    % Keep only the nbImages first
                    sortedVals = sortedVals(1:min(nbImages, end));
                    sortDates  = sortDates(1:min(nbImages, end));
            
                    % Store as cell array of [date, value] pairs
                    distMapQd{pixY, pixX} = [sortDates, double(sortedVals)];
                end
            end
        else
            fprintf('\n')
            disp('    Not enough learning dates. Query date skipped...')
            fprintf('\b')
            continue
        end
        distanceAll(:,:,qd) = distMapQd;
    end    
    fprintf('\n')
end

sortedDates.data = distanceAll;
sortedDates.date = queryDatesDate;

if optimPrep == false %&& saveMats == true
    disp('Saving KNNSorting.mat file...')
    save(fullfile(inputDir,'KNNSorting.mat'),'sortedDates', '-v7.3','-nocompression'); % Save Ranked Learning Dates per Query Date
else
    if saveOptimPrep == true
        disp('Saving KNNDistances.mat file for optimisation. May take a while...')
        save(fullfile(inputDir,'KNNDistances.mat'),'sortedDates', '-v7.3','-nocompression');
    end
end

end
