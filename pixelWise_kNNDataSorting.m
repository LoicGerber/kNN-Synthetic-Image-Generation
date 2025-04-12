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

climateDates      = table2array(climateData(:,'date'));
climateMaps       = table2array(removevars(climateData,'date'));

queryDatesDate = table2array(queryDates(:,1));
learningDatesDate = table2array(learningDates(:,1));

if optimPrep == false
    % Define learningDates as itself minus the query dates
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
        queryClimate = cell(longWindow, numel(climateVars));
        idx = find(climateDates == currentQDate);
        if idx > longWindow
            kj = 1;
            for k = (longWindow-1):-1:0
                queryClimate(kj,:) = climateMaps(idx-k,:);
                kj = kj+1;
            end

            % Compute the distances between the query climate and the climate for each learning date
            climateDistance = cell(totLDates,2);
            % Display progress - only for serial computing
            for ld = 1:totLDates
                learningClimate = cell(longWindow, numel(climateVars));
                currentLDate    = learningDatesDate(ld);
                dayOfYearL      = day(datetime(currentLDate,'ConvertFrom','yyyyMMdd'),'dayofyear');
                if dayOfYearL == 366
                    dayOfYearL = 1;
                end
                idx = find(climateDates == currentLDate);
                if ismember(dayOfYearL,rangeQ) % if learning date is not within 3 months of the query date, it is skipped
                    if idx >= longWindow % skips learning dates that are in the longWindow
                        % Learning dates climate
                        kj = 1;
                        for k = (longWindow-1):-1:0
                            learningClimate(kj,:) = climateMaps(idx-k,:);
                            kj = kj+1;
                        end

                        % Concatenate time dimension (rows) into 3D arrays per variable
                        learningStack = cellfun(@(col) cat(3, col{:}), num2cell(learningClimate, 1), 'UniformOutput', false);
                        queryStack    = cellfun(@(col) cat(3, col{:}), num2cell(queryClimate, 1), 'UniformOutput', false);
                        
                        % Compute absolute difference and mean across time (3rd dimension)
                        if metricKNN == 1 % RMSE
                            climateDistAll = cellfun(@(a, b) sqrt(mean((a - b).^2, 3, 'omitnan')), ...
                                                     learningStack, queryStack, ...
                                                     'UniformOutput', false);
                        elseif metricKNN == 2 % MAE
                            climateDistAll = cellfun(@(a, b) mean(abs(a - b), 3, 'omitnan'), ...
                                                     learningStack, queryStack, ...
                                                     'UniformOutput', false);
                        else
                            error('Bad metricKNN parameter')
                        end
                        climateDistance{ld,1} = currentLDate;
                        climateDistance{ld,2} = sum(cat(3, climateDistAll{:}), 3, 'omitnan');  
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
        queryClimate = cell(longWindow, numel(climateVars));
        idx = find(climateDates == currentQDate);
        if idx > longWindow
            kj = 1;
            for k = (longWindow-1):-1:0
                queryClimate(kj,:) = climateMaps(idx-k,:);
                kj = kj+1;
            end

            % Compute the distances between the query climate and the climate for each learning date
            climateDistance = cell(totLDates,2);
            % Display progress - only for serial computing
            progress = 0;
            fprintf(1,'\n    Progress for current query date: %3.0f%%\n',progress);
            for ld = 1:totLDates
                learningClimate = cell(longWindow, numel(climateVars));
                currentLDate    = learningDatesDate(ld);
                dayOfYearL      = day(datetime(currentLDate,'ConvertFrom','yyyyMMdd'),'dayofyear');
                if dayOfYearL == 366
                    dayOfYearL = 1;
                end
                idx = find(climateDates == currentLDate);
                if ismember(dayOfYearL,rangeQ) % if learning date is not within 3 months of the query date, it is skipped
                    if idx >= longWindow % skips learning dates that are in the longWindow
                        % Learning dates climate
                        kj = 1;
                        for k = (longWindow-1):-1:0
                            learningClimate(kj,:) = climateMaps(idx-k,:);
                            kj = kj+1;
                        end

                        % Concatenate time dimension (rows) into 3D arrays per variable
                        learningStack = cellfun(@(col) cat(3, col{:}), num2cell(learningClimate, 1), 'UniformOutput', false);
                        queryStack    = cellfun(@(col) cat(3, col{:}), num2cell(queryClimate, 1), 'UniformOutput', false);
                        
                        % Compute absolute difference and mean across time (3rd dimension)
                        if metricKNN == 1 % RMSE
                            climateDistAll = cellfun(@(a, b) sqrt(mean((a - b).^2, 3, 'omitnan')), ...
                                                     learningStack, queryStack, ...
                                                     'UniformOutput', false);
                        elseif metricKNN == 2 % MAE
                            climateDistAll = cellfun(@(a, b) mean(abs(a - b), 3, 'omitnan'), ...
                                                     learningStack, queryStack, ...
                                                     'UniformOutput', false);
                        else
                            error('Bad metricKNN parameter')
                        end
                        climateDistance{ld,1} = currentLDate;
                        climateDistance{ld,2} = sum(cat(3, climateDistAll{:}), 3, 'omitnan');  
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
