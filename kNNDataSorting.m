function sortedDates = kNNDataSorting(climateVars,queryDates,learningDates,climateData,normMethods,shortWindow,longWindow,daysRange,Weights,nbImages,metricKNN,optimPrep,saveOptimPrep,parallelComputing,inputDir,useDOY)

%
%
%
% REDO DOCUMENTATION
%
%
%

%tic

% checks that at least one learning and query dates are present
if any(size(learningDates)==0)
    error('At least one dimension of LearningDates is 0! Code exited...')
elseif any(size(queryDates)==0)
    error('At least one dimension of QueryDates is 0! Code exited...')
end

% the learning dates are ranked based on a criterion that quantifies their distance to a given query date

% adaptation to have only one big climateDataAll file, one learningDates
% file containing date and target and same for queryDates
climateDates = table2array(climateData(:,'date'));
climateCells = table2array(removevars(climateData,'date'));
[~, nVar]    = size(climateCells);
climateMaps  = cell(1, nVar);
% Fill each cell with a [x y time] matrix
for v = 1:nVar
    % Extract the column for variable v and stack along the 3rd dimension
    climateMaps{v} = cat(3, climateCells{:, v});
end

climateVarsNames  = string(removevars(climateData,'date').Properties.VariableNames);

queryDatesDate = queryDates;

learningDatesDate = table2array(learningDates(:,1));

if optimPrep == false
    % Define learningDates as itself minus the query dates
    ismem = ismember(learningDatesDate, queryDatesDate);
    learningDatesDate = learningDatesDate(~ismem);
end

totQDates = size(queryDatesDate,1);
totLDates = size(learningDatesDate,1);

sortedDates   = cell(totQDates, 1);
sortedData    = cell(totQDates, 1);
sortedDist    = cell(totQDates, 1);

% Assign different weights
idxShort      = contains(Weights.Properties.VariableNames,'Short');
weightsShort  = table2array(Weights(:,idxShort));
idxLong       = contains(Weights.Properties.VariableNames,'Long');
weightsLong   = table2array(Weights(:,idxLong));
if metricKNN == 5
    idxSpem    = contains(Weights.Properties.VariableNames,'SPEM');
    spemW      = table2array(Weights(:,idxSpem));
    idxHel     = contains(Weights.Properties.VariableNames,'Hellinger');
    helW       = table2array(Weights(:,idxHel));
    spemHelW   = cat(2, spemW, helW);
else
    spemHelW   = [];
end
if useDOY
    idxDOY    = contains(Weights.Properties.VariableNames,'DOY');
    weightDOY = table2array(Weights(:,idxDOY));
end

disp('Starting loop to sort learning dates for each query date...')

% Display progression - for parallel computing
%progress = 0;
%fprintf(1,'Progress: %3.0f%%\n',progress);
%fprintf(['\n' repmat('.',1,totQDates) '\n\n']);

if parallelComputing == true
    parfor qd = 1:totQDates % parallel computing
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

        disp(['  Processing day ' num2str(qd) '/' num2str(totQDates) ' (' num2str(currentQDate) ')'])

        % Extract the longWindow climate for the current query date
        idx = find(climateDates == currentQDate);
        if idx > longWindow
            queryClimate = cellfun(@(M) M(:,:,idx-(longWindow-1):idx), climateMaps, 'UniformOutput', false);

           % Compute the distances between the query climate and the climate for each learning date
            %targetDistance  = cell(totLDates,1);
            climateDistance = cell(totLDates,2);
            % Display progress - only for serial computing
            %fprintf(1,'    Progress for current query date: %3.0f%%\n',progress);
            for ld = 1:totLDates
                currentLDate    = learningDatesDate(ld);
                dayOfYearL      = day(datetime(currentLDate,'ConvertFrom','yyyyMMdd'),'dayofyear');
                if dayOfYearL == 366
                    dayOfYearL = 1;
                end
                idx = find(climateDates == currentLDate);
                %disp(['    Computing distance to day ' num2str(l) '/' num2str(totLDates) ' (' num2str(currentLDate) ')'])
                if ismember(dayOfYearL,rangeQ) % if learning date is not within 3 months of the query date, it is skipped
                    if idx >= longWindow % skips learning dates that are in the longWindow
                        % Learning dates climate
                        learningClimate = cellfun(@(M) M(:,:,idx-(longWindow-1):idx), climateMaps, 'UniformOutput', false);

                        % Climate distance
                        % 1 date, 2 distance
                        climateDistAll = cell(ld,1);
                        climVarIdx     = zeros(1,numel(climateVarsNames));
                        otherIdx       = ~climVarIdx;

                        % DOY distance
                        if useDOY
                            absDiff = abs(dayOfYearQ - dayOfYearL);
                            doyDist = min(absDiff, 365-absDiff)/daysRange;
                        end

                        if sum(ismember(normMethods,4))>=1
                            normMeIdx  = find(normMethods==4);
                            climVarIdx = strcmpi(climateVars(normMeIdx),climateVarsNames);
                            otherIdx   = ~climVarIdx;
                            binaryLClim = cellfun(@(x) (x>0)+isnan(x).*x, learningClimate(:,climVarIdx), 'UniformOutput', false);
                            binaryQClim = cellfun(@(x) (x>0)+isnan(x).*x, queryClimate(:,climVarIdx), 'UniformOutput', false);
                            hammingDist = cellfun(@(x,y) single(mean(x(:) ~= y(:))),binaryLClim,binaryQClim,'UniformOutput',false); % Hamming distance
                            hellingDist = cellfun(@(x,y) single(hellingerDist(x(x>0), y(y>0))), ...
                                learningClimate(:,climVarIdx),queryClimate(:,climVarIdx),'UniformOutput',false); % Hellinger distance
                            climateDistHH = cellfun(@(x,y) (x/2)+(y/2),hammingDist,hellingDist,'UniformOutput',false); % Total Hamming + Hellinger distance
                            climateDistAll{ld,1}(:,climVarIdx) = climateDistHH(:);
                        end
%                         learningSubset = learningClimate(:, otherIdx);
%                         querySubset    = queryClimate(:, otherIdx);
                        if metricKNN == 1 % RMSE
                            climateDistAll{ld,1}(:, otherIdx) = cell2mat(cellfun(@(x, y) squeeze(sqrt(mean(mean((x - y).^2, 1, 'omitnan'), 2, 'omitnan'))), ...
                                learningClimate(:, otherIdx), queryClimate(:, otherIdx), 'UniformOutput', false));
                        elseif metricKNN == 2 % MAE
                            climateDistAll{ld,1}(:, otherIdx) = cell2mat(cellfun(@(x, y) squeeze(mean(mean(abs(x - y), 1, 'omitnan'), 2, 'omitnan')), ...
                                learningClimate(:, otherIdx), queryClimate(:, otherIdx), 'UniformOutput', false));
                        elseif metricKNN == 3 % 1-bSPEM
                            climateDistAll{ld,1}(:, otherIdx) = cell2mat(cellfun(@(x,y) squeeze((1 - spem(x, y))), ...
                                learningClimate(:, otherIdx), queryClimate(:, otherIdx), 'UniformOutput', false));
                        elseif metricKNN == 4 % Hellinger
                            climateDistAll{ld,1}(:, otherIdx) = cell2mat(cellfun(@(x,y) hellingerDist(x, y), ...
                                learningClimate(:, otherIdx), queryClimate(:, otherIdx), 'UniformOutput', false));
                        elseif metricKNN == 5 % 0.5*(1-bSPEM) + 0.5*Hellinger
                            if optimPrep == false
                                climateDistAll{ld,1}(:, otherIdx) = cell2mat(cellfun(@(x, y) ...
                                    computeHellingerSPEM(x, y, @spem, @hellingerDist, spemHelW), ...
                                    learningClimate(:, otherIdx), queryClimate(:, otherIdx), 'UniformOutput', false));
                            else
                                climateDistAll{ld,1}(:, otherIdx) = cell2mat(cellfun(@(x,y) (1 - spem(x, y)), ...
                                    learningClimate(:, otherIdx), queryClimate(:, otherIdx), 'UniformOutput', false));
                                climateDistAll{ld,2}(:, otherIdx) = cell2mat(cellfun(@(x,y) hellingerDist(x, y), ...
                                    learningClimate(:, otherIdx), queryClimate(:, otherIdx), 'UniformOutput', false));
                            end
                        elseif metricKNN == 6 % SPAEF
                            climateDistAll{ld,1}(:, otherIdx) = cell2mat(cellfun(@(x,y) (1 - spaef(x, y)), ...
                                learningClimate(:, otherIdx), queryClimate(:, otherIdx), 'UniformOutput', false));
                        else
                            error('Bad metricKNN parameter')
                        end
                        climateDistance{ld,1} = currentLDate;
                        if shortWindow > 0
                            climateDistance{ld,2}(1,:) = sum(climateDistAll{ld,1}(1:shortWindow,:),1);
                            if optimPrep == true && metricKNN == 5
                                climateDistance{ld,3}(1,:) = sum(climateDistAll{ld,2}(1:shortWindow,:),1);
                            end
                        else
                            climateDistance{ld,2}(1,:) = zeros(1,size(climateDistAll{ld,1},2));
                            if optimPrep == true && metricKNN == 5
                                climateDistance{ld,3}(1,:) = zeros(1,size(climateDistAll{ld,2},2));
                            end
                        end
                        if optimPrep == false 
                            climateDistance{ld,2}(2,:) = sum(climateDistAll{ld,1}(shortWindow+1:end,:),1);
                        elseif optimPrep == true && metricKNN == 5
                            climateDistance{ld,2}(2,:) = sum(climateDistAll{ld,1}(shortWindow+1:end,:),1);
                            climateDistance{ld,3}(2,:) = sum(climateDistAll{ld,2}(shortWindow+1:end,:),1);
                            if useDOY
                                climateDistance{ld,4} = weightDOY * doyDist;
                            end
                        end
                        % Assign weights to corresponding index
                        if optimPrep == false
                            climateDistance{ld,2}(1,:) = climateDistance{ld,2}(1,:) .* weightsShort;
                            climateDistance{ld,2}(2,:) = climateDistance{ld,2}(2,:) .* weightsLong;
                            climateDistance{ld,2} = sum(climateDistance{ld,2},'all');
                            if useDOY
                                climateDistance{ld,2} = climateDistance{ld,2} + (weightDOY * doyDist);
                            end
                        end
                    else
                        % If not enough climate days available, skip until loop reaches longWindow
                        %warning(['Climate data available is shorter than longWindow, ' num2str(currentLDate) ' skipped.'])
                        continue
                    end
                else
                    % If learning date not in query date range, skip it
                    %disp(['    Learning day ',num2str(currentLDate),' not in query date range, skipped'])
                    continue
                end
            end
        else
            fprintf('\n')
            disp('    Not enough learning dates. Query date skipped...')
            fprintf('\b')
            continue
        end

        % Learning dates distance: 1 date, 2 distance
        distance = climateDistance(~cellfun('isempty',climateDistance(:,1)),:);
        %targetDistance  = targetDistance(~cellfun('isempty',targetDistance),:);
        if optimPrep == false
            distancesSort   = sortrows(distance,2); % Sort rows in ascending order according to column 2
            distancesDates  = distancesSort(1:nbImages,1);
            distSorted      = distancesSort(1:nbImages,2);
            sortedDates{qd} = currentQDate;
            sortedData{qd}  = cell2mat(distancesDates);
            sortedDist{qd}  = cell2mat(distSorted);
        else
            distancesDates  = distance(:,1);
            distSorted      = distance(:,2:end);
            sortedDates{qd} = currentQDate;
            sortedData{qd}  = distancesDates;
            sortedDist{qd}  = distSorted;
        end
    end

    sortedDatesAll = [sortedDates sortedData sortedDist];
    sortedDates    = sortedDatesAll;

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

        % Extract the longWindow climate for the current query date
        queryClimate = cell(longWindow, numel(climateVars));
        idx = find(climateDates == currentQDate);
        if idx > longWindow
            queryClimate = cellfun(@(M) M(:,:,idx-(longWindow-1):idx), climateMaps, 'UniformOutput', false);
            % Compute the distances between the query climate and the climate for each learning date
            %targetDistance  = cell(totLDates,1);
            climateDistance = cell(totLDates,2);
            if optimPrep == false
                climateDistAll  = inf(longWindow,numel(climateVars));
            elseif optimPrep == true && metricKNN == 5
                climateDistAll  = cell(1,2);
            end
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
                %disp(['    Computing distance to day ' num2str(l) '/' num2str(totLDates) ' (' num2str(currentLDate) ')'])
                if ismember(dayOfYearL,rangeQ) % if learning date is not within 3 months of the query date, it is skipped
                    %disp(['    Processing learning day ', num2str(currentLDate)])
                    if idx >= longWindow % skips learning dates that are in the longWindow
                        % Learning dates climate
                        learningClimate = cellfun(@(M) M(:,:,idx-(longWindow-1):idx), climateMaps, 'UniformOutput', false);

                        % DOY distance
                        if useDOY
                            absDiff = abs(dayOfYearQ - dayOfYearL);
                            doyDist = min(absDiff, 365-absDiff)/daysRange;
                        end

                        % Climate distance
                        % 1 date, 2 distance
                        climVarIdx = zeros(1,numel(climateVarsNames));
                        otherIdx   = ~climVarIdx;
                        if sum(ismember(normMethods,4))>=1
                            normMeIdx  = find(normMethods==4);
                            climVarIdx = strcmpi(climateVars(normMeIdx),climateVarsNames);
                            otherIdx   = ~climVarIdx;
                            binaryLClim = cellfun(@(x) (x>0)+isnan(x).*x, learningClimate(:,climVarIdx), 'UniformOutput', false);
                            binaryQClim = cellfun(@(x) (x>0)+isnan(x).*x, queryClimate(:,climVarIdx), 'UniformOutput', false);
                            hammingDist = cellfun(@(x,y) single(mean(x(:) ~= y(:))),binaryLClim,binaryQClim,'UniformOutput',false); % Hamming distance
                            %hellingDist = cellfun(@(x, y) 1/sqrt(2)*sqrt(sum((sqrt(x(x > 0)) - sqrt(y(y > 0))).^2)), ...
                            %    learningClimate(:,climVarIdx),queryClimate(:,climVarIdx),'UniformOutput',false); % Hellinger distance
                            hellingDist = cellfun(@(x,y) single(hellingerDist(x(x>0), y(y>0))), ...
                                learningClimate(:,climVarIdx),queryClimate(:,climVarIdx),'UniformOutput',false); % Hellinger distance
                            climateDistHH = cellfun(@(x,y) (x/2)+(y/2),hammingDist,hellingDist,'UniformOutput',false); % Total Hamming + Hellinger distance
                            climateDistAll(:,climVarIdx) = climateDistHH(:);
                        end
%                         learningSubset = learningClimate(:, otherIdx);
%                         querySubset    = queryClimate(:, otherIdx);
                        if metricKNN == 1 % RMSE
                            climateDistAll(:, otherIdx) = cell2mat(cellfun(@(x, y) squeeze(sqrt(mean(mean((x - y).^2, 1, 'omitnan'), 2, 'omitnan'))), ...
                                learningClimate(:, otherIdx), queryClimate(:, otherIdx), 'UniformOutput', false));
                        elseif metricKNN == 2 % MAE
                            climateDistAll(:, otherIdx) = cell2mat(cellfun(@(x, y) squeeze(mean(mean(abs(x - y), 1, 'omitnan'), 2, 'omitnan')), ...
                                learningClimate(:, otherIdx), queryClimate(:, otherIdx), 'UniformOutput', false));
                        elseif metricKNN == 3 % 1-bSPEM
                            climateDistAll(:, otherIdx) = cell2mat(cellfun(@(x, y) (1 - spem(x, y)), ...
                                    learningClimate(:, otherIdx), queryClimate(:, otherIdx), 'UniformOutput', false));
                        elseif metricKNN == 4 % Hellinger
                            climateDistAll(:, otherIdx) = cell2mat(cellfun(@(x, y) hellingerDist(x, y), ...
                                    learningClimate(:, otherIdx), queryClimate(:, otherIdx), 'UniformOutput', false));
                        elseif metricKNN == 5 % 0.5*(1-bSPEM) + 0.5*Hellinger
                            if optimPrep == false
                                climateDistAll(:, otherIdx) = cell2mat(cellfun(@(x, y) ...
                                    computeHellingerSPEM(x, y, @spem, @hellingerDist, spemHelW), ...
                                    learningClimate(:, otherIdx), queryClimate(:, otherIdx), 'UniformOutput', false));
                            else
                                climateDistAll{:,1}(:, otherIdx) = cell2mat(cellfun(@(x,y) (1 - spem(x, y)), ...
                                    learningClimate(:, otherIdx), queryClimate(:, otherIdx), 'UniformOutput', false));
                                climateDistAll{:,2}(:, otherIdx) = cell2mat(cellfun(@(x,y) hellingerDist(x, y), ...
                                    learningClimate(:, otherIdx), queryClimate(:, otherIdx), 'UniformOutput', false));
                            end
                        elseif metricKNN == 6 % SPAEF
                            climateDistAll(:, otherIdx) = cell2mat(cellfun(@(x, y) (1 - spaef(x, y)), ...
                                    learningClimate(:, otherIdx), queryClimate(:, otherIdx), 'UniformOutput', false));
                        else
                            error('Bad metricKNN parameter')
                        end
                        climateDistance{ld,1} = currentLDate;
           
                        if shortWindow > 0
                            if optimPrep == false
                                climateDistance{ld,2}(1,:) = sum(climateDistAll(1:shortWindow,:),1);
                            elseif optimPrep == true && metricKNN == 5
                                climateDistance{ld,2}(1,:) = sum(climateDistAll{:,1}(1:shortWindow,:),1);
                                climateDistance{ld,3}(1,:) = sum(climateDistAll{:,2}(1:shortWindow,:),1);
                            end
                        else
                            if optimPrep == false
                                climateDistance{ld,2}(1,:) = zeros(1,size(climateDistAll,2));
                            elseif optimPrep == true && metricKNN == 5
                                climateDistance{ld,2}(1,:) = zeros(1,size(climateDistAll{:,1},2));
                                climateDistance{ld,3}(1,:) = zeros(1,size(climateDistAll{:,2},2));
                            end
                        end
                        if optimPrep == false 
                            climateDistance{ld,2}(2,:) = sum(climateDistAll(shortWindow+1:end,:),1);
                        elseif optimPrep == true && metricKNN == 5
                            climateDistance{ld,2}(2,:) = sum(climateDistAll{:,1}(shortWindow+1:end,:),1);
                            climateDistance{ld,3}(2,:) = sum(climateDistAll{:,2}(shortWindow+1:end,:),1);
                            if useDOY
                                climateDistance{ld,4} = weightDOY * doyDist;
                            end
                        end

                        % Assign weights to corresponding index
                        if optimPrep == false
                            climateDistance{ld,2}(1,:) = climateDistance{ld,2}(1,:) .* weightsShort;
                            climateDistance{ld,2}(2,:) = climateDistance{ld,2}(2,:) .* weightsLong;
                            climateDistance{ld,2} = sum(climateDistance{ld,2},'all');
                            if useDOY
                                climateDistance{ld,2} = climateDistance{ld,2} + (weightDOY * doyDist);
                            end
                        end
                    else
                        % If not enough climate days available, skip until loop reaches longWindow
                        %warning(['Climate data available is shorter than longWindow, ' num2str(currentLDate) ' skipped.'])
                        continue
                    end
                else
                    % If learning date not in query date range, skip it
                    % Display computation progress - only for serial computing
                    progress = (100*(ld/totLDates));
                    fprintf(1,'\b\b\b\b%3.0f%%',progress);
                    %disp(['    Learning day ',num2str(currentLDate),' not in query date range, skipped'])
                    continue
                end
                % Display computation progress - only for serial computing
                progress = (100*(ld/totLDates));
                fprintf(1,'\b\b\b\b%3.0f%%',progress);
            end
        else
            fprintf('\n')
            disp('    Not enough learning dates. Query date skipped...')
            fprintf('\b')
            continue
        end

        % Learning dates distance: 1 date, 2 distance
        distance = climateDistance(~cellfun('isempty',climateDistance(:,1)),:);
        %targetDistance  = targetDistance(~cellfun('isempty',targetDistance),:);
        if optimPrep == false
            distancesSort   = sortrows(distance,2); % Sort rows in ascending order according to column 2
            distancesDates  = distancesSort(1:nbImages,1);
            distSorted      = distancesSort(1:nbImages,2);
            sortedDates{qd} = currentQDate;
            sortedData{qd}  = cell2mat(distancesDates);
            sortedDist{qd}  = cell2mat(distSorted);
        else
            distancesDates  = distance(:,1);
            distSorted      = distance(:,2:end);
            sortedDates{qd} = currentQDate;
            sortedData{qd}  = distancesDates;
            sortedDist{qd}  = distSorted;
        end

        % Display progression - for parallel computing
        %progress = (100*(l/totLDates));
        %fprintf(1,'\b\b\b\b%3.0f%%',progress);

        %toc
    end

    sortedDatesAll = [sortedDates sortedData sortedDist];
    sortedDates    = sortedDatesAll;
    
    fprintf('\n')
end

sortedDates = sortedDates(~cellfun('isempty',sortedDates(:,1)),:);

if optimPrep == false %&& saveMats == true
    disp('Saving KNNSorting.mat file...')
    save(fullfile(inputDir,'KNNSorting.mat'),'sortedDates', '-v7.3','-nocompression'); % Save Ranked Learning Dates per Query Date
else
    if saveOptimPrep == true
        disp('Saving KNNDistances.mat file for optimisation. May take a while...')
        save(fullfile(inputDir,'KNNDistances.mat'),'sortedDates', '-v7.3','-nocompression');
    end
end

%toc

end
