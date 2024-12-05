function sortedDates = kNNDataSorting(targetVar,climateVars,addVars,queryDates,learningDates,climateData,additionalVars,normMethods,shortWindow,longWindow,daysRange,Weights,nbImages,metricKNN,optimPrep,saveOptimPrep,parallelComputing,inputDir,saveMats)

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
climateDates      = table2array(climateData(:,'date'));
climateMaps       = table2array(removevars(climateData,'date'));
climateVarsNames  = string(removevars(climateData,'date').Properties.VariableNames);

queryDatesDate = table2array(queryDates(:,1));
queryDatesData = table2array(queryDates(:,2:end));

learningDatesDate = table2array(learningDates(:,1));
learningDatesData = table2array(learningDates(:,2:end));

if optimPrep == false
    % Define learningDates as itself minus the query dates
    ismem = ismember(learningDatesDate, queryDatesDate);
    learningDatesDate = learningDatesDate(~ismem);
    learningDatesData = learningDatesData(~ismem);
end

totQDates = size(queryDatesDate,1);
totLDates = size(learningDatesDate,1);

sortedDates   = cell(totQDates, 1);
sortedData    = cell(totQDates, 1);
sortedTarget  = cell(totQDates, 1);
sortedAddVars = cell(totQDates, 1);
sortedDist    = cell(totQDates, 1);

if ~isempty(additionalVars)
    addVarsDates = table2array(additionalVars(:,'date'));
    addVarsData  = table2array(removevars(additionalVars,'date'));
else
    addVarsDates = [];
    addVarsData  = [];
end

% Assign different weights
idxTarget     = contains(Weights.Properties.VariableNames,targetVar);
weightsTarget = table2array(Weights(:,idxTarget));
idxShort      = contains(Weights.Properties.VariableNames,'Short');
weightsShort  = table2cell(Weights(:,idxShort));
idxLong       = contains(Weights.Properties.VariableNames,'Long');
weightsLong   = table2cell(Weights(:,idxLong));
if ~isempty(additionalVars)
    idxAddVars     = contains(Weights.Properties.VariableNames,addVars);
    weightsAddVars = table2cell(Weights(:,idxAddVars));
else
    weightsAddVars = [];
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
        queryClimate = cell(longWindow, numel(climateVars));
        idx = find(climateDates == currentQDate);
        if idx > longWindow
            kj = 1;
            for k = (longWindow-1):-1:0
                queryClimate(kj,:) = climateMaps(idx-k,:);
                kj = kj+1;
            end

            % Extract the additional data for the current query date
            if ~isempty(additionalVars)
                queryAddVars = cell(1, numel(addVars));
                idx = find(addVarsDates == currentQDate);
                queryAddVars(1,:) = addVarsData(idx,:);
            else
                queryAddVars = [];
            end

            % Compute the distances between the query climate and the climate for each learning date
            %targetDistance  = cell(totLDates,1);
            addVarsDistance = cell(totLDates,1);
            climateDistance = cell(totLDates,2);
            % Display progress - only for serial computing
            %fprintf(1,'    Progress for current query date: %3.0f%%\n',progress);
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
                    if idx >= longWindow % skips learning dates that are in the longWindow
                        % Learning dates climate
                        kj = 1;
                        for k = (longWindow-1):-1:0
                            learningClimate(kj,:) = climateMaps(idx-k,:);
                            kj = kj+1;
                        end

                        % Extract the additional data for the current query date
                        if ~isempty(additionalVars)
                            learningAddVars = cell(1, numel(addVars));
                            idx = find(addVarsDates == currentLDate);
                            learningAddVars(1,:) = addVarsData(idx,:);
                        else
                            learningAddVars = [];
                        end

%                         % Target variable comparison
%                         if ~(isempty(cell2mat(queryDatesData(qd,:))) || unique(isnan(cell2mat(queryDatesData(qd,:))))) %&& sum(sum(cell2mat(queryDatesData(qd,:))))~=0
%                             if metricKNN == 1 % RMSE
%                                 targetDistance{ld} = cellfun(@(x, y) sqrt(mean((x - y).^2, 'all', 'omitnan')), ...
%                                     queryDatesData(qd,:), learningDatesData(ld,:), 'UniformOutput', false); % RMSE
%                             elseif metricKNN == 2 % MAE
%                                 targetDistance{ld} = cellfun(@(x, y) mean(abs(x - y), 'all', 'omitnan'), ...
%                                     queryDatesData(qd,:), learningDatesData(ld,:), 'UniformOutput', false); % MAE
%                             elseif metricKNN == 3 % Manhattan
%                                 targetDistance{ld} = cellfun(@(x, y) sum(abs(x - y), 'all', 'omitnan'), ...
%                                     queryDatesData(qd,:), learningDatesData(ld,:), 'UniformOutput', false); % Manhattan
%                             elseif metricKNN == 4 % Euclidean
%                                 targetDistance{ld} = cellfun(@(x, y) sqrt(sum((x - y).^2, 'all', 'omitnan')), ...
%                                     queryDatesData(qd,:), learningDatesData(ld,:), 'UniformOutput', false); % Euclidean
%                             elseif metricKNN == 5 % SPEM
%                                 targetDistance{ld} = cellfun(@(x, y) spem(x, y), ...
%                                     queryDatesData(qd,:), learningDatesData(ld,:));
%                             else
%                                 error('Bad metricKNN parameter')
%                             end
%                             targetDistance{ld} = sum(cell2mat(targetDistance{ld}),1,'omitnan');
%                             if optimPrep == false
%                                 targetDistance{ld} = targetDistance{ld}.*weightsTarget;
%                                 targetDistance{ld} = sum(cell2mat(targetDistance(ld)),2,'omitnan');
%                             end
%                         else
%                             targetDistance{ld} = 0;
%                         end

                        % Additional variable comparison
                        % 1 distance
                        if ~isempty(addVars) && ~isempty(addVarsData)
                            if ~isempty(addVarsData(qd,:))
                                if metricKNN == 1 % RMSE
                                    addVarsDistance{ld,1} = cellfun(@(x, y) sqrt(mean((x - y).^2, 'all', 'omitnan')), ...
                                        queryAddVars, learningAddVars, 'UniformOutput', false); % RMSE
                                elseif metricKNN == 2 % MAE
                                    addVarsDistance{ld,1} = cellfun(@(x, y) mean(abs(x - y), 'all', 'omitnan'), ...
                                        queryAddVars, learningAddVars, 'UniformOutput', false); % MAE
                                elseif metricKNN == 3 % Manhattan
                                    addVarsDistance{ld,1} = cellfun(@(x, y) sum(abs(x - y), 'all', 'omitnan'), ...
                                        queryAddVars, learningAddVars, 'UniformOutput', false); % Manhattan
                                elseif metricKNN == 4 % Euclidean
                                    addVarsDistance{ld,1} = cellfun(@(x, y) sqrt(sum((x - y).^2, 'all', 'omitnan')), ...
                                        queryAddVars, learningAddVars, 'UniformOutput', false); % Euclidean
                                elseif metricKNN == 5 % SPEM
                                    addVarsDistance{ld,1} = cellfun(@(x, y) spem(x, y), ...
                                        queryAddVars, learningAddVars);
                                else
                                    error('Bad metricKNN parameter')
                                end
                                addVarsDistance{ld,1} = sum(cell2mat(addVarsDistance{ld,1}),1,'omitnan');
                                if optimPrep == false
                                    if numel(addVars) == 1
                                        addVarsDistance{ld,1} = addVarsDistance{ld,1} .* cell2mat(weightsAddVars);
                                    else
                                        addVarsDistance{ld,1} = num2cell(cell2mat(addVarsDistance{ld,1}) .* cell2mat(weightsAddVars));
                                    end
                                end
                            end
                        else
                            addVarsDistance{ld,1} = 0;
                        end

                        % Climate distance
                        % 1 date, 2 distance
                        climateDistAll = cell(ld,1);
                        spemDist       = cell(ld,1);
                        combinedResults = cell(ld,1);
                        hammingDist    = cell(ld,1);
                        hellingDist    = cell(ld,1);
                        climateDistHH  = cell(ld,1);
                        climVarIdx     = zeros(1,numel(climateVarsNames));
                        otherIdx       = ~climVarIdx;
                        if sum(ismember(normMethods,4))>=1
                            normMeIdx  = find(normMethods==4);
                            climVarIdx = strcmpi(climateVars(normMeIdx),climateVarsNames);
                            otherIdx   = ~climVarIdx;
                            binaryLClim = cellfun(@(x) (x>0)+isnan(x).*x, learningClimate(:,climVarIdx), 'UniformOutput', false);
                            binaryQClim = cellfun(@(x) (x>0)+isnan(x).*x, queryClimate(:,climVarIdx), 'UniformOutput', false);
                            hammingDist{ld,1} = cellfun(@(x,y) single(mean(x(:) ~= y(:))),binaryLClim,binaryQClim,'UniformOutput',false); % Hamming distance
                            hellingDist{ld,1} = cellfun(@(x,y) single(hellingerDist(x(x>0), y(y>0))), ...
                                learningClimate(:,climVarIdx),queryClimate(:,climVarIdx),'UniformOutput',false); % Hellinger distance
                            climateDistHH{ld,1} = cellfun(@(x,y) (x/2)+(y/2),hammingDist{ld,1},hellingDist{ld,1},'UniformOutput',false); % Total Hamming + Hellinger distance
                            climateDistAll{ld,1}(:,climVarIdx) = climateDistHH{ld,1}(:);
                        end
                        learningSubset = learningClimate(:, otherIdx);
                        querySubset    = queryClimate(:, otherIdx);
                        if metricKNN == 1 % RMSE
                            climateDistAll{ld,1}(:,otherIdx) = cellfun(@(x, y) sqrt(mean((x - y).^2, 'all', 'omitnan')), ...
                                learningSubset, querySubset, 'UniformOutput', false); % RMSE
                        elseif metricKNN == 2 % MAE
                            climateDistAll{ld,1}(:,otherIdx) = cellfun(@(x, y) mean(abs(x - y), 'all', 'omitnan'), ...
                                learningSubset, querySubset, 'UniformOutput', false); % MAE
                        elseif metricKNN == 3 % Manhattan
                            climateDistAll{ld,1}(:,otherIdx) = cellfun(@(x, y) sum(abs(x - y), 'all', 'omitnan'), ...
                                learningSubset, querySubset, 'UniformOutput', false); % Manhattan
                        elseif metricKNN == 4 % Euclidean
                            climateDistAll{ld,1}(:,otherIdx) = cellfun(@(x, y) sqrt(sum((x - y).^2, 'all', 'omitnan')), ...
                                learningSubset, querySubset, 'UniformOutput', false); % Euclidean
                        elseif metricKNN == 5 % SPEM
                            climateDistAll{ld,1}(:,otherIdx) = cellfun(@(x,y) (1 - spem(x, y)), ...
                                learningSubset, querySubset);
                        elseif metricKNN == 6 % 0.5*SPEM + 0.5*Hellinger
                            %spemDist{ld,1}    = num2cell(cellfun(@(x,y) spem(x, y), learningSubset, querySubset));
                            %hellingDist{ld,1} = cellfun(@(x,y) single(hellingerDist(x,y)), ...
                            %    learningSubset,querySubset,'UniformOutput',false);
                            %climateDistAll{ld,1}(:,otherIdx) = cellfun(@(x,y) ((1-x)/2)+(y/2),spemDist{ld,1},hellingDist{ld,1},'UniformOutput',false);
                            combinedResults{ld,1} = cellfun(@(x, y) ...
                                computeHellingerSPEM(x, y, @spem, @hellingerDist), ...
                                learningSubset, querySubset, 'UniformOutput', false);
                            climateDistAll{ld,1}(:, otherIdx) = combinedResults{ld,1};
                        else
                            error('Bad metricKNN parameter')
                        end
                        climateDistance{ld,1} = currentLDate;
                        if shortWindow > 0
                            climateDistance{ld,2}(1,:) = sum(cell2mat(climateDistAll{ld,1}(1:shortWindow,:)),1,'omitnan');
                        else
                            climateDistance{ld,2}(1,:) = single(zeros(1,size(climateDistAll{ld,1},2)));
                        end
                        climateDistance{ld,2}(2,:) = sum(cell2mat(climateDistAll{ld,1}(shortWindow+1:end,:)),1,'omitnan');
                        % Assign weights to corresponding index
                        if optimPrep == false
                            climateDistance{ld,2}(1,:) = climateDistance{ld,2}(1,:) .* cell2mat(weightsShort);
                            climateDistance{ld,2}(2,:) = climateDistance{ld,2}(2,:) .* cell2mat(weightsLong);
                            climateDistance{ld,2} = sum(climateDistance{ld,2},1,'omitnan');
                            climateDistance{ld,2} = sum(climateDistance{ld,2},2,'omitnan')+addVarsDistance{ld,1}; %targetDistance{ld,1}+
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
        distance        = climateDistance(~cellfun('isempty',climateDistance(:,1)),:);
        %targetDistance  = targetDistance(~cellfun('isempty',targetDistance),:);
        addVarsDistance = addVarsDistance(~cellfun('isempty',addVarsDistance),:);
        if optimPrep == false
            distancesSort   = sortrows(distance,2); % Sort rows in ascending order according to column 2
            distancesBest   = distancesSort(1:nbImages,1);
            distSorted      = distancesSort(1:nbImages,2);
            sortedDates{qd} = currentQDate;
            sortedData{qd}  = cell2mat(distancesBest);
            sortedDist{qd}  = cell2mat(distSorted);
        else
            distancesBest     = distance(:,1);
            distSorted        = distance(:,2);
            sortedDates{qd}   = currentQDate;
            sortedData{qd}    = distancesBest;
            %sortedTarget{qd}  = targetDistance;
            sortedAddVars{qd} = addVarsDistance;
            sortedDist{qd}    = distSorted;
        end
    end

    if optimPrep == false
        sortedDatesAll = [sortedDates sortedData sortedDist];
        sortedDates    = sortedDatesAll;
    else
        sortedDatesAll = [sortedDates sortedData sortedTarget sortedAddVars sortedDist];
        sortedDates    = sortedDatesAll;
    end

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
            kj = 1;
            for k = (longWindow-1):-1:0
                queryClimate(kj,:) = climateMaps(idx-k,:);
                %queryClimate(kj,j) = climateMaps.(varClimate)(:,:,idx-k);
                kj = kj+1;
            end

            % Extract the additional data for the current query date
            if ~isempty(additionalVars)
                queryAddVars = cell(1, numel(addVars));
                idx = find(addVarsDates == currentQDate);
                queryAddVars(1,:) = addVarsData(idx,:);
            else
                queryAddVars = [];
            end

            % Compute the distances between the query climate and the climate for each learning date
            %targetDistance  = cell(totLDates,1);
            addVarsDistance = cell(totLDates,1);
            climateDistance = cell(totLDates,2);
            climateDistAll  = cell(longWindow,numel(climateVars));
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
                        kj = 1;
                        for k = (longWindow-1):-1:0
                            learningClimate(kj,:) = climateMaps(idx-k,:);
                            kj = kj+1;
                        end

                        % Extract the additional data for the current query date
                        if ~isempty(additionalVars)
                            learningAddVars = cell(1, numel(addVars));
                            idx = find(addVarsDates == currentLDate);
                            learningAddVars(1,:) = addVarsData(idx,:);
                        else
                            learningAddVars = [];
                        end

                        % Target variable comparison
%                         if ~(isempty(cell2mat(queryDatesData(qd,:))) || unique(isnan(cell2mat(queryDatesData(qd,:))))) %&& sum(sum(cell2mat(queryDatesData(qd,:))))~=0
%                             if metricKNN == 1 % RMSE
%                                 targetDistance{ld} = cellfun(@(x, y) sqrt(mean((x - y).^2, 'all', 'omitnan')), ...
%                                     queryDatesData(qd,:), learningDatesData(ld,:), 'UniformOutput', false); % RMSE
%                             elseif metricKNN == 2 % MAE
%                                 targetDistance{ld} = cellfun(@(x, y) mean(abs(x - y), 'all', 'omitnan'), ...
%                                     queryDatesData(qd,:), learningDatesData(ld,:), 'UniformOutput', false); % MAE
%                             elseif metricKNN == 3 % Manhattan
%                                 targetDistance{ld} = cellfun(@(x, y) sum(abs(x - y), 'all', 'omitnan'), ...
%                                     queryDatesData(qd,:), learningDatesData(ld,:), 'UniformOutput', false); % Manhattan
%                             elseif metricKNN == 4 % Euclidean
%                                 targetDistance{ld} = cellfun(@(x, y) sqrt(sum((x - y).^2, 'all', 'omitnan')), ...
%                                     queryDatesData(qd,:), learningDatesData(ld,:), 'UniformOutput', false); % Euclidean
%                             elseif metricKNN == 5 % SPEM
%                                 targetDistance{ld} = cellfun(@(x, y) spem(x, y), ...
%                                     queryDatesData(qd,:), learningDatesData(ld,:));
%                             else
%                                 error('Bad metricKNN parameter')
%                             end
%                             targetDistance{ld} = sum(cell2mat(targetDistance{ld}),1,'omitnan');
%                             if optimPrep == false
%                                 targetDistance{ld} = targetDistance{ld}.*weightsTarget;
%                                 targetDistance{ld} = sum(cell2mat(targetDistance(ld)),2,'omitnan');
%                             end
%                         else
%                             targetDistance{ld} = 0;
%                         end

                        % Additional variable comparison
                        % 1 distance
                        if ~isempty(additionalVars) && ~isempty(cell2mat(addVarsData(qd,:)))
                            if metricKNN == 1 % RMSE
                                addVarsDistance{ld} = cellfun(@(x, y) sqrt(mean((x - y).^2, 'all', 'omitnan')), ...
                                    queryAddVars, learningAddVars, 'UniformOutput', false); % RMSE
                            elseif metricKNN == 2 % MAE
                                addVarsDistance{ld} = cellfun(@(x, y) mean(abs(x - y), 'all', 'omitnan'), ...
                                    queryAddVars, learningAddVars, 'UniformOutput', false); % MAE
                            elseif metricKNN == 3 % Manhattan
                                addVarsDistance{ld} = cellfun(@(x, y) sum(abs(x - y), 'all', 'omitnan'), ...
                                    queryAddVars, learningAddVars, 'UniformOutput', false); % Manhattan
                            elseif metricKNN == 4 % Euclidean
                                addVarsDistance{ld} = cellfun(@(x, y) sqrt(sum((x - y).^2, 'all', 'omitnan')), ...
                                    queryAddVars, learningAddVars, 'UniformOutput', false); % Euclidean
                            elseif metricKNN == 5 % SPEM
                                addVarsDistance{ld} = cellfun(@(x, y) spem(x, y), ...
                                    queryAddVars, learningAddVars);
                            else
                                error('Bad metricKNN parameter')
                            end
                            addVarsDistance{ld} = sum(cell2mat(addVarsDistance{ld}),1,'omitnan');
                            if optimPrep == false
                                if numel(addVars) == 1
                                    addVarsDistance{ld} = addVarsDistance{ld} .* cell2mat(weightsAddVars);
                                else
                                    addVarsDistance{ld} = cell2mat(addVarsDistance{ld}) .* cell2mat(weightsAddVars);
                                end
                            end
                        else
                            addVarsDistance{ld} = 0;
                        end

                        % Climate distance
                        % 1 date, 2 distance

                        %                         rmse      = cellfun(@(x, y) sqrt(mean((x - y).^2, 'all', 'omitnan')), ...
                        %                             learningClimate, queryClimate, 'UniformOutput', false); % root mean square error
                        %                         mae       = cellfun(@(x, y) mean(abs(x - y), 'all', 'omitnan'), ...
                        %                             learningClimate, queryClimate, 'UniformOutput', false); % mean absolute error
                        %                         manhattan = cellfun(@(x, y) sum(abs(x - y), 'all', 'omitnan'), ...
                        %                             learningClimate, queryClimate, 'UniformOutput', false); % Manhattan
                        %                         euclidean = cellfun(@(x, y) sqrt(sum((x - y).^2, 'all', 'omitnan')), ...
                        %                             learningClimate, queryClimate, 'UniformOutput', false); % Euclidean
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
                        learningSubset = learningClimate(:, otherIdx);
                        querySubset    = queryClimate(:, otherIdx);
                        if metricKNN == 1 % RMSE
                            climateDistAll(:,otherIdx) = cellfun(@(x, y) sqrt(mean((x - y).^2, 'all', 'omitnan')), ...
                                learningSubset, querySubset, 'UniformOutput', false); % RMSE
                        elseif metricKNN == 2 % MAE
                            climateDistAll(:,otherIdx) = cellfun(@(x, y) mean(abs(x - y), 'all', 'omitnan'), ...
                                learningSubset, querySubset, 'UniformOutput', false); % MAE
                        elseif metricKNN == 3 % Manhattan
                            climateDistAll(:,otherIdx) = cellfun(@(x, y) sum(abs(x - y), 'all', 'omitnan'), ...
                                learningSubset, querySubset, 'UniformOutput', false); % Manhattan
                        elseif metricKNN == 4 % Euclidean
                            climateDistAll(:,otherIdx) = cellfun(@(x, y) sqrt(sum((x - y).^2, 'all', 'omitnan')), ...
                                learningSubset, querySubset, 'UniformOutput', false); % Euclidean
                        elseif metricKNN == 5 % SPEM
                            climateDistAll(:,otherIdx) = num2cell(cellfun(@(x, y) (1 - spem(x, y)), ...
                                    learningSubset, querySubset));
                        elseif metricKNN == 6 % 0.5*SPEM + 0.5*Hellinger
                            %spemDist    = num2cell(cellfun(@(x,y) spem(x, y), learningSubset, querySubset));
                            %hellingDist = cellfun(@(x,y) single(hellingerDist(x, y)), ...
                            %    learningSubset,querySubset,'UniformOutput',false);
                            %climateDistAll(:,otherIdx) = cellfun(@(x,y) ((1-x)/2)+(y/2),spemDist,hellingDist,'UniformOutput',false);
                            combinedResults = cellfun(@(x, y) ...
                                computeHellingerSPEM(x, y, @spem, @hellingerDist), ...
                                learningSubset, querySubset, 'UniformOutput', false);
                            climateDistAll(:, otherIdx) = combinedResults;
                        else
                            error('Bad metricKNN parameter')
                        end
                        climateDistance{ld,1} = currentLDate;
                        if shortWindow > 0
                            climateDistance{ld,2}(1,:) = sum(cell2mat(climateDistAll(1:shortWindow,:)),1,'omitnan');
                        else
                            climateDistance{ld,2}(1,:) = single(zeros(1,size(climateDistAll,2)));
                        end
                        climateDistance{ld,2}(2,:) = sum(cell2mat(climateDistAll(shortWindow+1:end,:)),1,'omitnan');
                        % Assign weights to corresponding index
                        if optimPrep == false
                            climateDistance{ld,2}(1,:) = climateDistance{ld,2}(1,:) .* cell2mat(weightsShort);
                            climateDistance{ld,2}(2,:) = climateDistance{ld,2}(2,:) .* cell2mat(weightsLong);
                            climateDistance{ld,2} = sum(climateDistance{ld,2},1,'omitnan');
                            climateDistance{ld,2} = sum(climateDistance{ld,2},2,'omitnan')+addVarsDistance{ld}; %+targetDistance{ld}
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
        distance        = climateDistance(~cellfun('isempty',climateDistance(:,1)),:);
        %targetDistance  = targetDistance(~cellfun('isempty',targetDistance),:);
        addVarsDistance = addVarsDistance(~cellfun('isempty',addVarsDistance),:);
        if optimPrep == false
            distancesSort   = sortrows(distance,2); % Sort rows in ascending order according to column 2
            distancesBest   = distancesSort(1:nbImages,1);
            distSorted      = distancesSort(1:nbImages,2);
            sortedDates{qd} = currentQDate;
            sortedData{qd}  = cell2mat(distancesBest);
            sortedDist{qd}  = cell2mat(distSorted);
        else
            distancesBest     = distance(:,1);
            distSorted        = distance(:,2);
            sortedDates{qd}   = currentQDate;
            sortedData{qd}    = distancesBest;
            %sortedTarget{qd}  = targetDistance;
            sortedAddVars{qd} = addVarsDistance;
            sortedDist{qd}    = distSorted;
        end

        % Display progression - for parallel computing
        %progress = (100*(l/totLDates));
        %fprintf(1,'\b\b\b\b%3.0f%%',progress);

        %toc
    end

    if optimPrep == false
        sortedDatesAll = [sortedDates sortedData sortedDist];
        sortedDates    = sortedDatesAll;
    else
        sortedDatesAll = [sortedDates sortedData sortedTarget sortedAddVars sortedDist];
        sortedDates    = sortedDatesAll;
    end
    fprintf('\n')
end

sortedDates = sortedDates(~cellfun('isempty',sortedDates(:,1)),:);

if optimPrep == false && saveMats == true
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
