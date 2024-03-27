function sortedDates = pixelWise_kNNDataSorting(maskDir,targetVar,climateVars,addVars,queryDates,learningDates,climateData,additionalVars,normMethods,shortWindow,longWindow,daysRange,Weights,nbImages,metricKNN,optimPrep,saveOptimPrep,parallelComputing,inputDir)

%
%
%
% REDO DOCUMENTATION
%
%
%

%tic

maskData = readgeoraster(maskDir);

% checks that at least one learning and query dates are present
if any(size(learningDates)==0)
    error('At least one dimension of LearningDates is 0! Code exited...')
elseif any(size(queryDates)==0)
    error('At least one dimension of QueryDates is 0! Code exited...')
end

climateDates      = table2array(climateData(:,'date'));
climateMaps       = table2array(removevars(climateData,'date'));
%climateVarsNames  = string(removevars(climateData,'date').Properties.VariableNames);

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
%totPixels = size(maskData,1) * size(maskData,2);
totGPixels = sum(sum(maskData));
distanceAll = cell(size(maskData,1),size(maskData,2),totQDates);

%sortedDates   = cell(totQDates, 1);
%sortedData    = cell(totQDates, 1);
%sortedTarget  = cell(totQDates, 1);
%sortedAddVars = cell(totQDates, 1);
%sortedDist    = cell(totQDates, 1);

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
    %% parallel computing
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

        disp(['  Processing day ' num2str(qd) '/' num2str(totQDates) ' (' num2str(currentQDate) ')'])

        % Extract the longWindow climate for the current query date
        distMapQd = cell(size(maskData));

        queryClimate = nan(longWindow, numel(climateVars));
        idxQ = find(climateDates == currentQDate);
        if idxQ > longWindow
            for qxPix = 1:size(maskData,2)
                for qyPix = 1:size(maskData,1)
                    if maskData(qyPix,qxPix) == 1
                        for k = (longWindow-1):-1:0
                            queryClimate(k+1,:) = cellfun(@(x) x(qyPix,qxPix), climateMaps(idxQ-k,:));
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
                        targetDistance  = nan; %nan(totLDates,1);
                        addVarsDistance = nan; %nan(totLDates,1);
                        %climateDistance = nan(1,numel(climateVarsNames)); %nan(totLDates,2);
                        distancePix     = nan(totLDates,2);
                        %climateDistAll  = cell(longWindow,numel(climateVars));
                        % Display progress - only for serial computing
                        %fprintf(1,'    Progress for current query date: %3.0f%%\n',progress);
                        for ld = 1:totLDates
                            learningClimate = nan(longWindow, numel(climateVars));
                            currentLDate    = learningDatesDate(ld);
                            dayOfYearL      = day(datetime(currentLDate,'ConvertFrom','yyyyMMdd'),'dayofyear');
                            idxL             = find(climateDates == currentLDate);
                            %disp(['    Computing distance to day ' num2str(l) '/' num2str(totLDates) ' (' num2str(currentLDate) ')'])
                            if ismember(dayOfYearL,rangeQ) % if learning date is not within 3 months of the query date, it is skipped
                                if idxL >= longWindow % skips learning dates that are in the longWindow
                                    % Learning dates climate
                                    for k = (longWindow-1):-1:0
                                        learningClimate(k+1,:) = cellfun(@(x) x(qyPix,qxPix), climateMaps(idxL-k,:));
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
                                    if ~(isempty(cell2mat(queryDatesData(qd,:))) || unique(isnan(cell2mat(queryDatesData(qd,:))))) %&& sum(sum(cell2mat(queryDatesData(qd,:))))~=0
                                        if metricKNN == 1 % RMSE
                                            targetDistance(ld) = cellfun(@(x, y) sqrt(mean((x - y).^2, 'all', 'omitnan')), ...
                                                queryDatesData(qd,:), learningDatesData(ld,:), 'UniformOutput', false); % RMSE
                                        elseif metricKNN == 2 % MAE
                                            targetDistance(ld) = cellfun(@(x, y) mean(abs(x - y), 'all', 'omitnan'), ...
                                                queryDatesData(qd,:), learningDatesData(ld,:), 'UniformOutput', false); % MAE
                                        elseif metricKNN == 3 % Manhattan
                                            targetDistance(ld) = cellfun(@(x, y) sum(abs(x - y), 'all', 'omitnan'), ...
                                                queryDatesData(qd,:), learningDatesData(ld,:), 'UniformOutput', false); % Manhattan
                                        elseif metricKNN == 4 % Euclidean
                                            targetDistance(ld) = cellfun(@(x, y) sqrt(sum((x - y).^2, 'all', 'omitnan')), ...
                                                queryDatesData(qd,:), learningDatesData(ld,:), 'UniformOutput', false); % Euclidean
                                        else
                                            error('Bad metricKNN parameter')
                                        end
                                        targetDistance(ld) = sum(cell2mat(targetDistance(ld)),1,'omitnan');
                                        if optimPrep == false
                                            targetDistance(ld) = targetDistance(ld).*weightsTarget;
                                            targetDistance(ld) = sum(cell2mat(targetDistance(ld)),2,'omitnan');
                                        end
                                    else
                                        targetDistance = 0;
                                    end

                                    % Additional variable comparison
                                    % 1 distance
                                    if ~isempty(addVars) && ~isempty(addVarsData)
                                        if ~isempty(addVarsData(qd,:))
                                            if metricKNN == 1 % RMSE
                                                addVarsDistance(ld) = cellfun(@(x, y) sqrt(mean((x - y).^2, 'all', 'omitnan')), ...
                                                    queryAddVars, learningAddVars, 'UniformOutput', false); % RMSE
                                            elseif metricKNN == 2 % MAE
                                                addVarsDistance(ld) = cellfun(@(x, y) mean(abs(x - y), 'all', 'omitnan'), ...
                                                    queryAddVars, learningAddVars, 'UniformOutput', false); % MAE
                                            elseif metricKNN == 3 % Manhattan
                                                addVarsDistance(ld) = cellfun(@(x, y) sum(abs(x - y), 'all', 'omitnan'), ...
                                                    queryAddVars, learningAddVars, 'UniformOutput', false); % Manhattan
                                            elseif metricKNN == 4 % Euclidean
                                                addVarsDistance(ld) = cellfun(@(x, y) sqrt(sum((x - y).^2, 'all', 'omitnan')), ...
                                                    queryAddVars, learningAddVars, 'UniformOutput', false); % Euclidean
                                            else
                                                error('Bad metricKNN parameter')
                                            end
                                            addVarsDistance(ld) = sum(cell2mat(addVarsDistance(ld)),1,'omitnan');
                                            if optimPrep == false
                                                if numel(addVars) == 1
                                                    addVarsDistance(ld) = addVarsDistance(ld) .* cell2mat(weightsAddVars);
                                                else
                                                    addVarsDistance(ld) = num2cell(cell2mat(addVarsDistance(ld)) .* cell2mat(weightsAddVars));
                                                end
                                            end
                                        end
                                    else
                                        addVarsDistance = 0;
                                    end

                                    % Climate distance
%                                     climateDistance = 1 - corr(learningClimate,queryClimate);
                                    if metricKNN == 1 % RMSE
                                        climateDistance = sqrt(mean(learningClimate - queryClimate).^2);
                                    elseif metricKNN == 2 % MAE
                                        climateDistance = mean(abs(learningClimate - queryClimate));
                                    elseif metricKNN == 3 % Manhattan
                                        climateDistance = sum(abs(learningClimate - queryClimate));
                                    elseif metricKNN == 4 % Euclidean
                                        climateDistance = sqrt(sum((learningClimate - queryClimate).^2));
                                    else
                                        error('Bad metricKNN parameter')
                                    end

                                    distancePix(ld,1) = currentLDate;
                                    distancePix(ld,2) = sum(climateDistance) + targetDistance + addVarsDistance;
                                    if distancePix(ld,2) == 0
                                        distancePix(ld,2) = distancePix(ld,2) + eps; % VERY STRANGE, CHECK WHY 0 SOMETIMES <--------------------------------------------------------------------------------------------------------------
                                    end
                                else
                                    continue
                                end
                            else
                                % If not enough climate days available, skip until loop reaches longWindow
                                %warning(['Climate data available is shorter than longWindow, ' num2str(currentLDate) ' skipped.'])
                                continue
                            end
                        end
                        distancePix = sortrows(distancePix(~isnan(distancePix(:,2)),:),2);
                        distancePix = distancePix(1:nbImages,:);
                        distMapQd{qyPix,qxPix} = distancePix;
                    else
                        % If learning date not in query date range, skip it
                        %disp(['    Learning day ',num2str(currentLDate),' not in query date range, skipped'])
                        continue
                    end
                end
            end
        else
            fprintf('\n')
            disp('    Not enough learning dates. Query date skipped...')
            fprintf('\b')
            continue
        end

        distanceAll(:,:,qd) = distMapQd;

        %         % Learning dates distance: 1 date, 2 distance
        %         distance        = climateDistance(~cellfun('isempty',climateDistance(:,1)),:);
        %         targetDistance  = targetDistance(~cellfun('isempty',targetDistance),:);
        %         addVarsDistance = addVarsDistance(~cellfun('isempty',addVarsDistance),:);
        %         if optimPrep == false
        %             distancesSort   = sortrows(distance,2); % Sort rows in ascending order according to column 2
        %             distancesBest   = distancesSort(1:nbImages,1);
        %             distSorted      = distancesSort(1:nbImages,2);
        %             sortedDates{qd} = currentQDate;
        %             sortedData{qd}  = cell2mat(distancesBest);
        %             sortedDist{qd}  = cell2mat(distSorted);
        %         else
        %             distancesBest     = distance(:,1);
        %             distSorted        = distance(:,2);
        %             sortedDates{qd}   = currentQDate;
        %             sortedData{qd}    = distancesBest;
        %             sortedTarget{qd}  = targetDistance;
        %             sortedAddVars{qd} = addVarsDistance;
        %             sortedDist{qd}    = distSorted;
        %         end
    end

    %     if optimPrep == false
    %         sortedDatesAll = [sortedDates sortedData sortedDist];
    %         sortedDates    = sortedDatesAll;
    %     else
    %         sortedDatesAll = [sortedDates sortedData sortedTarget sortedAddVars sortedDist];
    %         sortedDates    = sortedDatesAll;
    %     end

    sortedDates.data = distanceAll;
    sortedDates.date = queryDatesDate;

    % Shut down parallel pool
    poolobj = gcp('nocreate');
    delete(poolobj);
else
    %% serial computing
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
        distMapQd = cell(size(maskData));

        queryClimate = nan(longWindow, numel(climateVars));
        idxQ = find(climateDates == currentQDate);
        if idxQ > longWindow
            % Display progress - only for serial computing
            cPix = 0;
            progress = 0;
            fprintf(1,'\n    Progress for current query date: %3.0f%%\n',progress);
            for qxPix = 1:size(maskData,2)
                for qyPix = 1:size(maskData,1)
                    if maskData(qyPix,qxPix) == 1
                        for k = (longWindow-1):-1:0
                            queryClimate(k+1,:) = cellfun(@(x) x(qyPix,qxPix), climateMaps(idxQ-k,:));
                        end

                        % Extract the additional data for the current query date
                        if ~isempty(additionalVars)
                            queryAddVars = cell(1, numel(addVars));
                            idxQ = find(addVarsDates == currentQDate);
                            queryAddVars(1,:) = addVarsData(idxQ,:);
                        else
                            queryAddVars = [];
                        end

                        % Compute the distances between the query climate and the climate for each learning date
                        %targetDistance  = cell(totLDates,1);
                        addVarsDistance = cell(totLDates,1);
                        %climateDistance = nan(totLDates,2);
                        distancePix     = nan(totLDates,2);
                        %climateDistAll  = cell(longWindow,numel(climateVars));
                        % Display progress - only for serial computing
                        %progress = 0;
                        %fprintf(1,'\n    Progress for current query date: %3.0f%%\n',progress);
                        for ld = 1:totLDates
                            learningClimate = nan(longWindow, numel(climateVars));
                            currentLDate    = learningDatesDate(ld);
                            dayOfYearL      = day(datetime(currentLDate,'ConvertFrom','yyyyMMdd'),'dayofyear');
                            idxL            = find(climateDates == currentLDate);
                            %disp(['    Computing distance to day ' num2str(l) '/' num2str(totLDates) ' (' num2str(currentLDate) ')'])
                            if ismember(dayOfYearL,rangeQ) % if learning date is not within 3 months of the query date, it is skipped
                                %disp(['    Processing learning day ', num2str(currentLDate)])
                                if idxL >= longWindow % skips learning dates that are in the longWindow
                                    % Learning dates climate
                                    for k = (longWindow-1):-1:0
                                        learningClimate(k+1,:) = cellfun(@(x) x(qyPix,qxPix), climateMaps(idxL-k,:));
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
                                    if ~(isempty(cell2mat(queryDatesData(qd,:))) || unique(isnan(cell2mat(queryDatesData(qd,:))))) %&& sum(sum(cell2mat(queryDatesData(qd,:))))~=0
                                        if metricKNN == 1 % RMSE
                                            targetDistance = cellfun(@(x, y) sqrt(mean((x - y).^2, 'all', 'omitnan')), ...
                                                queryDatesData(qd,:), learningDatesData(ld,:), 'UniformOutput', false); % RMSE
                                        elseif metricKNN == 2 % MAE
                                            targetDistance = cellfun(@(x, y) mean(abs(x - y), 'all', 'omitnan'), ...
                                                queryDatesData(qd,:), learningDatesData(ld,:), 'UniformOutput', false); % MAE
                                        elseif metricKNN == 3 % Manhattan
                                            targetDistance = cellfun(@(x, y) sum(abs(x - y), 'all', 'omitnan'), ...
                                                queryDatesData(qd,:), learningDatesData(ld,:), 'UniformOutput', false); % Manhattan
                                        elseif metricKNN == 4 % Euclidean
                                            targetDistance = cellfun(@(x, y) sqrt(sum((x - y).^2, 'all', 'omitnan')), ...
                                                queryDatesData(qd,:), learningDatesData(ld,:), 'UniformOutput', false); % Euclidean
                                        else
                                            error('Bad metricKNN parameter')
                                        end
                                        targetDistance = sum(cell2mat(targetDistance),1,'omitnan');
                                        if optimPrep == false
                                            targetDistance = targetDistance.*weightsTarget;
                                            targetDistance = sum(cell2mat(targetDistance(ld)),2,'omitnan');
                                        end
                                    else
                                        targetDistance = 0;
                                    end

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
                                        addVarsDistance = 0;
                                    end

                                    % Climate distance
%                                     for climVar = 1:numel(climateVars)
%                                         climateDistance(:,climVar) = corr(learningClimate(:,climVar),queryClimate(:,climVar));
%                                     end
                                    if metricKNN == 1 % RMSE
                                        climateDistance = sqrt(mean(learningClimate - queryClimate).^2);
                                    elseif metricKNN == 2 % MAE
                                        climateDistance = mean(abs(learningClimate - queryClimate));
                                    elseif metricKNN == 3 % Manhattan
                                        climateDistance = sum(abs(learningClimate - queryClimate));
                                    elseif metricKNN == 4 % Euclidean
                                        climateDistance = sqrt(sum((learningClimate - queryClimate).^2));
                                    else
                                        error('Bad metricKNN parameter')
                                    end
                                    distancePix(ld,1) = currentLDate;
                                    distancePix(ld,2) = sum(climateDistance) + targetDistance + addVarsDistance;
                                    if distancePix(ld,2) == 0
                                        distancePix(ld,2) = distancePix(ld,2) + eps; % VERY STRANGE, CHECK WHY 0 SOMETIMES <--------------------------------------------------------------------------------------------------------------
                                    end
                                    %                                     climateDistance{ld,1} = currentLDate;
                                    %                                     if shortWindow > 0
                                    %                                         climateDistance{ld,2}(1,:) = sum(cell2mat(climateDistAll(1:shortWindow,:)),1,'omitnan');
                                    %                                     else
                                    %                                         climateDistance{ld,2}(1,:) = single(zeros(1,size(climateDistAll,2)));
                                    %                                     end
                                    %                                     climateDistance{ld,2}(2,:) = sum(cell2mat(climateDistAll(shortWindow+1:end,:)),1,'omitnan');
                                    %                                     % Assign weights to corresponding index
                                    %                                     if optimPrep == false
                                    %                                         climateDistance{ld,2}(1,:) = climateDistance{ld,2}(1,:) .* cell2mat(weightsShort);
                                    %                                         climateDistance{ld,2}(2,:) = climateDistance{ld,2}(2,:) .* cell2mat(weightsLong);
                                    %                                         climateDistance{ld,2} = sum(climateDistance{ld,2},1,'omitnan');
                                    %                                         climateDistance{ld,2} = sum(climateDistance{ld,2},2,'omitnan')+targetDistance{ld}+addVarsDistance{ld};
                                    %                                     end
                                else
                                    % If not enough climate days available, skip until loop reaches longWindow
                                    %warning(['Climate data available is shorter than longWindow, ' num2str(currentLDate) ' skipped.'])
                                    continue
                                end
                            else
                                % If learning date not in query date range, skip it
                                % Display computation progress - only for serial computing
                                %progress = (100*(ld/totLDates));
                                %fprintf(1,'\b\b\b\b%3.0f%%',progress);
                                %disp(['    Learning day ',num2str(currentLDate),' not in query date range, skipped'])
                                continue
                            end
                        end
                        distancePix = sortrows(distancePix(~isnan(distancePix(:,2)),:),2);
                        distancePix = distancePix(1:nbImages,:);
                        distMapQd{qyPix,qxPix} = distancePix;
                    else
                        continue
                    end
                    cPix = cPix + 1;
                end
                % Display computation progress - only for serial computing
                %progress = (100*(ld/totLDates));
                progress = 100*(cPix/totGPixels);
                fprintf(1,'\b\b\b\b%3.0f%%',progress);
            end
        else
            fprintf('\n')
            disp('    Not enough learning dates. Query date skipped...')
            fprintf('\b')
            continue
        end

        distanceAll(:,:,qd) = distMapQd;

        %         % Learning dates distance: 1 date, 2 distance
        %         distance        = climateDistance(~cellfun('isempty',climateDistance(:,1)),:);
        %         targetDistance  = targetDistance(~cellfun('isempty',targetDistance),:);
        %         addVarsDistance = addVarsDistance(~cellfun('isempty',addVarsDistance),:);
        %         if optimPrep == false
        %             distancesSort   = sortrows(distance,2); % Sort rows in ascending order according to column 2
        %             distancesBest   = distancesSort(1:nbImages,1);
        %             distSorted      = distancesSort(1:nbImages,2);
        %             sortedDates{qd} = currentQDate;
        %             sortedData{qd}  = cell2mat(distancesBest);
        %             sortedDist{qd}  = cell2mat(distSorted);
        %         else
        %             distancesBest     = distance(:,1);
        %             distSorted        = distance(:,2);
        %             sortedDates{qd}   = currentQDate;
        %             sortedData{qd}    = distancesBest;
        %             sortedTarget{qd}  = targetDistance;
        %             sortedAddVars{qd} = addVarsDistance;
        %             sortedDist{qd}    = distSorted;
        %         end

        % Display progression - for parallel computing
        %progress = (100*(l/totLDates));
        %fprintf(1,'\b\b\b\b%3.0f%%',progress);

        %toc
    end

    %     if optimPrep == false
    %         sortedDatesAll = [sortedDates sortedData sortedDist];
    %         sortedDates    = sortedDatesAll;
    %     else
    %         sortedDatesAll = [sortedDates sortedData sortedTarget sortedAddVars sortedDist];
    %         sortedDates    = sortedDatesAll;
    %     end
    sortedDates.data = distanceAll;
    sortedDates.date = queryDatesDate;

    fprintf('\n')
end

%sortedDates = sortedDates(~cellfun('isempty',sortedDates(:,1)),:);

if optimPrep == false
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
