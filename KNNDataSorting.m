function sortedDates = KNNDataSorting(var,vars,addVars,queryDates,learningDates,climateData,additionalVars,shortWindow,longWindow,Weights,nbImages,optimPrep,saveOptimPrep,parallelComputing,inputDir)

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
climateMaps  = table2array(removevars(climateData,'date'));

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

if ~isempty(addVars)
    addVarsDates = table2array(additionalVars(:,'date'));
    addVarsData  = table2array(removevars(additionalVars,'date'));
else
    addVarsDates = [];
    addVarsData  = [];
end

% Assign different weights
idxTarget     = contains(Weights.Properties.VariableNames,var);
weightsTarget = table2array(Weights(:,idxTarget));
idxShort      = contains(Weights.Properties.VariableNames,'Short');
weightsShort  = table2cell(Weights(:,idxShort));
idxLong       = contains(Weights.Properties.VariableNames,'Long');
weightsLong   = table2cell(Weights(:,idxLong));
if ~isempty(addVars)
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
        minRangeQ = dayOfYearQ - 90;
        if minRangeQ <= 0, minRangeQ = 365 + minRangeQ; end
        maxRangeQ = dayOfYearQ + 90;
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
        queryClimate = cell(longWindow, numel(vars));
        idx = find(climateDates == currentQDate);
        if idx > longWindow
            for j = 1:numel(vars)
                kj = 1;
                for k = (longWindow-1):-1:0
                    queryClimate(kj,j) = climateMaps(idx-k,j);
                    kj = kj+1;
                end
            end

            % Extract the additional data for the current query date
            if ~isempty(addVars)
                queryAddVars = cell(1, numel(addVars));
                idx = find(addVarsDates == currentQDate);
                for j = 1:numel(addVars)
                    queryAddVars(1,j) = addVarsData(idx,j);
                end
            else
                queryAddVars = [];
            end

            % Compute the distances between the query climate and the climate for each learning date
            targetDistance  = cell(totLDates,1);
            addVarsDistance = cell(totLDates,1);
            climateDistance = cell(totLDates,2);
            % Display progress - only for serial computing
            %fprintf(1,'    Progress for current query date: %3.0f%%\n',progress);
            for ld = 1:totLDates
                learningClimate = cell(longWindow, numel(vars));
                currentLDate    = learningDatesDate(ld);
                dayOfYearL      = day(datetime(currentLDate,'ConvertFrom','yyyyMMdd'),'dayofyear');
                idx             = find(climateDates == currentLDate);
                %disp(['    Computing distance to day ' num2str(l) '/' num2str(totLDates) ' (' num2str(currentLDate) ')'])
                if ismember(dayOfYearL,rangeQ) % if learning date is not within 3 months of the query date, it is skipped
                    if idx >= longWindow % skips learning dates that are in the longWindow
                        % Learning dates climate
                        for j = 1:numel(vars)
                            kj = 1;
                            for k = (longWindow-1):-1:0
                                learningClimate(kj,j) = climateMaps(idx-k,j);
                                kj = kj+1;
                            end
                        end

                        % Extract the additional data for the current query date
                        if ~isempty(addVars)
                            learningAddVars = cell(1, numel(addVars));
                            idx = find(addVarsDates == currentLDate);
                            for j = 1:numel(addVars)
                                learningAddVars(1,j) = addVarsData(idx,j);
                            end
                        else
                            learningAddVars = [];
                        end

                        % Target variable comparison
                        if ~isnan(sum(sum(cell2mat(queryDatesData(qd,:))))) %&& ~isnan(sum(sum(cell2mat(queryDatesData(qd,2)))))
                            targetDistance{ld} = cellfun(@(x, y) mean(abs(x - y), 'all', 'omitnan'), ...
                                queryDatesData(qd,:), learningDatesData(ld,:), 'UniformOutput', false);
                            targetDistance{ld} = sum(cell2mat(targetDistance{ld}),1,'omitnan');
                            if optimPrep == false
                                targetDistance{ld} = targetDistance{ld}.*weightsTarget;
                                targetDistance{ld} = sum(cell2mat(targetDistance(ld)),2,'omitnan');
                            end
                        else
                            targetDistance{ld} = 0;
                        end

                        % Additional variable comparison
                        % 1 distance
                        if ~isempty(addVars)
                            addVarsDistance{ld,1} = cellfun(@(x, y) mean(abs(x - y), 'all', 'omitnan'), ...
                                queryAddVars, learningAddVars, 'UniformOutput', false);
                            addVarsDistance{ld,1} = sum(cell2mat(addVarsDistance{ld,1}),1,'omitnan');
                            if optimPrep == false
                                if numel(addVars) == 1
                                    addVarsDistance{ld,1} = addVarsDistance{ld,1} .* cell2mat(weightsAddVars);
                                else
                                    addVarsDistance{ld,1} = num2cell(cell2mat(addVarsDistance{ld,1}) .* cell2mat(weightsAddVars));
                                end
                            end
                        else
                            addVarsDistance{ld,1} = 0;
                        end

                        % Climate distance
                        % 1 date, 2 distance
                        climateDistAll = cell(ld,1);
                        climateDistAll{ld,1} = cellfun(@(x, y) mean(abs(x - y), 'all', 'omitnan'), ...
                            learningClimate, queryClimate, 'UniformOutput', false); % mean absolute error
                        climateDistance{ld,1} = currentLDate;
                        climateDistance{ld,2}(1,:) = sum(cell2mat(climateDistAll{ld,1}(1:shortWindow,:)),1,'omitnan');
                        climateDistance{ld,2}(2,:) = sum(cell2mat(climateDistAll{ld,1}(shortWindow+1:end,:)),1,'omitnan');
                        % Assign weights to corresponding index
                        if optimPrep == false
                            climateDistance{ld,2}(1,:) = climateDistance{ld,2}(1,:) .* cell2mat(weightsShort);
                            climateDistance{ld,2}(2,:) = climateDistance{ld,2}(2,:) .* cell2mat(weightsLong);
                            climateDistance{ld,2} = sum(climateDistance{ld,2},1,'omitnan');
                            climateDistance{ld,2} = sum(climateDistance{ld,2},2,'omitnan')+targetDistance{ld,1}+addVarsDistance{ld,1};
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
        targetDistance  = targetDistance(~cellfun('isempty',targetDistance),:);
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
            sortedTarget{qd}  = targetDistance;
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
        minRangeQ = dayOfYearQ - 90;
        if minRangeQ <= 0, minRangeQ = 365 + minRangeQ; end
        maxRangeQ = dayOfYearQ + 90;
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
        queryClimate = cell(longWindow, numel(vars));
        idx = find(climateDates == currentQDate);
        if idx > longWindow
            for j = 1:numel(vars)
                kj = 1;
                for k = (longWindow-1):-1:0
                    queryClimate(kj,j) = climateMaps(idx-k,j);
                    %queryClimate(kj,j) = climateMaps.(varClimate)(:,:,idx-k);
                    kj = kj+1;
                end
            end

            % Extract the additional data for the current query date
            if ~isempty(addVars)
                queryAddVars = cell(1, numel(addVars));
                idx = find(addVarsDates == currentQDate);
                for j = 1:numel(addVars)
                    queryAddVars(1,j) = addVarsData(idx,j);
                end
            else
                queryAddVars = [];
            end

            % Compute the distances between the query climate and the climate for each learning date
            targetDistance  = cell(totLDates,1);
            addVarsDistance = cell(totLDates,1);
            climateDistance = cell(totLDates,2);
            climateDistAll  = cell(1,3);
            % Display progress - only for serial computing
            progress = 0;
            fprintf(1,'\n    Progress for current query date: %3.0f%%\n',progress);
            for ld = 1:totLDates
                learningClimate = cell(longWindow, numel(vars));
                currentLDate    = learningDatesDate(ld);
                dayOfYearL      = day(datetime(currentLDate,'ConvertFrom','yyyyMMdd'),'dayofyear');
                idx             = find(climateDates == currentLDate);
                %disp(['    Computing distance to day ' num2str(l) '/' num2str(totLDates) ' (' num2str(currentLDate) ')'])
                if ismember(dayOfYearL,rangeQ) % if learning date is not within 3 months of the query date, it is skipped
                    %disp(['    Processing learning day ', num2str(currentLDate)])
                    if idx >= longWindow % skips learning dates that are in the longWindow
                        % Learning dates climate
                        for j = 1:numel(vars)
                            kj = 1;
                            for k = (longWindow-1):-1:0
                                learningClimate(kj,j) = climateMaps(idx-k,j);
                                kj = kj+1;
                            end
                        end

                        % Extract the additional data for the current query date
                        if ~isempty(addVars)
                            learningAddVars = cell(1, numel(addVars));
                            idx = find(addVarsDates == currentLDate);
                            for j = 1:numel(addVars)
                                learningAddVars(1,j) = addVarsData(idx,j);
                            end
                        else
                            learningAddVars = [];
                        end

                        % Target variable comparison
                        if ~isnan(sum(sum(cell2mat(queryDatesData(qd,1))))) && ~isnan(sum(sum(cell2mat(queryDatesData(qd,2)))))
                            targetDistance{ld} = cellfun(@(x, y) mean(abs(x - y), 'all', 'omitnan'), ...
                                queryDatesData(qd,:), learningDatesData(ld,:), 'UniformOutput', false);
                            targetDistance{ld} = sum(cell2mat(targetDistance{ld}),1,'omitnan');
                            if optimPrep == false
                                targetDistance{ld} = targetDistance{ld}.*weightsTarget;
                                targetDistance{ld} = sum(cell2mat(targetDistance(ld)),2,'omitnan');
                            end
                        else
                            targetDistance{ld} = 0;
                        end

                        % Additional variable comparison
                        % 1 distance
                        if ~isempty(addVars)
                            addVarsDistance{ld} = cellfun(@(x, y) mean(abs(x - y), 'all', 'omitnan'), ...
                                queryAddVars, learningAddVars, 'UniformOutput', false);
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
%                         manhattan = cellfun(@(x, y) sum(abs(x - y), 'all', 'omitnan'), ...
%                             learningClimate, queryClimate, 'UniformOutput', false); % Manhattan norm
%                         euclidean = cellfun(@(x, y) sqrt(sum((x - y).^2, 'all', 'omitnan')), ...
%                             learningClimate, queryClimate, 'UniformOutput', false); % Euclidean norm
%                         mae       = cellfun(@(x, y) mean(abs(x - y), 'all', 'omitnan'), ...
%                             learningClimate, queryClimate, 'UniformOutput', false); % mean absolute error
                        climateDistAll{1,1} = cellfun(@(x, y) mean(abs(x - y), 'all', 'omitnan'), ...
                            learningClimate, queryClimate, 'UniformOutput', false); % mean absolute error
                        climateDistance{ld,1} = currentLDate;
                        climateDistance{ld,2}(1,:) = sum(cell2mat(climateDistAll{1,1}(1:shortWindow,:)),1,'omitnan');
                        climateDistance{ld,2}(2,:) = sum(cell2mat(climateDistAll{1,1}(shortWindow+1:end,:)),1,'omitnan');
                        % Assign weights to corresponding index
                        if optimPrep == false
                            climateDistance{ld,2}(1,:) = climateDistance{ld,2}(1,:) .* cell2mat(weightsShort);
                            climateDistance{ld,2}(2,:) = climateDistance{ld,2}(2,:) .* cell2mat(weightsLong);
                            climateDistance{ld,2} = sum(climateDistance{ld,2},1,'omitnan');
                            climateDistance{ld,2} = sum(climateDistance{ld,2},2,'omitnan')+targetDistance{ld}+addVarsDistance{ld};
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
        targetDistance  = targetDistance(~cellfun('isempty',targetDistance),:);
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
            sortedTarget{qd}  = targetDistance;
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
