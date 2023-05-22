function sortedDates = KNNDataSorting(var,vars,addVars,queryDates,learningDates,climateData,additionalVars,shortWindow,longWindow,Weights,nbImages,optimisation,parallelComputing,inputDir)

%
%
%
% REDO DOCUMENTATION
%
%
%

tic

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

totQDates = size(queryDates,1);
totLDates = size(learningDates,1);

queryDatesDate = table2array(queryDates(:,1));
queryDatesData = table2array(queryDates(:,2));

learningDatesDate = table2array(learningDates(:,1));
learningDatesData = table2array(learningDates(:,2));

sortedDates   = cell(totQDates, 1);
sortedData    = cell(totQDates, 1);
sortedTarget  = cell(totQDates, 1);
sortedAddVars = cell(totQDates, 1);
sortedDist    = cell(totQDates, 1);
sortedStd     = cell(totQDates, 1);

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

if parallelComputing == 0
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
        addVarsDistance = cell(totLDates,2);
        climateDistance = cell(totLDates,3);
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
                    targetDistance{ld} = cellfun(@minus, queryDatesData(qd), learningDatesData(ld), 'UniformOutput', false);
                    targetDistance{ld} = cellfun(@abs,targetDistance{ld},'UniformOutput',false);
                    targetDistance{ld} = cellfun(@(x) mean(x,'all','omitnan'),targetDistance{ld},'UniformOutput',false);
                    targetDistance{ld} = sum(cellfun(@double,targetDistance{ld}),1,'omitnan');
                    if optimisation == 1
                        targetDistance{ld} = targetDistance{ld}.*weightsTarget;
                    end

                    % Additional variable comparison
                    % 1 distance, 2 std
                    if ~isempty(addVars)
                        addVarsDistance{ld,1} = cellfun(@minus, queryAddVars, learningAddVars, 'UniformOutput', false);
                        addVarsDistance{ld,1} = cellfun(@abs,addVarsDistance{ld,1},'UniformOutput',false);
                        addVarsDistance{ld,1} = cellfun(@(x) mean(x,'all','omitnan'),addVarsDistance{ld,1},'UniformOutput',false);
                        addVarsDistance{ld,2} = cellfun(@(x) std(x,0,'all','omitnan'),addVarsDistance{ld,1},'UniformOutput',false);
                        addVarsDistance{ld,1} = sum(cellfun(@double,addVarsDistance{ld,1}),1,'omitnan');
                        addVarsDistance{ld,2} = sum(cellfun(@double,addVarsDistance{ld,2}),1,'omitnan');
                        if optimisation == 1
                            if numel(addVars) == 1
                                addVarsDistance{ld,1} = addVarsDistance{ld,1} .* cell2mat(weightsAddVars);
                                addVarsDistance{ld,2} = addVarsDistance{ld,2} .* cell2mat(weightsAddVars);
                            else
                                addVarsDistance{ld,1} = num2cell(cell2mat(addVarsDistance{ld,1}) .* cell2mat(weightsAddVars));
                                addVarsDistance{ld,2} = num2cell(cell2mat(addVarsDistance{ld,2}) .* cell2mat(weightsAddVars));
                            end
                        end
                    else
                        addVarsDistance{ld,1} = 0;
                        addVarsDistance{ld,2} = 0;
                    end

                    % Climate distance
                    % 1 date, 2 distance, 3 std
                    climateDistance{ld,1} = currentLDate;
                    climateDistance{ld,2} = cellfun(@minus, learningClimate, queryClimate, 'UniformOutput', false);
                    climateDistance{ld,2} = cellfun(@abs,climateDistance{ld,2},'UniformOutput',false);
                    climateDistance{ld,3} = cellfun(@(x) std(x,0,'all','omitnan'),climateDistance{ld,2},'UniformOutput',false);
                    climateDistance{ld,2} = cellfun(@(x) mean(x,'all','omitnan'),climateDistance{ld,2},'UniformOutput',false);
                    % Assign weights to corresponding index
                    if optimisation == 1
                        climateDistance{ld,2}(1:shortWindow,:)     = num2cell(cell2mat(climateDistance{ld,2}(1:shortWindow,:))     .* cell2mat(weightsShort));
                        climateDistance{ld,2}(shortWindow+1:end,:) = num2cell(cell2mat(climateDistance{ld,2}(shortWindow+1:end,:)) .* cell2mat(weightsLong));
                        climateDistance{ld,3}(1:shortWindow,:)     = num2cell(cell2mat(climateDistance{ld,3}(1:shortWindow,:))     .* cell2mat(weightsShort));
                        climateDistance{ld,3}(shortWindow+1:end,:) = num2cell(cell2mat(climateDistance{ld,3}(shortWindow+1:end,:)) .* cell2mat(weightsLong));
                        climateDistance{ld,2} = sum(cellfun(@double,climateDistance{ld,2}),1,'omitnan');
                        climateDistance{ld,3} = sum(cellfun(@double,climateDistance{ld,3}),1,'omitnan');
                        climateDistance{ld,2} = sum(climateDistance{ld,2},2,'omitnan')+targetDistance{ld,1}+addVarsDistance{ld,1};
                        climateDistance{ld,3} = sum(climateDistance{ld,3},2,'omitnan')+addVarsDistance{ld,2};
                    end
                else
                    % If not enough climate days available, skip until loop reaches longWindow
                    %warning(['Climate data available is shorter than longWindow, ' num2str(currentLDate) ' skipped.'])
                    continue
                end
            else
                % If learning date not in query date range, skip it
                disp('    Learning date not in query date range, skipped')
                continue
            end
            % Display computation progress - only for serial computing
            %progress = (100*(l/totLDates));
    	    %fprintf(1,'\b\b\b\b%3.0f%%',progress);
        end

        % Learning dates distance: 1 date, 2 distance, 3 std
        distance = climateDistance(~cellfun('isempty',climateDistance(:,1)),:);
        if optimisation == 1
            distancesSort   = sortrows(distance,2); % Sort rows in ascending order according to column 2
            distancesBest   = distancesSort(1:nbImages,1);
            distSorted      = distancesSort(1:nbImages,2);
            stdSorted       = distancesSort(1:nbImages,3);
            sortedDates{qd} = currentQDate;
            sortedData{qd}  = cellfun(@double,distancesBest);
            sortedDist{qd}  = cellfun(@double,distSorted);
            sortedStd{qd}   = cellfun(@double,stdSorted);
        else
            distancesBest     = distance(:,1);
            distSorted        = distance(:,2);
            stdSorted         = distance(:,3);
            sortedDates{qd}   = currentQDate;
            sortedData{qd}    = distancesBest;
            sortedTarget{qd}  = targetDistance;
            sortedAddVars{qd} = addVarsDistance(:,1);
            sortedDist{qd}    = distSorted;
            sortedStd{qd}     = stdSorted;
        end

        % Display progression - for parallel computing
        %progress = (100*(l/totLDates));
    	%fprintf(1,'\b\b\b\b%3.0f%%',progress);

        %toc
    end

    if optimisation == 1
        sortedDatesAll = [sortedDates sortedData sortedDist sortedStd];
        sortedDates    = sortedDatesAll;
    else
        sortedDatesAll = [sortedDates sortedData sortedTarget sortedAddVars sortedDist sortedStd];
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
        addVarsDistance = cell(totLDates,2);
        climateDistance = cell(totLDates,3);
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
                    targetDistance{ld} = cellfun(@minus, queryDatesData(qd), learningDatesData(ld), 'UniformOutput', false);
                    targetDistance{ld} = cellfun(@abs,targetDistance{ld},'UniformOutput',false);
                    targetDistance{ld} = cellfun(@(x) mean(x,'all','omitnan'),targetDistance{ld},'UniformOutput',false);
                    targetDistance{ld} = sum(cellfun(@double,targetDistance{ld}),1,'omitnan');
                    if optimisation == 1
                        targetDistance{ld} = targetDistance{ld}.*weightsTarget;
                    end

                    % Additional variable comparison
                    % 1 distance, 2 std
                    if ~isempty(addVars)
                        addVarsDistance{ld,1} = cellfun(@minus, queryAddVars, learningAddVars, 'UniformOutput', false);
                        addVarsDistance{ld,1} = cellfun(@abs,addVarsDistance{ld,1},'UniformOutput',false);
                        addVarsDistance{ld,1} = cellfun(@(x) mean(x,'all','omitnan'),addVarsDistance{ld,1},'UniformOutput',false);
                        addVarsDistance{ld,2} = cellfun(@(x) std(x,0,'all','omitnan'),addVarsDistance{ld,1},'UniformOutput',false);
                        addVarsDistance{ld,1} = sum(cellfun(@double,addVarsDistance{ld,1}),1,'omitnan');
                        addVarsDistance{ld,2} = sum(cellfun(@double,addVarsDistance{ld,2}),1,'omitnan');
                        if optimisation == 1
                            if numel(addVars) == 1
                                addVarsDistance{ld,1} = addVarsDistance{ld,1} .* cell2mat(weightsAddVars);
                                addVarsDistance{ld,2} = addVarsDistance{ld,2} .* cell2mat(weightsAddVars);
                            else
                                addVarsDistance{ld,1} = num2cell(cell2mat(addVarsDistance{ld,1}) .* cell2mat(weightsAddVars));
                                addVarsDistance{ld,2} = num2cell(cell2mat(addVarsDistance{ld,2}) .* cell2mat(weightsAddVars));
                            end
                        end
                    else
                        addVarsDistance{ld,1} = 0;
                        addVarsDistance{ld,2} = 0;
                    end

                    % Climate distance
                    % 1 date, 2 distance, 3 std
                    climateDistance{ld,1} = currentLDate;
                    climateDistance{ld,2} = cellfun(@minus, learningClimate, queryClimate, 'UniformOutput', false);
                    climateDistance{ld,2} = cellfun(@abs,climateDistance{ld,2},'UniformOutput',false);
                    climateDistance{ld,3} = cellfun(@(x) std(x,0,'all','omitnan'),climateDistance{ld,2},'UniformOutput',false);
                    climateDistance{ld,2} = cellfun(@(x) mean(x,'all','omitnan'),climateDistance{ld,2},'UniformOutput',false);
                    % Assign weights to corresponding index
                    if optimisation == 1
                        climateDistance{ld,2}(1:shortWindow,:)     = num2cell(cell2mat(climateDistance{ld,2}(1:shortWindow,:))     .* cell2mat(weightsShort));
                        climateDistance{ld,2}(shortWindow+1:end,:) = num2cell(cell2mat(climateDistance{ld,2}(shortWindow+1:end,:)) .* cell2mat(weightsLong));
                        climateDistance{ld,3}(1:shortWindow,:)     = num2cell(cell2mat(climateDistance{ld,3}(1:shortWindow,:))     .* cell2mat(weightsShort));
                        climateDistance{ld,3}(shortWindow+1:end,:) = num2cell(cell2mat(climateDistance{ld,3}(shortWindow+1:end,:)) .* cell2mat(weightsLong));
                        climateDistance{ld,2} = sum(cellfun(@double,climateDistance{ld,2}),1,'omitnan');
                        climateDistance{ld,3} = sum(cellfun(@double,climateDistance{ld,3}),1,'omitnan');
                        climateDistance{ld,2} = sum(climateDistance{ld,2},2,'omitnan')+targetDistance{ld,1}+addVarsDistance{ld,1};
                        climateDistance{ld,3} = sum(climateDistance{ld,3},2,'omitnan')+addVarsDistance{ld,2};
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
            % Display computation progress - only for serial computing
            progress = (100*(ld/totLDates));
    	    fprintf(1,'\b\b\b\b%3.0f%%',progress);
        end

        % Learning dates distance: 1 date, 2 distance, 3 std
        distance = climateDistance(~cellfun('isempty',climateDistance(:,1)),:);
        if optimisation == 1
            distancesSort   = sortrows(distance,2); % Sort rows in ascending order according to column 2
            distancesBest   = distancesSort(1:nbImages,1);
            distSorted      = distancesSort(1:nbImages,2);
            stdSorted       = distancesSort(1:nbImages,3);
            sortedDates{qd} = currentQDate;
            sortedData{qd}  = cellfun(@double,distancesBest);
            sortedDist{qd}  = cellfun(@double,distSorted);
            sortedStd{qd}   = cellfun(@double,stdSorted);
        else
            distancesBest     = distance(:,1);
            distSorted        = distance(:,2);
            stdSorted         = distance(:,3);
            sortedDates{qd}   = currentQDate;
            sortedData{qd}    = distancesBest;
            sortedTarget{qd}  = targetDistance;
            sortedAddVars{qd} = addVarsDistance(:,1);
            sortedDist{qd}    = distSorted;
            sortedStd{qd}     = stdSorted;
        end

        % Display progression - for parallel computing
        %progress = (100*(l/totLDates));
    	%fprintf(1,'\b\b\b\b%3.0f%%',progress);

        %toc
    end

    if optimisation == 1
        sortedDatesAll = [sortedDates sortedData sortedDist sortedStd];
        sortedDates    = sortedDatesAll;
    else
        sortedDatesAll = [sortedDates sortedData sortedTarget sortedAddVars sortedDist sortedStd];
        sortedDates    = sortedDatesAll;
    end
    fprintf('\n')
end

if optimisation == 1
    disp('Saving KNNSorting.mat file...')
    save(fullfile(inputDir,'KNNSorting.mat'),'sortedDates', '-v7.3','-nocompression'); % Save Ranked Learning Dates per Query Date
else
    disp('Saving KNNDistances.mat file for optimisation...')
    save(fullfile(inputDir,'KNNDistances.mat'),'sortedDates', '-v7.3','-nocompression');
end

disp('KNN sorting done! Exiting function...')

toc

end
