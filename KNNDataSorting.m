function sortedDates = KNNDataSorting(var,vars,addVars,queryDates,learningDates,climateData,additionalVars,shortWindow,longWindow,Weights,nbImages,inputDir)

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
climateData  = table2array(removevars(climateData,'date'));
% the learning dates are ranked based on a criterion that quantifies their distance to a given query date

totQDates = size(queryDates,1);
totLDates = size(learningDates,1);

queryDatesDate = table2array(queryDates(:,1));
queryDatesData = table2array(queryDates(:,2));

learningDatesDate = table2array(learningDates(:,1));
learningDatesData = table2array(learningDates(:,2));

sortedDates = cell(totQDates, 1);
sortedData  = cell(totQDates, 1);

if ~isempty(addVars)
    addVarsDates = table2array(additionalVars(:,'date'));
    addVarsData  = table2array(removevars(additionalVars,'date'));
else
    addVarsDates = [];
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

parfor qd = 1:totQDates
%for i = 1:totQDates
    %tic

    currentQDate = queryDatesDate(qd);
    disp(['  Processing day ' num2str(qd) '/' num2str(totQDates) ' (' num2str(currentQDate) ')'])
    
    % Extract the longWindow climate for the current query date
    queryClimate = cell(longWindow, numel(vars));
    idx = find(climateDates == currentQDate);
    for j = 1:numel(vars)
        kj = 1;
        for k = (longWindow-1):-1:0
            queryClimate(kj,j) = climateData(idx-k,j);
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
    addVarsDistance = cell(totLDates,2);
    distances       = cell(totLDates, 3);
    % Display progress - only for serial computing
    %fprintf(1,'    Progress for current query date: %3.0f%%\n',progress);
    for ld = 1:totLDates
        learningClimate = cell(longWindow, numel(vars));
        currentLDate    = learningDatesDate(ld);
        idx             = find(climateDates == currentLDate);
        %disp(['    Computing distance to day ' num2str(l) '/' num2str(totLDates) ' (' num2str(currentLDate) ')'])
        if idx >= longWindow
            % Learning dates climate
            for j = 1:numel(vars)
                kj = 1;
                for k = (longWindow-1):-1:0
                    learningClimate(kj,j) = climateData(idx-k,j);
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
            targetDistance = cellfun(@minus, queryDatesData(qd), learningDatesData(ld), 'UniformOutput', false);
            targetDistance = cellfun(@abs,targetDistance,'UniformOutput',false);
            targetDistance = cellfun(@(x) mean(x,'all','omitnan'),targetDistance,'UniformOutput',false);
            targetDistance = sum(cellfun(@double,targetDistance),1,'omitnan');
            targetDistance = targetDistance.*weightsTarget;

            % Additional variable comparison
            % 1 distance, 2 std
            if ~isempty(addVars)
                addVarsDistance{ld,1} = cellfun(@minus, queryAddVars, learningAddVars, 'UniformOutput', false);
                addVarsDistance{ld,1} = cellfun(@abs,addVarsDistance{ld,1},'UniformOutput',false);
                addVarsDistance{ld,1} = cellfun(@(x) mean(x,'all','omitnan'),addVarsDistance{ld,1},'UniformOutput',false);
                addVarsDistance{ld,2} = cellfun(@(x) std(x,0,'all','omitnan'),addVarsDistance{ld,1},'UniformOutput',false);
                addVarsDistance{ld,1} = sum(cellfun(@double,addVarsDistance{ld,1}),1,'omitnan');
                addVarsDistance{ld,2} = sum(cellfun(@double,addVarsDistance{ld,2}),1,'omitnan');
                if numel(addVars) == 1
                    addVarsDistance{ld,1} = addVarsDistance{ld,1} .* cell2mat(weightsAddVars);
                    addVarsDistance{ld,2} = addVarsDistance{ld,2} .* cell2mat(weightsAddVars);
                else
                    addVarsDistance{ld,1} = num2cell(cell2mat(addVarsDistance{ld,1}) .* cell2mat(weightsAddVars));
                    addVarsDistance{ld,2} = num2cell(cell2mat(addVarsDistance{ld,2}) .* cell2mat(weightsAddVars));
                end
            else
                addVarsDistance{ld,1} = 0;
                addVarsDistance{ld,2} = 0;
            end
            
            % Climate distance
            % 1 date, 2 distance, 3 std
            distances{ld,1} = currentLDate;
            distances{ld,2} = cellfun(@minus, learningClimate, queryClimate, 'UniformOutput', false);
            distances{ld,2} = cellfun(@abs,distances{ld,2},'UniformOutput',false);
            distances{ld,3} = cellfun(@(x) std(x,0,'all','omitnan'),distances{ld,2},'UniformOutput',false);
            distances{ld,2} = cellfun(@(x) mean(x,'all','omitnan'),distances{ld,2},'UniformOutput',false);
            % Assign weights to corresponding index
            distances{ld,2}(1:shortWindow,:)     = num2cell(cell2mat(distances{ld,2}(1:shortWindow,:))     .* cell2mat(weightsShort));
            distances{ld,2}(shortWindow+1:end,:) = num2cell(cell2mat(distances{ld,2}(shortWindow+1:end,:)) .* cell2mat(weightsLong));
            distances{ld,3}(1:shortWindow,:)     = num2cell(cell2mat(distances{ld,3}(1:shortWindow,:))     .* cell2mat(weightsShort));
            distances{ld,3}(shortWindow+1:end,:) = num2cell(cell2mat(distances{ld,3}(shortWindow+1:end,:)) .* cell2mat(weightsLong));
            distances{ld,2} = sum(cellfun(@double,distances{ld,2}),1,'omitnan');
            distances{ld,3} = sum(cellfun(@double,distances{ld,3}),1,'omitnan');
            % !!! if not all query dates have a close target map, bias towards those without map because distance will always be smaller !!!
            distances{ld,2} = sum(distances{ld,2},2,'omitnan')+targetDistance+addVarsDistance{ld,1};
            distances{ld,3} = sum(distances{ld,3},2,'omitnan')+addVarsDistance{ld,2};
        else
            % If not enough climate days available, skip until loop reaches longWindow
            %warning(['Climate data available is shorter than longWindow, ' num2str(currentLDate) ' skipped.'])
            continue
        end
        % Display computation progress - only for serial computing
        %progress = (100*(l/totLDates));
	    %fprintf(1,'\b\b\b\b%3.0f%%',progress);
    end
    
    distances       = distances(~cellfun('isempty',distances(:,1)),:);
    distancesSort   = sortrows(distances,2);
    distancesBest   = distancesSort(1:nbImages,1);
    sortedDates{qd} = currentQDate;
    sortedData{qd}  = cellfun(@double,distancesBest);
    sortedDist{qd}  = cellfun(@double,distancesSort(1:nbImages,2));
    sortedStd{qd}   = cellfun(@double,distancesSort(1:nbImages,3));
    
    % Display progression - for parallel computing
    %progress = (100*(l/totLDates));
	%fprintf(1,'\b\b\b\b%3.0f%%',progress);

    %toc
end

sortedDatesAll = [sortedDates sortedData]; % sortedDist sortedStd];
sortedDates    = sortedDatesAll;

% Shut down parallel pool
poolobj = gcp('nocreate');
delete(poolobj);

disp('Saving results...')
save(fullfile(inputDir,'KNNSorting.mat'),'sortedDates', '-v7.3','-nocompression'); % Save Ranked Learning Dates per Query Date

disp('KNN sorting done! Exiting function...')

toc

end
