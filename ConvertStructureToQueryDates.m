function [queryDates,learningDates] = ConvertStructureToQueryDates(var,QdateStart,QdateEnd,learningDates,climateData,longWindow,validation,outputTime,inputDir,outputDir)

%
%
%
% REDO DOCUMENTATION
%
%
%

tic

% Query dates - variable to be generated
disp('  Processing ' + var + ' for queryDates...')

datesAll = climateData.date;
learningDatesDate = learningDates.date;

imgLength = size(learningDates{1,2}{1,1},1);
imgWidth  = size(learningDates{1,2}{1,1},2);

% Query dates and adapt Learning dates if Validation ON
if validation == 0 % VALIDATION OFF
    % Query dates are all dates in query window, without dates in Learning dates
    if outputTime == 1 % daily
        % Select the dates that are not in learningDates
        [r,~] = find(datesAll>=QdateStart & datesAll<=QdateEnd);
        queryDates = datesAll(r);
        queryDates = setdiff(queryDates, learningDatesDate);
    elseif outputTime == 2 % monthly
        % Convert dailyDates to a matrix of year, month, and day components
        dateVec    = datevec(datetime(datesAll,'ConvertFrom','yyyyMMdd'));
        % Find the indices of the dates where the day component is the last day of the month
        lastDays   = find(dateVec(:,3) == eomday(dateVec(:,1), dateVec(:,2)));
        % Select the dates that are not in learningDates
        [r,~]      = find(datesAll>=QdateStart & datesAll<=QdateEnd);
        queryDates = datesAll(r);
        queryDates = setdiff(queryDates(lastDays), learningDatesDate);
    else
        error('Invalid outputTime value')
    end
elseif validation == 1 % VALIDATION ON
    % Query dates are all dates in query window, replacing dates in Learning dates
    if ~exist(outputDir,'dir')
        mkdir(outputDir)
    end
    if ~exist(fullfile(outputDir,'referenceImages'),'dir')
        mkdir(fullfile(outputDir,'referenceImages'))
    end
    if outputTime == 1 % daily
        % Select the dates that are not in learningDates
        [r,~] = find(datesAll>=QdateStart & datesAll<=QdateEnd);
        queryDates = datesAll(r);
        % Define learningDates as itself minus the query dates
        ismem = ismember(learningDatesDate, queryDates);
        learningDataValidation  = learningDates(~ismem,:);
        learningDatesValidation = learningDatesDate(ismem);
        referenceValidation     = table2cell(learningDates(ismem,var));
        for i = 1:size(learningDatesValidation,1)
            disp(['  Saving ',num2str(learningDatesValidation(i)),' reference image for validation...'])
            t = Tiff(fullfile(outputDir,'referenceImages',strcat(num2str(learningDatesValidation(i)),'.tif')), 'w');
            tagstruct.ImageLength         = imgLength;
            tagstruct.ImageWidth          = imgWidth;
            tagstruct.Compression         = Tiff.Compression.None;
            tagstruct.SampleFormat        = Tiff.SampleFormat.IEEEFP;
            tagstruct.Photometric         = Tiff.Photometric.MinIsBlack;
            tagstruct.BitsPerSample       = 32;
            tagstruct.SamplesPerPixel     = 1;
            tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
            t.setTag(tagstruct);
            t.write(single(referenceValidation{i,1}));
            t.close();
        end
        learningDates = learningDataValidation;
    elseif outputTime == 2 % monthly
        [r,~] = find(datesAll>=QdateStart & datesAll<=QdateEnd);
        queryDates = datesAll(r);
        % Convert dailyDates to a matrix of year, month, and day components
        dateVec    = datevec(datetime(queryDates,'ConvertFrom','yyyyMMdd'));
        % Find the indices of the dates where the day component is the last day of the month
        lastDays   = find(dateVec(:,3) == eomday(dateVec(:,1), dateVec(:,2)));
        % Select the dates that are not in learningDates
        queryDates = queryDates(lastDays);
        % Define learningDates as itself minus the query dates
        ismem      = ismember(learningDatesDate, queryDates);
        learningDataValidation  = learningDates(~ismem,:);
        learningDatesValidation = learningDatesDate(ismem);
        referenceValidation     = table2cell(learningDates(ismem,var));
        for i = 1:size(learningDatesValidation,1)
            disp(['  Saving ',num2str(learningDatesValidation(i)),' reference image for validation...'])
            t = Tiff(fullfile(outputDir,'referenceImages',strcat(num2str(learningDatesValidation(i)),'.tif')), 'w');
            tagstruct.ImageLength         = imgLength;
            tagstruct.ImageWidth          = imgWidth;
            tagstruct.Compression         = Tiff.Compression.None;
            tagstruct.SampleFormat        = Tiff.SampleFormat.IEEEFP;
            tagstruct.Photometric         = Tiff.Photometric.MinIsBlack;
            tagstruct.BitsPerSample       = 32;
            tagstruct.SamplesPerPixel     = 1;
            tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
            t.setTag(tagstruct);
            t.write(single(referenceValidation{i,1}));
            t.close();
        end
        learningDates = learningDataValidation;
    else
        error('Invalid outputTime value')
    end
end

% Select closest varGen index for each Query date
nearestIdx = nan(size(queryDates));
if validation == 0
    for i = 1:numel(queryDates)
        [nearest, nearestIdx(i)] = min(abs(learningDatesDate - queryDates(i)));  % find index of closest date
        if nearest > longWindow
            nearestIdx(i) = nan;
        end
    end
    % Match the dates
    matchedTargetVarDates = [queryDates, nan(size(queryDates))];
    for i = 1:length(nearestIdx)
        if ~isnan(nearestIdx(i))
            matchedTargetVarDates(i,2) = learningDatesDate(nearestIdx(i));
        end
    end
    % Assign closest varGen map to each Query date
    matchedTargetVarTable = table('Size',size(matchedTargetVarDates), 'VariableTypes',{'double', 'cell'});
    targetVarData = learningDates.(var);
    % Loop through the matched dates
    for i = 1:size(matchedTargetVarDates, 1)
        % Get the date to match
        matchDate = matchedTargetVarDates(i, 2);
        % If a match was found, add the date and data to the output table
        if ~isnan(matchDate)
            matchedTargetVarTable{i, 1} = queryDates(i);
            matchedTargetVarTable{i, 2} = targetVarData(nearestIdx(i));
        else
            matchedTargetVarTable{i, 1} = queryDates(i);
            matchedTargetVarTable{i, 2} = {nan(size(targetVarData{1,1}))};
        end
    end
elseif validation == 1
    matchedTargetVarDates = [queryDates, nan(size(queryDates))];
    matchedTargetVarTable = table('Size',size(matchedTargetVarDates), 'VariableTypes',{'double', 'cell'});
    targetVarData = learningDates.(var);
    % Loop through the matched dates
    for i = 1:size(matchedTargetVarDates, 1)
        % fill the table with NaNs the size of the variable to be generated
        matchedTargetVarTable{i, 1} = queryDates(i);
        matchedTargetVarTable{i, 2} = {nan(size(targetVarData{1,1}))};
    end
end
% Rename the columns
matchedTargetVarTable.Properties.VariableNames = {'Date', convertStringsToChars(var)};
queryDates = matchedTargetVarTable;
disp('  Saving Query dates, may take a while depending on input size...')
save(fullfile(inputDir,'queryDates.mat'), 'queryDates', '-v7.3','-nocompression');

toc

end
