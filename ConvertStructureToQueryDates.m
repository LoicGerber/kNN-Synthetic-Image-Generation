function [queryDates,learningDates,refValidation] = ConvertStructureToQueryDates(var,QdateStart,QdateEnd,learningDates,climateData,longWindow,GeoRef,validationPrep,optimPrep,outputTime,inputDir,outputDir)

%
%
%
% REDO DOCUMENTATION
%
%
%

% Query dates - variable to be generated
disp("  Processing '" + var + "' for queryDates...")

datesAll = climateData.date;
learningDatesDate = learningDates.date;

%imgLength = size(learningDates{1,2}{1,1},1);
%imgWidth  = size(learningDates{1,2}{1,1},2);

% Query dates and adapt Learning dates if Validation ON
if validationPrep == false && optimPrep == false % VALIDATION OFF
    % Query dates are all dates in query window, without dates in Learning dates
    if outputTime == 1 % daily
        % Select the dates that are not in learningDates
        [r,~] = find(datesAll>=QdateStart & datesAll<=QdateEnd);
        queryDates = datesAll(r);
        queryDates = setdiff(queryDates, learningDatesDate);
        if isempty(queryDates)
            error('Query dates match with learning dates, nothing to generate')
        end
    elseif outputTime == 2 % monthly
        [r,~]      = find(datesAll>=QdateStart & datesAll<=QdateEnd);
        queryDates = datesAll(r);
        % Convert dailyDates to a matrix of year, month, and day components
        dateVec    = datevec(datetime(datesAll,'ConvertFrom','yyyyMMdd'));
        % Find the indices of the dates where the day component is the last day of the month
        lastDays   = find(dateVec(:,3) == eomday(dateVec(:,1), dateVec(:,2)));
        % Select the dates that are not in learningDates
        queryDates = setdiff(queryDates(lastDays), learningDatesDate);
        if isempty(queryDates)
            error('Query dates match with learning dates, nothing to generate')
        end
    else
        error('Invalid outputTime value')
    end
    refValidation = [];
elseif validationPrep == true || optimPrep == true % validation or optimPrep ON
    % Query dates are all dates in query window, replacing dates in Learning dates
    if ~exist(outputDir,'dir')
        mkdir(outputDir)
    end
    if ~exist(fullfile(outputDir,'referenceImages'),'dir')
        mkdir(fullfile(outputDir,'referenceImages'))
    end
    delete(fullfile(outputDir,'referenceImages','*'));
    if outputTime == 1 % daily
        % Select the dates that are not in learningDates
        [r,~] = find(datesAll>=QdateStart & datesAll<=QdateEnd);
        queryDates = datesAll(r);
    elseif outputTime == 2 % monthly
        [r,~] = find(datesAll>=QdateStart & datesAll<=QdateEnd);
        queryDates = datesAll(r);
        % Convert dailyDates to a matrix of year, month, and day components
        dateVec    = datevec(datetime(queryDates,'ConvertFrom','yyyyMMdd'));
        % Find the indices of the dates where the day component is the last day of the month
        lastDays   = find(dateVec(:,3) == eomday(dateVec(:,1), dateVec(:,2)));
        % Select the dates that are not in learningDates
        queryDates = queryDates(lastDays);
    else
        error('Invalid outputTime value')
    end
    % Define learningDates as itself minus the query dates
    ismem = ismember(learningDatesDate, queryDates);
    learningDataValidation  = learningDates(~ismem,:);
    learningDatesValidation = learningDatesDate(ismem);
    referenceValidation     = table2cell(learningDates(ismem,var));
    imagesRefValidation     = nan(size(referenceValidation{1,1},1),size(referenceValidation{1,1},2),size(learningDatesValidation,1));
    % Display progress
%     progress = 0;
%     fprintf(1,'  Reference images for validation download progress: %3.0f%%\n',progress);
%     if isempty(GeoRef)
%         for i = 1:size(learningDatesValidation,1)
%             %disp(['  Saving ',num2str(learningDatesValidation(i)),' reference image for validation...'])
%             t = Tiff(fullfile(outputDir,'referenceImages',strcat(num2str(learningDatesValidation(i)),'.tif')), 'w');
%             tagstruct.ImageLength         = imgLength;
%             tagstruct.ImageWidth          = imgWidth;
%             tagstruct.Compression         = Tiff.Compression.None;
%             tagstruct.SampleFormat        = Tiff.SampleFormat.IEEEFP;
%             tagstruct.Photometric         = Tiff.Photometric.MinIsBlack;
%             tagstruct.BitsPerSample       = 32;
%             tagstruct.SamplesPerPixel     = 1;
%             tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
%             t.setTag(tagstruct);
%             t.write(single(referenceValidation{i,1}));
%             t.close();
%             % Display computation progress
%             progress = (100*(i/size(learningDatesValidation,1)));
%             fprintf(1,'\b\b\b\b%3.0f%%',progress);
%             imagesRefValidation(:,:,i) = referenceValidation{i,1};
%         end
%     else
%         for i = 1:size(learningDatesValidation,1)
%             %disp(['  Saving ',num2str(learningDatesValidation(i)),' reference image for validation...'])
%             geotiffwrite(fullfile(outputDir,'referenceImages',strcat(num2str(learningDatesValidation(i)),'.tif')), ...
%                 single(referenceValidation{i,1}),GeoRef,'TiffTags',struct('Compression',Tiff.Compression.None));
%             imagesRefValidation(:,:,i) = referenceValidation{i,1};
%             % Display computation progress
%             progress = (100*(i/size(learningDatesValidation,1)));
%             fprintf(1,'\b\b\b\b%3.0f%%',progress);
%         end
%     end
    for i = 1:size(learningDatesValidation,1)
        imagesRefValidation(:,:,i) = referenceValidation{i,1};
    end
    refValidation.date = learningDatesValidation;
    refValidation.maps = imagesRefValidation;
    disp('Saving refValidation.mat file...')
    save(fullfile(inputDir,'refValidation.mat'), 'refValidation', '-v7.3','-nocompression');
    learningDates = learningDataValidation;
end

% Select closest targetVar index for each Query date
nearestIdx = nan(size(queryDates));
if validationPrep == false && optimPrep == false % validation OFF
    for i = 1:numel(queryDates)
        %[nearest, nearestIdx(i)] = min(abs(learningDatesDate - queryDates(i)));  % find index of closest date
        [nearest, nearestIdx(i)] = min(abs(datetime(learningDatesDate,'ConvertFrom','yyyymmdd') - datetime(queryDates(i),'ConvertFrom','yyyyMMdd')));
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
    % Assign closest targetVar map to each Query date
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
elseif validationPrep == true || optimPrep == true % validation or optimPrep ON
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
fprintf('\n')
disp('  Saving Query dates, may take a while depending on input size...')
save(fullfile(inputDir,'queryDates.mat'), 'queryDates', '-v7.3','-nocompression');

end
