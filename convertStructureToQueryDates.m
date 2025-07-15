function [queryDates, learningDates, refValidation] = convertStructureToQueryDates( ...
    targetVar, targetDim, QdateStart, QdateEnd, learningDates, ...
    climateData, validationPrep, optimPrep, outputTime, ...
    inputDir, saveMats)

%
%
%
% REDO DOCUMENTATION
%
%
%

disp('  Processing queryDates for all target variables...')

targetVarL = lower(targetVar);
datesAll = climateData.date;
learningDatesDate = learningDates.date;

% === Validate inputs ===
if numel(QdateStart) ~= numel(QdateEnd)
    error('QdateStart and QdateEnd must be the same length.');
end

% === Get indices within Qdate ranges ===
r = [];
for i = 1:numel(QdateStart)
    qStart = QdateStart(i);
    qEnd = QdateEnd(i);
    if min(datesAll) > qStart
        error('Climate data first date > Query period start for period %d', i);
    elseif max(datesAll) < qEnd
        error('Climate data last date < Query period end for period %d', i);
    end
    r = [r; find(datesAll >= qStart & datesAll <= qEnd)];
end
r = unique(r);
queryDates = datesAll(r);

% === Apply temporal filtering ===
dateVec = datevec(datetime(queryDates, 'ConvertFrom', 'yyyyMMdd'));
switch outputTime
    case 1  % Daily — no filtering
        % nothing to do
    case 2  % Monthly — keep last day of month
        lastDays = dateVec(:,3) == eomday(dateVec(:,1), dateVec(:,2));
        queryDates = queryDates(lastDays);
    case 3  % Dekadal — 10th, 20th, and end of month
        dekadalDays = ismember(dateVec(:,3), [10, 20]) | ...
                      dateVec(:,3) == eomday(dateVec(:,1), dateVec(:,2));
        queryDates = queryDates(dekadalDays);
    otherwise
        error('Invalid outputTime value');
end

% === Validation OFF ===
if ~validationPrep && ~optimPrep
    queryDates = setdiff(queryDates, learningDatesDate);
    if isempty(queryDates)
        error('Query dates match with learning dates, nothing to generate');
    end
    refValidation = [];

% === Validation ON ===
else
    % Remove queryDates from learningDates
    ismem = ismember(learningDatesDate, queryDates);
    learningDataValidation = learningDates(~ismem,:);
    learningDatesValidation = learningDatesDate(ismem);

    % Build refValidation struct
    refValidation = struct();
    for j = 1:numel(targetVarL)
        refData = table2cell(learningDates(ismem, targetVarL(j)));
        if targetDim ~= 1
            % 2D: assume matrix time series
            nT = numel(refData);
            sz = size(refData{1});
            imagesRefValidation = nan(sz(1), sz(2), nT);
            for i = 1:nT
                imagesRefValidation(:,:,i) = refData{i};
            end
        else
            % 1D time series
            imagesRefValidation = nan(numel(refData), 1);
            for i = 1:numel(refData)
                imagesRefValidation(i) = refData{i};
            end
        end
        refValidation.(targetVarL{j}) = single(imagesRefValidation);
    end
    refValidation.date = learningDatesValidation;
    learningDates = learningDataValidation;

    if saveMats
        disp('  Saving refValidation.mat file...')
        save(fullfile(inputDir, 'refValidation.mat'), 'refValidation', '-v7.3', '-nocompression');
    end
end

% === Save queryDates ===
if saveMats
    disp('  Saving queryDates.mat file...')
    save(fullfile(inputDir, 'queryDates.mat'), 'queryDates', '-v7.3', '-nocompression');
end

end
