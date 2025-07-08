function synImages = generateSynImgsOptim(targetVar,targetDim,learningDates,sortedDatesOptim,nbImages,generationType)

%
%
%
% REDO DOCUMENTATION
%
%
%

targetVarL = lower(targetVar);

% Check if output directories exist, if not create them
for i = 1:numel(targetVarL)
    % Preallocate variables for efficiency
    learningDatesDate = table2array(learningDates(:,'date'));
    learningData      = table2array(learningDates(:,i+1));

    if targetDim ~= 1
        imgLength = size(learningData{1},1);
        imgWidth  = size(learningData{1},2);
        selectedImages = NaN(imgLength, imgWidth, size(sortedDatesOptim{1,2}, 1));
        %resultImages   = cell(size(sortedDates, 1), 1);
        imagesSynAll = single(NaN(imgLength,imgWidth,size(sortedDatesOptim,1)));
    else
        selectedImages = NaN(size(sortedDatesOptim{1,2}, 1),1);
        imagesSynAll   = NaN(size(sortedDatesOptim,1),1);
    end
    map = imagesSynAll;

    for rowIndex = 1:size(sortedDatesOptim,1)
        % Find the index of the current image in the Dates variable
        [~, dateIndex] = ismember(sortedDatesOptim{rowIndex,2},learningDatesDate);
        % Select the K best image from the Learning dataset and add it to selectedImages
        for imageIndex = 1:nbImages %length(sortedDates{rowIndex,2})
            if nbImages ~= length(sortedDatesOptim{rowIndex,2}) && imageIndex == 1
                warning(['nbImages .ne. number of available analogues (' num2str(nbImages) ' vs ' num2str(length(sortedDatesOptim{rowIndex,2})) ')'])
            end
            if targetDim ~= 1
                selectedImages(:,:,imageIndex) = learningData{dateIndex(imageIndex)};
            else
                selectedImages(imageIndex) = learningData(dateIndex(imageIndex));
            end
        end
        if generationType == 1
            % Calculate the mode and save it to resultImages
            if targetDim ~= 1
                resultImages = mode(selectedImages,3);
            else
                resultImages = mode(selectedImages);
            end
        elseif generationType == 2
            % Calculate the mean and save it to resultImages
            sortedDatesOptim{rowIndex,3}(sortedDatesOptim{rowIndex,3} == 0) = eps;
            selectedDist = 1./sortedDatesOptim{rowIndex,3}(1:nbImages);
            % Normalize the selectedDist values
            normalizedWeights = selectedDist / sum(selectedDist);
            % Perform element-wise multiplication with the weights
            if targetDim ~= 1
                weightedImages = bsxfun(@times, selectedImages, reshape(normalizedWeights, 1, 1, nbImages)); %length(sortedDates{rowIndex,2})
                resultImages = sum(weightedImages,3);
            else
                weightedImages = selectedImages .* normalizedWeights;
                resultImages = sum(weightedImages);
            end
        elseif generationType == 3
            % Calculate the mean and save it to resultImages
            if targetDim ~= 1
                resultImages = mean(selectedImages,3);
            else
                resultImages = mean(selectedImages);
            end
        elseif generationType == 4
            % Calculate the median and save it to resultImages
            if targetDim ~= 1
                resultImages = median(selectedImages,3);
            else
                resultImages = median(selectedImages);
            end
        else
            error('Generation type not defined or not appropriate for parameter optimisation!')
        end
        if targetDim ~= 1
            map(:,:,rowIndex) = resultImages;
        else
            map(rowIndex) = resultImages;
        end
    end
    if i == 1
        synImages.date = cell2mat(sortedDatesOptim(:,1));
    end
    synImages.(targetVarL(i)) = map;
end

end
