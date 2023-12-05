function rawData = convertGeoTIFFtoStructure(targetVar,climateVars,addVars,rawDir,inputDir)
    
    %
    %
    %
    % REDO DOCUMENTATION
    %
    %
    %
    
    % Convert GeoTIFF data to structure
    % Input directory containing GeoTIFF files
    % files must be named with the variable. Example: precipitation data must
    % be called pre_%DATE%.tiff

    % Check if output directories exist, if not create them
    if ~exist(inputDir, 'dir')
        mkdir(inputDir)
    end

    varsAll = [targetVar climateVars addVars];

    % Extract contents of the zip folder
    unzip(fullfile(rawDir, 'data.zip'),rawDir);

    % Get a list of all GeoTIFF files in the input directory
    files = dir(fullfile(rawDir, 'data', '*.tif'));

    rawData  = struct();
    rawData  = structfun(@(x) [], rawData,  'UniformOutput', false);
    
    data  = cell(length(files),1);
    dates = nan(length(files),1);

    for i = 1:length(files)
        % Open the GeoTIFF file
        tiffFile = fullfile(rawDir, 'data', files(i).name);
        [~, filename, ~] = fileparts(tiffFile);
        
        % Split the filename into variable and date
        parts = split(filename, '_');
        varname = parts{1};
        dateStr = parts{2};

        if ismember(lower(string(varname)), lower(varsAll))
            varnameID = strcat(varname, 'Index');
            disp(strcat("  Processing '", varname, "' data..."))

            % Read the variable data
            data{i,1} = single(readgeoraster(tiffFile));

            % Convert date string to numeric format
            dates(i,1) = str2double(dateStr);
        else
           disp(strcat("  '", varname, "' not in variables of interest, file not processed..."))
        end
    end

    % Create the output structure
    rawData.(varname)   = data;
    rawData.(varnameID) = dates;

    % disp('Saving allVariables.mat file...')
    % allVarsSave = fullfile(inputDir, 'allVariables.mat');
    % save(allVarsSave, '-struct', 'output_data');
    % disp('All variables processed and saved to allVariables.mat...')

end
