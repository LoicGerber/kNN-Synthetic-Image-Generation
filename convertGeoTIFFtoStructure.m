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
    unzip(dir(fullfile(rawDir, 'data.zip')));

    % Get a list of all GeoTIFF files in the input directory
    files = dir(fullfile(rawDir, 'data', '*.tif'));

    rawData  = struct();
    rawData  = structfun(@(x) [], rawData,  'UniformOutput', false);

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
            [data, R] = readgeoraster(tiffFile);

            % Convert date string to numeric format
            dates = str2double(dateStr);

            % Create the output structure
            rawData.(varname)   = struct('data', data, 'R', R);
            rawData.(varnameID) = dates;
        else
            disp(strcat("  '", varname, "' not in variables of interest, file not processed..."))
        end
    end

    % disp('Saving allVariables.mat file...')
    % allVarsSave = fullfile(inputDir, 'allVariables.mat');
    % save(allVarsSave, '-struct', 'output_data');
    % disp('All variables processed and saved to allVariables.mat...')

end
