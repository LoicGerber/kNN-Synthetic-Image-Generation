function rawData = convertGeoTIFFtoStructure(targetVar,climateVars,addVars,rawDir,inputDir)

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
        tiffFile = fullfile(inputDir, files(i).name);
        [~, varname, ~] = fileparts(tiffFile);
        if ismember(lower(string(varname)), lower(varsAll))
            varnameID = strcat(varname, 'Index');
            disp(strcat("  Processing '", varname, "' data..."))

            % Read GeoTIFF information
            info = geotiffinfo(tiffFile);

            % Read the variable data
            data = geotiffread(tiffFile);

            % Get the dates (assuming file names are in yyyymmdd format)
            dates = str2double(regexp(files(i).name, '\d{8}', 'match'));

            % Create the output structure
            rawData.(varname)   = data;
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
