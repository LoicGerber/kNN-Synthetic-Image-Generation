function rawData = convertRawDataToStructure(targetVar,climateVars,addVars,rawDir,inputDir)

%
%
%
% REDO DOCUMENTATION
%
%
%

% Convert NetCDF data to structure
% Input directory containing NetCDF files
% files must be named with the variable. Example: precipitation data must
% be called pre.nc

% Check if output directories exist, if not create them
if ~exist(inputDir, 'dir')
    mkdir(inputDir)
end

varsAll = [targetVar climateVars addVars];
fields  = [];

rawData  = struct();
rawData  = structfun(@(x) [], rawData,  'UniformOutput', false);

% Get a list of all NetCDF files in the input directory
files = dir(fullfile(rawDir, '*.nc'));
if ~isempty(files)
    for i = 1:length(files)
        % Open the NetCDF file
        ncfile = fullfile(rawDir, files(i).name);
        [~, varName, ~] = fileparts(ncfile);
        if ismember(lower(string(varName)),lower(varsAll))
            varnameID = strcat(varName,'Index');
            disp(strcat("  Processing '", varName,"' data..."))
            ncid = netcdf.open(ncfile,'NC_NOWRITE'); % Open NetCDF in read-only

            % Get the variable IDs and dimensions
            varid = netcdf.inqVarID(ncid,varName); % Return ID associated with variable name
            [~,~,dimids,~] = netcdf.inqVar(ncid,varid); % returns 4 outputs, 1-2-4 are ignored (1-var name,2-data type,3-dimension IDs,4-num of var attributes)
            dims = cell(length(dimids),1);
            for j=1:length(dimids)
                [dims{j,1}, dims{j,2}] = netcdf.inqDim(ncid,dimids(j));
            end

            % Get the time variable and convert to MATLAB datenum
            time_varid   = netcdf.inqVarID(ncid,'time');
            time_data    = netcdf.getVar(ncid,time_varid);
            time_units   = netcdf.getAtt(ncid,time_varid,'units');
            if strlength(time_units) > 21
                time_origin = datetime(strrep(time_units,'days since ',''),'InputFormat','yyyy-MM-dd hh:mm:ss');
            else
                time_origin = datetime(strrep(time_units,'days since ',''),'InputFormat','yyyy-MM-dd');
            end
            time_datenum = datenum(time_origin + days(time_data));

            % Get the variable data for each day and store in a cell array
            num_days = length(time_datenum);
            data     = cell(num_days,1);
            dates    = cell(num_days,1);
            for j=1:num_days
                day_start = time_datenum(j);
                day_end   = day_start + 1;
                time_idx  = find(time_datenum >= day_start & time_datenum < day_end);
                data{j}   = netcdf.getVar(ncid,varid,[0 0 time_idx-1],[dims{1,2} dims{2,2} 1],'single');
                dates{j}  = datestr(day_start,'yyyymmdd');
            end

            data  = cellfun(@transpose,data,'UniformOutput',false);
            dates = str2double(dates);

            % Create the output structure
            rawData.(varName)   = data;
            rawData.(varnameID) = dates;
        else
            disp(strcat("  '",varName,"' not in variables of interest, file not processed..."))
        end
    end
    fields = string(fieldnames(rawData))';
    fields = fields(ismember(fields,lower(varsAll)));
end
if isempty(files) || numel(fields) < numel(varsAll)
    % Get a list of all GeoTIFF files in the input directory
    files = dir(fullfile(rawDir, 'data', '*.tif'));
    if isempty(files)
        error('No GeoTIFF files in rawDir... Missing variable files or unknown formt.')
    end

    data  = cell(length(files),1);
    dates = nan(length(files),1);
    
    varPrevious = [];

    for i = 1:length(files)
        % Open the GeoTIFF file
        tiffFile = fullfile(rawDir, 'data', files(i).name);
        [~, filename, ~] = fileparts(tiffFile);

        % Split the filename into variable and date
        parts = split(filename, '_');
        varName = parts{1};
        dateStr = parts{2};

        if ismember(lower(string(varName)), lower(varsAll))
            varnameID = strcat(varName, 'Index');
            if ~strcmp(varName,varPrevious)
                disp(strcat("  Processing '", varName, "' data..."))
            end
            % Read the variable data
            data{i,1} = single(readgeoraster(tiffFile));

            % Convert date string to numeric format
            dates(i,1) = str2double(dateStr);
        else
            disp(strcat("  '", varName, "' not in variables of interest, file not processed..."))
        end
        varPrevious = varName;
    end
    % Create the output structure
    rawData.(varName)   = data;
    rawData.(varnameID) = dates;
end

end
