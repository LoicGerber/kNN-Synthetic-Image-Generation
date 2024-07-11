function rawData = convertRawDataToStructure(targetVar,targetDim,climateVars,addVars,rawDir,inputDir)

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
            varnameID = strcat(lower(varName),'Index');
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
                time_units   = strrep(time_units,'hours since ','');
                try
                    time_origin = datetime(time_units,'InputFormat','yyyy-MM-dd hh:mm:ss');
                catch
                    time_origin = datetime(time_units,'InputFormat','yyyy-MM-dd hh:mm:ss.s');
                end
                time_datenum = datenum(time_origin + hours(time_data));
            else
                time_units   = strrep(time_units,'days since ','');
                time_origin = datetime(time_units,'InputFormat','yyyy-MM-dd');
                time_datenum = datenum(time_origin + days(time_data));
            end

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
            rawData.(lower(varName))   = data;
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
    varsTiff = varsAll(~ismember(lower(varsAll),lower(fields)));
    for k = 1:numel(varsTiff)
        varName = varsTiff(k);
        files = dir(fullfile(rawDir, varName, '*.tif'));
        if ~isempty(files) && targetDim ~= 1
            data  = cell(length(files),1);
            dates = nan(length(files),1);
            varPrevious = [];
            j = 0;
            for i = 1:length(files)
                j = j + 1;
                % Open the GeoTIFF file
                tiffFile = fullfile(rawDir, varName, files(i).name);
                %[~, filename, ~] = fileparts(tiffFile);
                [~, dateStr, ~] = fileparts(tiffFile);

                % Split the filename into variable and date
                %parts = split(filename, '_');
                %varName = parts{1};
                %dateStr = parts{2};

                if ismember(lower(string(varName)), lower(varsAll))
                    if ~strcmp(varName,varPrevious)
                        if ~isempty(varPrevious)
                            % Create the output structure
                            varnameID = strcat(varPrevious, 'Index');
                            data  = data(~cellfun('isempty',data));
                            dates = dates(~isnan(dates));
                            rawData.(varPrevious) = data;
                            rawData.(varnameID)   = dates;
                        end
                        disp(strcat("  Processing '", varName, "' data..."))
                        j = 1;
                    end
                    % Read the variable data
                    data{j,1} = single(readgeoraster(tiffFile));
                    % Convert date string to numeric format
                    dates(j,1) = str2double(dateStr);
                elseif ~strcmp(varName,varPrevious)
                    disp(strcat("  '", varName, "' not in variables of interest, file not processed..."))
                end
                varPrevious = varName;
            end
            % Create the output structure
            varnameID = strcat(lower(varName), 'Index');
            data  = data(~cellfun('isempty',data));
            dates = dates(~isnan(dates));
            rawData.(lower(varName))   = data;
            rawData.(varnameID) = dates;
        elseif targetDim == 1
            % Looking for 1D txt file
            txtFile = fullfile(rawDir, strcat(targetVar, '.txt'));
            if exist(txtFile, 'file')
                disp(strcat("  Processing '", targetVar, "' 1D data from text file..."))

                % Read the text file
                fid = fopen(txtFile, 'r');
                % Skip the first 5 header lines
                for i = 1:5
                    fgetl(fid);
                end

                % Initialize arrays for dates and data
                dates = [];
                data = [];

                % Read the rest of the file
                while ~feof(fid)
                    line = fgetl(fid);
                    if ischar(line)
                        parts = strsplit(line);
                        value = str2double(parts{6});
                        if value ~= -9999
                            dateVec = parts(1:5);
                            dates = [dates; strcat(dateVec(1), dateVec(2), dateVec(3))];
                            data = [data; value];
                        end
                    end
                end
                fclose(fid);

                % Create the output structure
                rawData.(lower(targetVar)) = data;
                rawData.(strcat(lower(targetVar), 'Index')) = str2double(dates);
            else
                error(strcat("No text file found for variable '", targetVar, "' in rawDir..."))
            end
        else
            error('No GeoTIFF or txt files in rawDir... Missing variable files or unknown format.')
        end
    end
end

end
