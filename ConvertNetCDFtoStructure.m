function rawData = ConvertNetCDFtoStructure(var,vars,addVars,rawDir,inputDir)

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

tic

% Check if output directories exist, if not create them
if ~exist(inputDir, 'dir')
    mkdir(inputDir)
end

varsAll = [var vars addVars];

% Get a list of all NetCDF files in the input directory
files = dir(fullfile(rawDir, '*.nc'));

rawData  = struct();
rawData  = structfun(@(x) [], rawData,  'UniformOutput', false);

for i = 1:length(files)

    % Open the NetCDF file
    ncfile = fullfile(rawDir, files(i).name);
    [~, varname, ~] = fileparts(ncfile);
    if ismember(lower(string(varname)),lower(varsAll))
        varnameID = strcat(varname,'Index');
        disp(strcat("  Processing '", varname,"' data..."))
        ncid = netcdf.open(ncfile,'NC_NOWRITE'); % Open NetCDF in read-only

        % Get the variable IDs and dimensions
        varid = netcdf.inqVarID(ncid,varname); % Return ID associated with variable name
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
            data{j}   = netcdf.getVar(ncid,varid,[0 0 time_idx-1],[dims{1,2} dims{2,2} 1],'double');
            dates{j}  = datestr(day_start,'yyyymmdd');
        end

        data  = cellfun(@transpose,data,'UniformOutput',false);
        dates = str2double(dates);

        % Create the output structure
        rawData.(varname)   = data;
        rawData.(varnameID) = dates;
    else
        disp(strcat("  '",varname,"' not in variables of interest, file not processed..."))
    end
end

% disp('Saving allVariables.mat file...')
% allVarsSave = fullfile(inputDir, 'allVariables.mat');
% save(allVarsSave, '-struct', 'output_data');
% disp('All variables processed and saved to allVariables.mat...')

toc

end
