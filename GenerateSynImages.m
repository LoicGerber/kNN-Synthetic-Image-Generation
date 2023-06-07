function GenerateSynImages(var,learningDates,sortedDates,GeoRef,outputDir,GenerationType,optimisation,bootstrap,ensemble,OutputType)

%
%
%
% REDO DOCUMENTATION
%
%
%

outputDirImages = [outputDir 'syntheticImages\'];
% Check if output directories exist, if not create them
if ~exist(outputDirImages,'dir')
    mkdir(outputDirImages)
end
delete(fullfile(outputDirImages,'*'));

% Preallocate variables for efficiency
learningDatesDate = table2array(learningDates(:,'date'));
learningData      = table2array(learningDates(:,2));

imgLength = size(learningData{1},1);
imgWidth  = size(learningData{1},2);

selectedImages = NaN(imgLength, imgWidth, size(sortedDates{1,2}, 1));
%resultImages   = cell(size(sortedDates, 1), 1);

if bootstrap == true
    disp(['Bootstrap switch ON, using ' num2str(ensemble) ' ensembles'])
end
% Display progress
if optimisation == false
    progress = 0;
    if OutputType == 1
        fprintf(1,'Downloading synthetic GeoTiff images: %3.0f%%\n',progress);
    else
        fprintf(1,'Downloading synthetic images as NetCDF files: %3.0f%%\n',progress);
    end
end
% Loop through each row in sortedDates
for rowIndex = 1:size(sortedDates,1)
    if bootstrap == true
        % Find the index of the current image in the Dates variable
        [~, dateIndex] = ismember(sortedDates{rowIndex,2},learningDatesDate);
        % Select the Landsat image from the Landsat variable and add it to selectedImages
        for imageIndex = 1:length(sortedDates{rowIndex,2})
            selectedImages(:,:,imageIndex) = learningData{dateIndex(imageIndex)};
        end
        % Calculate either the mode or the mean of the selected Landsat images
        if GenerationType == 1
            % Calculate the mode and save it to resultImages
            resultImages = mode(selectedImages(:,:,:),3);
        elseif GenerationType == 2
            % Calculate the mean and save it to resultImages
            resultImages = mean(selectedImages(:,:,:),3);
        elseif GenerationType == 3
            % Calculate the median and save it to resultImages
            resultImages = median(selectedImages(:,:,:),3);
        else
            error('Generation type not defined!')
        end
        % Write the resulting image to a GeoTIFF file
        outputBaseName = string(sortedDates(rowIndex,1)) + '.tif';
        fullDestinationFileName = fullfile(outputDirImages, outputBaseName);
        %disp(['  Downlading image ' num2str(rowIndex) '/' num2str(size(sortedDates,1))])
        if isempty(GeoRef)
            %disp('    Georeferencing files missing! Unreferenced output...')
            t = Tiff(fullDestinationFileName, 'w');
            tagstruct.ImageLength         = imgLength;
            tagstruct.ImageWidth          = imgWidth;
            tagstruct.Compression         = Tiff.Compression.None;
            tagstruct.SampleFormat        = Tiff.SampleFormat.IEEEFP;
            tagstruct.Photometric         = Tiff.Photometric.MinIsBlack;
            tagstruct.BitsPerSample       = 32;
            tagstruct.SamplesPerPixel     = 1;
            tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
            t.setTag(tagstruct);
            t.write(single(resultImages));
            t.close();
        else
            geotiffwrite(fullDestinationFileName,single(resultImages),GeoRef,'TiffTags',struct('Compression',Tiff.Compression.None));
        end
        % bootstrap
        resultImages = NaN(imgLength, imgWidth, size(sortedDates{1,2}, 1));
        invDistance      = 1 ./ sortedDates{rowIndex,3};
        bootstrapWeights = normalize(invDistance,'range',[0.1 1]); % normalise distance (3) / std (4) to [0.1 1]
        %bootstrapWeights = invDistance/sum(invDistance);
        for bs = 1:ensemble
            bootstrapDates = randsample(sortedDates{rowIndex,2},numel(sortedDates{rowIndex,2}),true,bootstrapWeights);
            % Find the index of the current image in the Dates variable
            [~, dateIndex] = ismember(bootstrapDates,learningDatesDate);
            % Select the Landsat image from the Landsat variable and add it to selectedImages
            for imageIndex = 1:length(sortedDates{rowIndex,2})
                selectedImages(:,:,imageIndex) = learningData{dateIndex(imageIndex)};
            end
            % Calculate either the mode or the mean of the selected Landsat images
            if GenerationType == 1
                % Calculate the mode and save it to resultImages
                resultImages(:,:,bs) = mode(selectedImages(:,:,:),3);
            elseif GenerationType == 2
                % Calculate the mean and save it to resultImages
                resultImages(:,:,bs) = mean(selectedImages(:,:,:),3);
            elseif GenerationType == 3
                % Calculate the median and save it to resultImages
                resultImages(:,:,bs) = median(selectedImages(:,:,:),3);
            else
                error('Generation type not defined!')
            end
            % Write the resulting image to a GeoTIFF file
            outputBaseName = string(sortedDates(rowIndex,1)) + '_' + bs + '.tif';
            fullDestinationFileName = fullfile(outputDirImages, outputBaseName);
            %disp(['  Downlading image ' num2str(rowIndex) '/' num2str(size(sortedDates,1))])
            if isempty(GeoRef)
                %disp('    Georeferencing files missing! Unreferenced output...')
                t = Tiff(fullDestinationFileName, 'w');
                tagstruct.ImageLength         = imgLength;
                tagstruct.ImageWidth          = imgWidth;
                tagstruct.Compression         = Tiff.Compression.None;
                tagstruct.SampleFormat        = Tiff.SampleFormat.IEEEFP;
                tagstruct.Photometric         = Tiff.Photometric.MinIsBlack;
                tagstruct.BitsPerSample       = 32;
                tagstruct.SamplesPerPixel     = 1;
                tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
                t.setTag(tagstruct);
                t.write(single(resultImages(:,:,bs)));
                t.close();
            else
                geotiffwrite(fullDestinationFileName,single(resultImages(:,:,bs)),GeoRef,'TiffTags',struct('Compression',Tiff.Compression.None));
            end
        end
        resultImagesMean = mean(resultImages(:,:,:),3);
        % Write the resulting image to a GeoTIFF file
        outputBaseName = string(sortedDates(rowIndex,1)) + '_BootstrapAll.tif';
        fullDestinationFileName = fullfile(outputDirImages, outputBaseName);
        %disp(['  Downlading image ' num2str(rowIndex) '/' num2str(size(sortedDates,1))])
        if isempty(GeoRef)
            %disp('    Georeferencing files missing! Unreferenced output...')
            t = Tiff(fullDestinationFileName, 'w');
            tagstruct.ImageLength         = imgLength;
            tagstruct.ImageWidth          = imgWidth;
            tagstruct.Compression         = Tiff.Compression.None;
            tagstruct.SampleFormat        = Tiff.SampleFormat.IEEEFP;
            tagstruct.Photometric         = Tiff.Photometric.MinIsBlack;
            tagstruct.BitsPerSample       = 32;
            tagstruct.SamplesPerPixel     = 1;
            tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
            t.setTag(tagstruct);
            t.write(single(resultImagesMean));
            t.close();
        else
            geotiffwrite(fullDestinationFileName,single(resultImagesMean),GeoRef,'TiffTags',struct('Compression',Tiff.Compression.None));
        end
    else
        % Find the index of the current image in the Dates variable
        [~, dateIndex] = ismember(sortedDates{rowIndex,2},learningDatesDate);
        % Select the Landsat image from the Landsat variable and add it to selectedImages
        for imageIndex = 1:length(sortedDates{rowIndex,2})
            selectedImages(:,:,imageIndex) = learningData{dateIndex(imageIndex)};
        end
        % Calculate either the mode or the mean of the selected Landsat images
        if GenerationType == 1
            % Calculate the mode and save it to resultImages
            resultImages = mode(selectedImages(:,:,:),3);
        elseif GenerationType == 2
            % Calculate the mean and save it to resultImages
            resultImages = mean(selectedImages(:,:,:),3);
        elseif GenerationType == 3
            % Calculate the median and save it to resultImages
            resultImages = median(selectedImages(:,:,:),3);
        else
            error('Generation type not defined!')
        end
        if OutputType == 1
            % Write the resulting image to a GeoTIFF file
            outputBaseName = string(sortedDates(rowIndex,1)) + '.tif';
            fullDestinationFileName = fullfile(outputDirImages, outputBaseName);
            %disp(['  Downlading image ' num2str(rowIndex) '/' num2str(size(sortedDates,1))])
            if isempty(GeoRef)
                %disp('    Georeferencing files missing! Unreferenced output...')
                t = Tiff(fullDestinationFileName, 'w');
                tagstruct.ImageLength         = imgLength;
                tagstruct.ImageWidth          = imgWidth;
                tagstruct.Compression         = Tiff.Compression.None;
                tagstruct.SampleFormat        = Tiff.SampleFormat.IEEEFP;
                tagstruct.Photometric         = Tiff.Photometric.MinIsBlack;
                tagstruct.BitsPerSample       = 32;
                tagstruct.SamplesPerPixel     = 1;
                tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
                t.setTag(tagstruct);
                t.write(single(resultImages));
                t.close();
            else
                geotiffwrite(fullDestinationFileName,single(resultImages),GeoRef,'TiffTags',struct('Compression',Tiff.Compression.None));
            end
        elseif OutputType == 2
            % Assign the CRS value
            crs_wkt = wktstring(GeoRef.GeographicCRS);
            % Extract the EPSG code from the WKT string using regular expressions
            expression = 'ID\["EPSG",(\d+)\]';
            tokens = regexp(crs_wkt, expression, 'tokens');
            crs_value = tokens{1};
            % Store the resulting image in a geolocated netCDF file
            outputBaseName = string(sortedDates(rowIndex,1)) + '.nc';
            fullDestinationFileName = fullfile(outputDirImages, outputBaseName);
            %disp(['  Writing netCDF file ' num2str(rowIndex) '/' num2str(size(sortedDates,1))]);
            % Create a new netCDF file and define dimensions
            ncid       = netcdf.create(fullDestinationFileName,'NETCDF4');
            dimid_lat  = netcdf.defDim(ncid,'lat',GeoRef.RasterSize(1));
            dimid_lon  = netcdf.defDim(ncid,'lon',GeoRef.RasterSize(2));
            dimid_time = netcdf.defDim(ncid,'time',1);
            % Define variables
            %varid = netcdf.defVar(ncid,var,'double',[dimid_lat,dimid_lon]);
            varid  = netcdf.defVar(ncid,var,'double',[dimid_lon,dimid_lat]);
            timeid = netcdf.defVar(ncid,'time','double',dimid_time);
            latid  = netcdf.defVar(ncid,'lat','double',dimid_lat);
            lonid  = netcdf.defVar(ncid,'lon','double',dimid_lon);
            % Define attributes
            netcdf.putAtt(ncid,varid,'long_name',var);
            netcdf.putAtt(ncid,timeid,'long_name','time');
            netcdf.putAtt(ncid,timeid,'units','days since 1970-01-01 00:00:00');
            netcdf.putAtt(ncid,latid,'long_name','latitude');
            netcdf.putAtt(ncid,latid,'units','degrees_north');
            netcdf.putAtt(ncid,lonid,'long_name','longitude');
            netcdf.putAtt(ncid,lonid,'units','degrees_east');
            % Assign the CRS as a global attribute to the netCDF file
            netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'crs_wkt', crs_wkt);
            netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'crs', crs_value);
            % End definition mode
            netcdf.endDef(ncid);

            % Assign latitude and longitude values
            %lat = GeoRef.LatitudeLimits(2)  : -GeoRef.CellExtentInLatitude : GeoRef.LatitudeLimits(1);
            %lon = GeoRef.LongitudeLimits(1) : GeoRef.CellExtentInLongitude : GeoRef.LongitudeLimits(2);

            % Assign latitude and longitude values
            lat_start = GeoRef.LatitudeLimits(2);
            lat_end   = GeoRef.LatitudeLimits(1);
            lon_start = GeoRef.LongitudeLimits(1);
            lon_end   = GeoRef.LongitudeLimits(2);
            lat_size  = size(resultImages, 1);
            lon_size  = size(resultImages, 2);
            lat_step  = (lat_end - lat_start) / lat_size;
            lon_step  = (lon_end - lon_start) / lon_size;
            lat       = lat_start:lat_step:lat_end;
            lon       = lon_start:lon_step:lon_end;
            % Adjust the size of lat and lon vectors to match the image dimensions
            lat       = lat(1:lat_size);
            lon       = lon(1:lon_size);
            % Assign latitude and longitude values to the corresponding variables
            netcdf.putVar(ncid,latid,lat);
            netcdf.putVar(ncid,lonid,lon);
            % Assign date
            dateStr  = convertStringsToChars(string(sortedDates{rowIndex, 1}));
            yearStr  = dateStr(1:4);
            monthStr = dateStr(5:6);
            dayStr   = dateStr(7:8);
            dateStrFormatted = [yearStr '-' monthStr '-' dayStr];
            time = datenum(dateStrFormatted, 'yyyy-mm-dd');
            netcdf.putVar(ncid, timeid, time);
            % Define metadata attributes
            ncwriteatt(fullDestinationFileName, '/', 'crs', crs_value);
            ncwriteatt(fullDestinationFileName, '/', 'xllcorner', GeoRef.LongitudeLimits(1));
            ncwriteatt(fullDestinationFileName, '/', 'yllcorner', GeoRef.LatitudeLimits(2));
            ncwriteatt(fullDestinationFileName, '/', 'origin', [GeoRef.LongitudeLimits(1) GeoRef.LatitudeLimits(2)]);
            ncwriteatt(fullDestinationFileName, '/', 'cellsize', [GeoRef.CellExtentInLatitude GeoRef.CellExtentInLongitude]);
            ncwriteatt(fullDestinationFileName, '/', 'columnStart', GeoRef.ColumnsStartFrom);
            ncwriteatt(fullDestinationFileName, '/', 'rowStart', GeoRef.RowsStartFrom);
            ncwriteatt(fullDestinationFileName, '/', 'date', string(sortedDates{rowIndex, 1}),'Datatype','string');
            ncwriteatt(fullDestinationFileName, '/', 'nodata_value', -9999);
            % Write data to variable
            ncwrite(fullDestinationFileName, var, single(resultImages)');
            % Close the netCDF file
            netcdf.close(ncid);
        else
            error('Unknown output type! Choose 1 for GeoTiff or 2 for NetCDF...')
        end
    end
    if optimisation == false
        % Display computation progress
        progress = (100*(rowIndex/size(sortedDates,1)));
        fprintf(1,'\b\b\b\b%3.0f%%',progress);
    end
end

if optimisation == false
    fprintf('\n')
end

end
