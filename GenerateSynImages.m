function GenerateSynImages(var,learningDates,sortedDates,R,outputDir,GenerationType,OutputType)

%
%
%
% REDO DOCUMENTATION
%
%
%

tic

outputDirImages = [outputDir 'syntheticImages\'];
% Check if output directories exist, if not create them
if ~exist(outputDirImages,'dir')
    mkdir(outputDirImages)
end

% Preallocate variables for efficiency
learningDatesDate = table2array(learningDates(:,'date'));
learningData      = table2array(learningDates(:,2));

imgLength = size(learningData{1},1);
imgWidth  = size(learningData{1},2);

selectedImages = NaN(imgLength, imgWidth, size(sortedDates{1,2}, 1));
resultImages   = cell(size(sortedDates, 1), 1);

% Loop through each row in sortedDates
for rowIndex = 1:size(sortedDates,1)
    % Loop through the second to the nbKNNth column of the current row
    % Find the index of the current image in the Dates variable
    [~, dateIndex] = ismember(sortedDates{rowIndex,2},learningDatesDate);
    % Select the Landsat image from the Landsat variable and add it to selectedImages
    for imageIndex = 1:length(sortedDates{rowIndex,2})
        selectedImages(:,:,imageIndex) = learningData{dateIndex(imageIndex)};
    end
    % Calculate either the mode or the mean of the selected Landsat images
    if GenerationType == 1
        % Calculate the mode and save it to resultImages
        resultImages{rowIndex} = mode(selectedImages(:,:,:),3);
    elseif GenerationType == 2
        % Calculate the mean and save it to resultImages
        resultImages{rowIndex} = mean(selectedImages(:,:,:),3);
    elseif GenerationType == 3
        % Calculate the median and save it to resultImages
        resultImages{rowIndex} = median(selectedImages(:,:,:),3);
    else
        error('Generation type not defined!')
    end
    if OutputType == 1
        % Write the resulting image to a GeoTIFF file
        outputBaseName = string(sortedDates(rowIndex,1)) + '.tif';
        fullDestinationFileName = fullfile(outputDirImages, outputBaseName);
        disp(['  Downlading image ' num2str(rowIndex) '/' num2str(size(sortedDates,1))])
        if isempty(R)
            disp('    Georeferencing files missing! Unreferenced output...')
            t = Tiff(fullDestinationFileName, 'w');
            tagstruct.ImageLength  = imgLength;
            tagstruct.ImageWidth   = imgWidth;
            tagstruct.Compression  = Tiff.Compression.None;
            tagstruct.SampleFormat = Tiff.SampleFormat.IEEEFP;
            tagstruct.Photometric  = Tiff.Photometric.MinIsBlack;
            tagstruct.BitsPerSample = 32;
            tagstruct.SamplesPerPixel = 1;
            tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
            t.setTag(tagstruct);
            t.write(single(resultImages{rowIndex,1}));
            t.close();
        else
            geotiffwrite(fullDestinationFileName,single(resultImages{rowIndex,1}),R,'TiffTags',struct('Compression',Tiff.Compression.None));
        end
    elseif OutputType == 2
        % Store the resulting image in a geolocated netCDF file
        outputBaseName = string(sortedDates(rowIndex,1)) + '.nc';
        fullDestinationFileName = fullfile(outputDirImages, outputBaseName);
        disp(['  Writing netCDF file ' num2str(rowIndex) '/' num2str(size(sortedDates,1))]);
        % Create a new netCDF file and define dimensions
        ncid = netcdf.create(fullDestinationFileName,'NETCDF4');
        dimid_x = netcdf.defDim(ncid,'x',size(resultImages{rowIndex,1},2));
        dimid_y = netcdf.defDim(ncid,'y',size(resultImages{rowIndex,1},1));
        dimid_lat = netcdf.defDim(ncid,'lat',size(resultImages{rowIndex,1},1));
        dimid_lon = netcdf.defDim(ncid,'lon',size(resultImages{rowIndex,1},2));
        % Define variables
        varid = netcdf.defVar(ncid,'resultImage','double',[dimid_y,dimid_x]);
        latid = netcdf.defVar(ncid,'lat','double',[dimid_lat,dimid_lon]);
        lonid = netcdf.defVar(ncid,'lon','double',[dimid_lat,dimid_lon]);
        % Define attributes
        netcdf.putAtt(ncid,varid,'long_name',var);
        netcdf.putAtt(ncid,latid,'long_name','latitude');
        netcdf.putAtt(ncid,latid,'units','degrees_north');
        netcdf.putAtt(ncid,lonid,'long_name','longitude');
        netcdf.putAtt(ncid,lonid,'units','degrees_east');
        % End definition mode
        netcdf.endDef(ncid);
        % Define metadata attributes
        if ~isempty(R)
            if ~isempty(R.GeographicCRS)
                ncwriteatt(fullDestinationFileName, '/', 'grid_mapping', 'geographic');
                ncwriteatt(fullDestinationFileName, '/', 'crs', R.GeographicCRS.Name);
            else
                ncwriteatt(fullDestinationFileName, '/', 'grid_mapping', 'geographic');
                ncwriteatt(fullDestinationFileName, '/', 'geographic', 'crs', 'WGS 84');
            end
            ncwriteatt(fullDestinationFileName, '/', 'xllcorner', R.LongitudeLimits(1));
            ncwriteatt(fullDestinationFileName, '/', 'yllcorner', R.LatitudeLimits(1));
            ncwriteatt(fullDestinationFileName, '/', 'cellsize', [R.CellExtentInLatitude R.CellExtentInLongitude]);
            ncwriteatt(fullDestinationFileName, '/', 'columnStart', R.ColumnsStartFrom);
            ncwriteatt(fullDestinationFileName, '/', 'rowStart', R.RowsStartFrom);
        end
        ncwriteatt(fullDestinationFileName, '/', 'date', num2str(sortedDates(rowIndex,1)),'Datatype','string');
        ncwriteatt(fullDestinationFileName, '/', 'nodata_value', -9999);
        % Write data to variable
        ncwrite(fullDestinationFileName, 'resultImage', single(resultImages{rowIndex,1}));
        % Close the netCDF file
        netcdf.close(ncid);
    else
        error('Unknown output type! Choose 1 for GeoTiff or 2 for NetCDF...')
    end
end

toc

end
