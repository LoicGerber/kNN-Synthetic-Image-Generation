function synImages = generateSynImages(targetVar,targetDim,learningDates,sortedDates,mps,geoRef,outputDir,generationType,nanValue,validation,optimisation,bootstrap,bsSaveAll,nbImages,ensemble,outputType)

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
    disp(strcat("Processing variable '",convertStringsToChars(targetVar(i)),"'..."))

    % Preallocate variables for efficiency
    learningDatesDate = table2array(learningDates(:,'date'));
    learningData      = table2array(learningDates(:,i+1));

    if generationType == 5
        % Mask with shape of image
        mask = double(~isnan(learningData{1,1}));
        mask(mask==1) = nan;
        mask(mask==0) = nanValue;
        % Map with numbered pixels, to prevent pixel movement
        values = 1:(size(mask,1) * size(mask,2));
        loc = reshape(values, size(mask,1), size(mask,2));
        mask = cat(3, mask, loc);
        % Kernel creation
        ki = ones(mps.kernel_dims);
        for w = 1:length(mps.kernel_weights)
            ki(:,:,w) = ki(:,:,w) .* mps.kernel_weights(w);
        end
    end
    
    if targetDim ~= 1
        imgLength = size(learningData{1},1);
        imgWidth  = size(learningData{1},2);

        GeoRef = geoRef.(targetVar(i));

        selectedImages = single(NaN(imgLength, imgWidth, size(sortedDates{1,2}, 1)));
        %resultImages   = cell(size(sortedDates, 1), 1);
        imagesSynAll = single(NaN(imgLength,imgWidth,size(sortedDates,1)));
        map          = imagesSynAll;
        varMap       = imagesSynAll;
        availablePix = imagesSynAll;
        varianceBS   = imagesSynAll;

        if bootstrap == true
            imagesSynAll = cell(size(sortedDates,1),1);
            disp(['  Bootstrap switch ON, using ' num2str(ensemble) ' ensembles'])
        end
        % Display progress
        if optimisation == false
            progress = 0;
            if outputType == 1
                fprintf(1,'  Downloading synthetic GeoTiff images: %3.0f%%\n',progress);
            else
                fprintf(1,'  Downloading synthetic images as NetCDF file: %3.0f%%\n',progress);
            end
        end

        % netCDF file definition
        if outputType == 2 && bootstrap == false
            % Define the main netCDF file
            outputBaseName = strcat(targetVarL(i),'.nc');
            fullDestinationFileName = fullfile(outputDir, outputBaseName);
            % Assign the CRS value
            %try
            crs_wkt = wktstring(GeoRef.GeographicCRS);
            % Extract the EPSG code from the WKT string using regular expressions
            expression = 'ID\["EPSG",(\d+)\]';
            tokens = regexp(crs_wkt, expression, 'tokens');
            crs_value = tokens{1};
            %catch

            %end
            % Create the main netCDF file and define dimensions
            ncid = netcdf.create(fullDestinationFileName, 'NETCDF4');
            dimid_lat = netcdf.defDim(ncid, 'lat', GeoRef.RasterSize(1));
            dimid_lon = netcdf.defDim(ncid, 'lon', GeoRef.RasterSize(2));
            dimid_time = netcdf.defDim(ncid, 'time', netcdf.getConstant('NC_UNLIMITED'));
            % Define variables
            varid = netcdf.defVar(ncid, targetVarL(i), 'double', [dimid_lon, dimid_lat, dimid_time]);
            timeid = netcdf.defVar(ncid, 'time', 'double', dimid_time);
            latid = netcdf.defVar(ncid, 'lat', 'double', dimid_lat);
            lonid = netcdf.defVar(ncid, 'lon', 'double', dimid_lon);
            % Define attributes
            netcdf.putAtt(ncid, varid, 'long_name', targetVarL(i));
            netcdf.putAtt(ncid, varid, '_FillValue', -999);
            netcdf.putAtt(ncid, timeid, 'long_name', 'time');
            netcdf.putAtt(ncid, timeid, 'units', 'days since 1970-01-01');
            netcdf.putAtt(ncid, timeid, 'calendar', 'proleptic_gregorian');
            netcdf.putAtt(ncid, latid, 'long_name', 'latitude');
            netcdf.putAtt(ncid, latid, 'units', 'degrees_north');
            netcdf.putAtt(ncid, lonid, 'long_name', 'longitude');
            netcdf.putAtt(ncid, lonid, 'units', 'degrees_east');
            % Assign the CRS as a global attribute to the netCDF file
            netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'crs_wkt', crs_wkt);
            netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'crs', crs_value);
            % End definition mode
            netcdf.endDef(ncid);
            % Assign latitude and longitude values
            lat_start = GeoRef.LatitudeLimits(2);
            lat_end   = GeoRef.LatitudeLimits(1);
            lon_start = GeoRef.LongitudeLimits(1);
            lon_end   = GeoRef.LongitudeLimits(2);
            lat_size  = imgLength;
            lon_size  = imgWidth;
            lat_step  = (lat_end - lat_start) / lat_size;
            lon_step  = (lon_end - lon_start) / lon_size;
            lat       = lat_start:lat_step:lat_end;
            lon       = lon_start:lon_step:lon_end;
            % Adjust the size of lat and lon vectors to match the image dimensions
            lat       = lat(1:lat_size)+(lat_step/2);
            lon       = lon(1:lon_size)+(lon_step/2);
            % Assign latitude and longitude values to the corresponding variables
            netcdf.putVar(ncid,latid,lat);
            netcdf.putVar(ncid,lonid,lon);
        elseif bootstrap == true
            % ---- MININMAL ----
            % Define the main netCDF file
            outputBaseNameMin = strcat(targetVarL(i),'_min.nc');
            fullDestinationFileNameMin = fullfile(outputDir, outputBaseNameMin);
            % Assign the CRS value
            crs_wkt = wktstring(GeoRef.GeographicCRS);
            % Extract the EPSG code from the WKT string using regular expressions
            expression = 'ID\["EPSG",(\d+)\]';
            tokens = regexp(crs_wkt, expression, 'tokens');
            crs_value = tokens{1};
            % Create the main netCDF file and define dimensions
            ncid_min = netcdf.create(fullDestinationFileNameMin, 'NETCDF4');
            dimid_lat = netcdf.defDim(ncid_min, 'lat', GeoRef.RasterSize(1));
            dimid_lon = netcdf.defDim(ncid_min, 'lon', GeoRef.RasterSize(2));
            dimid_time = netcdf.defDim(ncid_min, 'time', netcdf.getConstant('NC_UNLIMITED'));
            % Define variables
            varid = netcdf.defVar(ncid_min, targetVarL(i), 'double', [dimid_lon, dimid_lat, dimid_time]);
            timeid = netcdf.defVar(ncid_min, 'time', 'double', dimid_time);
            latid = netcdf.defVar(ncid_min, 'lat', 'double', dimid_lat);
            lonid = netcdf.defVar(ncid_min, 'lon', 'double', dimid_lon);
            % Define attributes
            netcdf.putAtt(ncid_min, varid, 'long_name', targetVarL(i));
            netcdf.putAtt(ncid_min, varid, '_FillValue', -999);
            netcdf.putAtt(ncid_min, timeid, 'long_name', 'time');
            netcdf.putAtt(ncid_min, timeid, 'units', 'days since 1970-01-01');
            netcdf.putAtt(ncid_min, timeid, 'calendar', 'proleptic_gregorian');
            netcdf.putAtt(ncid_min, latid, 'long_name', 'latitude');
            netcdf.putAtt(ncid_min, latid, 'units', 'degrees_north');
            netcdf.putAtt(ncid_min, lonid, 'long_name', 'longitude');
            netcdf.putAtt(ncid_min, lonid, 'units', 'degrees_east');
            % Assign the CRS as a global attribute to the netCDF file
            netcdf.putAtt(ncid_min, netcdf.getConstant('NC_GLOBAL'), 'crs_wkt', crs_wkt);
            netcdf.putAtt(ncid_min, netcdf.getConstant('NC_GLOBAL'), 'crs', crs_value);
            % End definition mode
            netcdf.endDef(ncid_min);
            % Assign latitude and longitude values
            lat_start = GeoRef.LatitudeLimits(2);
            lat_end   = GeoRef.LatitudeLimits(1);
            lon_start = GeoRef.LongitudeLimits(1);
            lon_end   = GeoRef.LongitudeLimits(2);
            lat_size  = imgLength;
            lon_size  = imgWidth;
            lat_step  = (lat_end - lat_start) / lat_size;
            lon_step  = (lon_end - lon_start) / lon_size;
            lat       = lat_start:lat_step:lat_end;
            lon       = lon_start:lon_step:lon_end;
            % Adjust the size of lat and lon vectors to match the image dimensions
            lat       = lat(1:lat_size)+(lat_step/2);
            lon       = lon(1:lon_size)+(lon_step/2);
            % Assign latitude and longitude values to the corresponding variables
            netcdf.putVar(ncid_min,latid,lat);
            netcdf.putVar(ncid_min,lonid,lon);
            % ---- DETERMINISTIC ----
            % Define the main netCDF file
            outputBaseNameDet = strcat(targetVarL(i),'_det.nc');
            fullDestinationFileNameDet = fullfile(outputDir, outputBaseNameDet);
            % Assign the CRS value
            crs_wkt = wktstring(GeoRef.GeographicCRS);
            % Extract the EPSG code from the WKT string using regular expressions
            expression = 'ID\["EPSG",(\d+)\]';
            tokens = regexp(crs_wkt, expression, 'tokens');
            crs_value = tokens{1};
            % Create the main netCDF file and define dimensions
            ncid_det = netcdf.create(fullDestinationFileNameDet, 'NETCDF4');
            dimid_lat = netcdf.defDim(ncid_det, 'lat', GeoRef.RasterSize(1));
            dimid_lon = netcdf.defDim(ncid_det, 'lon', GeoRef.RasterSize(2));
            dimid_time = netcdf.defDim(ncid_det, 'time', netcdf.getConstant('NC_UNLIMITED'));
            % Define variables
            varid = netcdf.defVar(ncid_det, targetVarL(i), 'double', [dimid_lon, dimid_lat, dimid_time]);
            timeid = netcdf.defVar(ncid_det, 'time', 'double', dimid_time);
            latid = netcdf.defVar(ncid_det, 'lat', 'double', dimid_lat);
            lonid = netcdf.defVar(ncid_det, 'lon', 'double', dimid_lon);
            % Define attributes (similar to your existing code)
            netcdf.putAtt(ncid_det, varid, 'long_name', targetVarL(i));
            netcdf.putAtt(ncid_det, varid, '_FillValue', -999);
            netcdf.putAtt(ncid_det, timeid, 'long_name', 'time');
            netcdf.putAtt(ncid_det, timeid, 'units', 'days since 1970-01-01');
            netcdf.putAtt(ncid_det, timeid, 'calendar', 'proleptic_gregorian');
            netcdf.putAtt(ncid_det, latid, 'long_name', 'latitude');
            netcdf.putAtt(ncid_det, latid, 'units', 'degrees_north');
            netcdf.putAtt(ncid_det, lonid, 'long_name', 'longitude');
            netcdf.putAtt(ncid_det, lonid, 'units', 'degrees_east');
            % Assign the CRS as a global attribute to the netCDF file
            netcdf.putAtt(ncid_det, netcdf.getConstant('NC_GLOBAL'), 'crs_wkt', crs_wkt);
            netcdf.putAtt(ncid_det, netcdf.getConstant('NC_GLOBAL'), 'crs', crs_value);
            % End definition mode
            netcdf.endDef(ncid_det);
            % Assign latitude and longitude values
            lat_start = GeoRef.LatitudeLimits(2);
            lat_end   = GeoRef.LatitudeLimits(1);
            lon_start = GeoRef.LongitudeLimits(1);
            lon_end   = GeoRef.LongitudeLimits(2);
            lat_size  = imgLength;
            lon_size  = imgWidth;
            lat_step  = (lat_end - lat_start) / lat_size;
            lon_step  = (lon_end - lon_start) / lon_size;
            lat       = lat_start:lat_step:lat_end;
            lon       = lon_start:lon_step:lon_end;
            % Adjust the size of lat and lon vectors to match the image dimensions
            lat       = lat(1:lat_size)+(lat_step/2);
            lon       = lon(1:lon_size)+(lon_step/2);
            % Assign latitude and longitude values to the corresponding variables
            netcdf.putVar(ncid_det,latid,lat);
            netcdf.putVar(ncid_det,lonid,lon);
            % ---- MAXINMAL ----
            % Define the main netCDF file
            outputBaseNameMax = strcat(targetVarL(i),'_max.nc');
            fullDestinationFileNameMax = fullfile(outputDir, outputBaseNameMax);
            % Assign the CRS value
            crs_wkt = wktstring(GeoRef.GeographicCRS);
            % Extract the EPSG code from the WKT string using regular expressions
            expression = 'ID\["EPSG",(\d+)\]';
            tokens = regexp(crs_wkt, expression, 'tokens');
            crs_value = tokens{1};
            % Create the main netCDF file and define dimensions
            ncid_max = netcdf.create(fullDestinationFileNameMax, 'NETCDF4');
            dimid_lat = netcdf.defDim(ncid_max, 'lat', GeoRef.RasterSize(1));
            dimid_lon = netcdf.defDim(ncid_max, 'lon', GeoRef.RasterSize(2));
            dimid_time = netcdf.defDim(ncid_max, 'time', netcdf.getConstant('NC_UNLIMITED'));
            % Define variables
            varid = netcdf.defVar(ncid_max, targetVarL(i), 'double', [dimid_lon, dimid_lat, dimid_time]);
            timeid = netcdf.defVar(ncid_max, 'time', 'double', dimid_time);
            latid = netcdf.defVar(ncid_max, 'lat', 'double', dimid_lat);
            lonid = netcdf.defVar(ncid_max, 'lon', 'double', dimid_lon);
            % Define attributes (similar to your existing code)
            netcdf.putAtt(ncid_max, varid, 'long_name', targetVarL(i));
            netcdf.putAtt(ncid_max, varid, '_FillValue', -999);
            netcdf.putAtt(ncid_max, timeid, 'long_name', 'time');
            netcdf.putAtt(ncid_max, timeid, 'units', 'days since 1970-01-01');
            netcdf.putAtt(ncid_max, timeid, 'calendar', 'proleptic_gregorian');
            netcdf.putAtt(ncid_max, latid, 'long_name', 'latitude');
            netcdf.putAtt(ncid_max, latid, 'units', 'degrees_north');
            netcdf.putAtt(ncid_max, lonid, 'long_name', 'longitude');
            netcdf.putAtt(ncid_max, lonid, 'units', 'degrees_east');
            % Assign the CRS as a global attribute to the netCDF file
            netcdf.putAtt(ncid_max, netcdf.getConstant('NC_GLOBAL'), 'crs_wkt', crs_wkt);
            netcdf.putAtt(ncid_max, netcdf.getConstant('NC_GLOBAL'), 'crs', crs_value);
            % End definition mode
            netcdf.endDef(ncid_max);
            % Assign latitude and longitude values
            lat_start = GeoRef.LatitudeLimits(2);
            lat_end   = GeoRef.LatitudeLimits(1);
            lon_start = GeoRef.LongitudeLimits(1);
            lon_end   = GeoRef.LongitudeLimits(2);
            lat_size  = imgLength;
            lon_size  = imgWidth;
            lat_step  = (lat_end - lat_start) / lat_size;
            lon_step  = (lon_end - lon_start) / lon_size;
            lat       = lat_start:lat_step:lat_end;
            lon       = lon_start:lon_step:lon_end;
            % Adjust the size of lat and lon vectors to match the image dimensions
            lat       = lat(1:lat_size)+(lat_step/2);
            lon       = lon(1:lon_size)+(lon_step/2);
            % Assign latitude and longitude values to the corresponding variables
            netcdf.putVar(ncid_max,latid,lat);
            netcdf.putVar(ncid_max,lonid,lon);
        end
    else
        selectedImages = NaN(size(sortedDates{1,2}, 1),1);
        %resultImages   = NaN(size(sortedDates, 1), 1);
        imagesSynAll = NaN(size(sortedDates,1),1);
        map          = imagesSynAll;
        varMap       = imagesSynAll;
        availablePix = imagesSynAll;
        varianceBS   = imagesSynAll;
    end

    for rowIndex = 1:size(sortedDates,1)
        if bootstrap == true
            if bsSaveAll == true
                outputDirBootstrap = fullfile(outputDir, 'bootstrap', string(sortedDates(rowIndex,1)));
                if ~exist(outputDirBootstrap,'dir')
                    mkdir(outputDirBootstrap)
                end
            end
            % Find the index of the current image in the Dates variable
            [~, dateIndex] = ismember(sortedDates{rowIndex,2},learningDatesDate);
            % Select the K best image from the Learning dataset and add it to selectedImages
            for imageIndex = 1:nbImages %length(sortedDates{rowIndex,2})
                if nbImages ~= length(sortedDates{rowIndex,2}) && imageIndex == 1
                    warning(['nbImages .ne. number of available analogues (' num2str(nbImages) ' vs ' num2str(length(sortedDates{rowIndex,2})) ')'])
                end
                selectedImages(:,:,imageIndex) = learningData{dateIndex(imageIndex)};
            end
            % Calculate either the mode or the mean of the selected images
            if generationType == 1
                % Calculate the mode and save it to resultImages
                resultImages = mode(selectedImages,3);
            elseif generationType == 2
                % Calculate the mean and save it to resultImages
                selectedDist = 1./sortedDates{rowIndex,3}(1:nbImages);
                % Normalize the selectedDist values
                normalizedWeights = selectedDist / sum(selectedDist);
                % Perform element-wise multiplication with the weights
                weightedImages = bsxfun(@times, selectedImages, reshape(normalizedWeights, 1, 1, nbImages)); %length(sortedDates{rowIndex,2})
                varMap(:,:,rowIndex) = var(selectedImages,normalizedWeights,3);
                resultImages = sum(weightedImages,3);
            elseif generationType == 3
                % Calculate the mean and save it to resultImages
                resultImages = mean(selectedImages,3);
            elseif generationType == 4
                % Calculate the median and save it to resultImages
                resultImages = median(selectedImages,3);
            elseif generationType == 5 
                % Perform MPS simulation to generate result
                cellSelImages = cell(1, size(selectedImages, 3));
                selectedImages(isnan(selectedImages)) = nanValue;
                for img = 1:size(selectedImages, 3)
                    cellSelImages{img} = cat(3, selectedImages(:, :, img), loc);
                end
                resultImages = g2s('-a','qs', ...
                                   '-di',mask, ...
                                   '-ti',cellSelImages, ...
                                   '-ki',ki, ...
                                   '-dt',mps.dataType, ...
                                   '-k',mps.kValue, ...
                                   '-n',mps.neighbours, ...
                                   '-j',mps.processPwr);
                resultImages = resultImages(:,:,1);
                resultImages(resultImages==nanValue) = nan;
            else
                error('Generation type not defined!')
            end
            % Calculate the count of non-NaN values
            availablePix(:,:,rowIndex) = sum(~isnan(weightedImages), 3);
            if bsSaveAll == true
                % Write the resulting image to a GeoTIFF file
                outputBaseName = string(sortedDates(rowIndex,1)) + '_' + targetVarL(i) + '.tif';
                fullDestinationFileName = fullfile(outputDir, 'datesAll', outputBaseName);
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
            end
            map(:,:,rowIndex) = resultImages;
            resultImages(isnan(resultImages)) = -999;
            % bootstrap
            resultImagesBS     = NaN(imgLength, imgWidth, ensemble);
            %invDistance      = 1 ./ sortedDates{rowIndex,3};
            %bootstrapWeights = normalize(invDistance,'range',[0.1 1]); % normalise distance (3) / std (4) to [0.1 1]
            %bootstrapWeights = invDistance/sum(invDistance);
            for bs = 1:ensemble
                %bootstrapDates = randsample(sortedDates{rowIndex,2},numel(sortedDates{rowIndex,2}),true,bootstrapWeights);
                %bootstrapDates = randsample(sortedDates{rowIndex,2},numel(sortedDates{rowIndex,2}),true);
                bootstrapDates = randsample(sortedDates{rowIndex,2}(1:nbImages),nbImages,true);
                % Find the index of the current image in the Dates variable
                [~, dateIndex] = ismember(bootstrapDates,learningDatesDate);
                [~, distIndex] = ismember(bootstrapDates,sortedDates{rowIndex,2}(1:nbImages));
                % Select the K best image from the Learning dataset and add it to selectedImages
                for imageIndex = 1:nbImages %length(sortedDates{rowIndex,2})
                    if nbImages ~= length(sortedDates{rowIndex,2}) && imageIndex == 1
                        warning(['nbImages .ne. number of available analogues (' num2str(nbImages) ' vs ' num2str(length(sortedDates{rowIndex,2})) ')'])
                    end
                    selectedImages(:,:,imageIndex) = learningData{dateIndex(imageIndex)};
                end
                selectedDist = 1./sortedDates{rowIndex,3}(distIndex);
                % Normalize the selectedDist values
                normalizedWeights = selectedDist / sum(selectedDist);
                % Perform element-wise multiplication with the weights
                weightedImages = bsxfun(@times, selectedImages, reshape(normalizedWeights, 1, 1, nbImages)); %length(sortedDates{rowIndex,2})
                % Calculate either the mode or the mean of the selected images
                if generationType == 1
                    % Calculate the mode and save it to resultImagesBS
                    resultImagesBS(:,:,bs) = mode(selectedImages,3);
                elseif generationType == 2
                    % Calculate the mean and save it to resultImagesBS
                    resultImagesBS(:,:,bs) = sum(weightedImages,3);
                elseif generationType == 3
                    % Calculate the mean and save it to resultImagesBS
                    resultImagesBS(:,:,bs) = mean(selectedImages,3);
                elseif generationType == 4
                    % Calculate the median and save it to resultImagesBS
                    resultImagesBS(:,:,bs) = median(selectedImages,3);
                elseif generationType == 5
                    % Perform MPS simulation to generate result
                    cellSelImages = cell(1, size(selectedImages, 3));
                    selectedImages(isnan(selectedImages)) = nanValue;
                    for img = 1:size(selectedImages, 3)
                        cellSelImages{img} = cat(3, selectedImages(:, :, img), loc);
                    end
                    resultImages = g2s('-a','qs', ...
                                       '-di',mask, ...
                                       '-ti',cellSelImages, ...
                                       '-ki',ki, ...
                                       '-dt',mps.dataType, ...
                                       '-k',mps.kValue, ...
                                       '-n',mps.neighbours, ...
                                       '-j',mps.processPwr);
                    resultImages = resultImages(:,:,1);
                    resultImages(resultImages==nanValue) = nan;
                    resultImagesBS(:,:,bs) = resultImages;
                else
                    error('Generation type not defined!')
                end
            end
            % Calculate the count of non-NaN values
            availablePix(:,:,rowIndex) = sum(~isnan(selectedImages), 3);
            % Compute variance per pixel
            varianceBS(:,:,rowIndex) = var(resultImagesBS, 0, 3);
            % Compute mean of each day to determine quantile
            dayAvg = squeeze(mean(mean(resultImagesBS,'omitnan'),'omitnan'));
            dayAvg = sortrows([dayAvg (1:ensemble)']);
            %resultImagesMean = mean(resultImagesBS,3);
            % Store all bs days sorted according to mean of each day
            %imagesSynAll{rowIndex} = resultImagesBS;
            imagesSynAll{rowIndex} = resultImagesBS(:,:,dayAvg(:,2));
            bsMin = single(imagesSynAll{rowIndex}(:,:,1));
            bsMax = single(imagesSynAll{rowIndex}(:,:,end));
            bsMin(isnan(bsMin)) = -999;
            bsMax(isnan(bsMax)) = -999;
            if bsSaveAll == true
                for bs = 1:ensemble
                    % Write the resulting image to a GeoTIFF file
                    outputBaseName = string(sortedDates(rowIndex,1)) + '_' + num2str(bs) + '_' + targetVarL(i) + '.tif';
                    fullDestinationFileName = fullfile(outputDirBootstrap, outputBaseName);
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
                        t.write(single(imagesSynAll{rowIndex}(:,:,bs)));
                        t.close();
                    else
                        geotiffwrite(fullDestinationFileName,single(imagesSynAll{rowIndex}(:,:,bs)),GeoRef,'TiffTags',struct('Compression',Tiff.Compression.None));
                    end
                end
            end
            % Save min, deterministic and max in netCDF
            % Assign date
            dateStr  = convertStringsToChars(string(sortedDates{rowIndex, 1}));
            yearStr  = dateStr(1:4);
            monthStr = dateStr(5:6);
            dayStr   = dateStr(7:8);
            dateStrFormatted = [yearStr '-' monthStr '-' dayStr];
            % Write data for each date as a new time step along the 'time' dimension
            time = datenum(dateStrFormatted, 'yyyy-mm-dd');
            netcdf.putVar(ncid_min, timeid, rowIndex - 1, 1, time - 719529); % 719529 = 1970-01-01
            netcdf.putVar(ncid_det, timeid, rowIndex - 1, 1, time - 719529); % 719529 = 1970-01-01
            netcdf.putVar(ncid_max, timeid, rowIndex - 1, 1, time - 719529); % 719529 = 1970-01-01
            % Write data to the variable (hydrological map) for the current date
            ncwrite(fullDestinationFileNameMin, targetVarL(i), bsMin', [1, 1, rowIndex]); % <-----------------------------------------------------------------------------------
            ncwrite(fullDestinationFileNameDet, targetVarL(i), resultImages', [1, 1, rowIndex]);
            ncwrite(fullDestinationFileNameMax, targetVarL(i), bsMax', [1, 1, rowIndex]); % <-----------------------------------------------------------------------------------

            % Write the resulting image to a GeoTIFF file
            %outputBaseName = string(sortedDates(rowIndex,1)) + '_bsMean.tif';
            %fullDestinationFileName = fullfile(outputDir, var_low(i), outputBaseName);
            %disp(['  Downlading image ' num2str(rowIndex) '/' num2str(size(sortedDates,1))])
            %if isempty(GeoRef)
            %disp('    Georeferencing files missing! Unreferenced output...')
            %t = Tiff(fullDestinationFileName, 'w');
            %tagstruct.ImageLength         = imgLength;
            %tagstruct.ImageWidth          = imgWidth;
            %tagstruct.Compression         = Tiff.Compression.None;
            %tagstruct.SampleFormat        = Tiff.SampleFormat.IEEEFP;
            %tagstruct.Photometric         = Tiff.Photometric.MinIsBlack;
            %tagstruct.BitsPerSample       = 32;
            %tagstruct.SamplesPerPixel     = 1;
            %tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
            %t.setTag(tagstruct);
            %t.write(single(resultImagesMean));
            %t.close();
            %else
            %geotiffwrite(fullDestinationFileName,single(resultImagesMean),GeoRef,'TiffTags',struct('Compression',Tiff.Compression.None));
            %end
        else
            % Find the index of the current image in the Dates variable
            [~, dateIndex] = ismember(sortedDates{rowIndex,2},learningDatesDate);
            % Select the K best image from the Learning dataset and add it to selectedImages
            for imageIndex = 1:nbImages %length(sortedDates{rowIndex,2})
                if nbImages ~= length(sortedDates{rowIndex,2}) && imageIndex == 1
                    warning(['nbImages .ne. number of available analogues (' num2str(nbImages) ' vs ' num2str(length(sortedDates{rowIndex,2})) ')'])
                end
                if targetDim ~= 1
                    selectedImages(:,:,imageIndex) = learningData{dateIndex(imageIndex)};
                else
                    selectedImages(imageIndex) = learningData(dateIndex(imageIndex));
                end
            end
            % Calculate either the mode or the mean of the selected images
            if generationType == 1
                % Calculate the mode and save it to resultImages
                if targetDim ~= 1
                    resultImages = mode(selectedImages,3);
                else
                    resultImages = mode(selectedImages);
                end
            elseif generationType == 2
                % Calculate the mean and save it to resultImages
                sortedDates{rowIndex,3}(sortedDates{rowIndex,3} == 0) = 0.0001;
                selectedDist = 1./sortedDates{rowIndex,3}(1:nbImages);
                % Normalize the selectedDist values
                normalizedWeights = selectedDist / sum(selectedDist);
                % Perform element-wise multiplication with the weights
                if targetDim ~= 1
                    weightedImages = bsxfun(@times, selectedImages, reshape(normalizedWeights, 1, 1, nbImages)); %length(sortedDates{rowIndex,2})
                    varMap(:,:,rowIndex) = var(selectedImages,normalizedWeights,3);
                    resultImages = sum(weightedImages,3);
                else
                    weightedImages = selectedImages .* normalizedWeights;
                    varMap(rowIndex) = var(selectedImages,normalizedWeights);
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
            elseif generationType == 5
                % Perform MPS simulation to generate result
                if targetDim ~= 1
                    cellSelImages = cell(1, size(selectedImages, 3));
                    selectedImages(isnan(selectedImages)) = nanValue;
                    for img = 1:size(selectedImages, 3)
                        cellSelImages{img} = cat(3, selectedImages(:, :, img), loc);
                    end
                    resultImages = g2s('-a','qs', ...
                                       '-di',mask, ...
                                       '-ti',cellSelImages, ...
                                       '-ki',ki, ...
                                       '-dt',mps.dataType, ...
                                       '-k',mps.kValue, ...
                                       '-n',mps.neighbours, ...
                                       '-j',mps.processPwr);
                    resultImages = resultImages(:,:,1);
                    resultImages(resultImages==nanValue) = nan;
                else
                    error('MPS not implemented for 1D data generation...')
                end
            else
                error('Generation type not defined!')
            end
            if targetDim ~= 1
                map(:,:,rowIndex) = resultImages;
                % Calculate the count of non-NaN values
                availablePix(:,:,rowIndex) = sum(~isnan(selectedImages), 3);
                if outputType == 1
                    % Write the resulting image to a GeoTIFF file
                    outputBaseName = string(sortedDates(rowIndex,1)) + targetVarL(i) + '.tif';
                    fullDestinationFileName = fullfile(outputDir, targetVarL(i), outputBaseName);
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
                elseif outputType == 2
                    % Assign date
                    dateStr  = convertStringsToChars(string(sortedDates{rowIndex, 1}));
                    yearStr  = dateStr(1:4);
                    monthStr = dateStr(5:6);
                    dayStr   = dateStr(7:8);
                    dateStrFormatted = [yearStr '-' monthStr '-' dayStr];
                    % Write data for each date as a new time step along the 'time' dimension
                    time = datenum(dateStrFormatted, 'yyyy-mm-dd');
                    netcdf.putVar(ncid, timeid, rowIndex - 1, 1, time - 719529); % 719529 = 1970-01-01
                    % Write data to the variable for the current date
                    resultImages(isnan(resultImages)) = -999;
                    ncwrite(fullDestinationFileName, targetVarL(i), single(resultImages)', [1, 1, rowIndex]);
                else
                    error('Unknown output type. Choose 1 for GeoTiff or 2 for NetCDF...')
                end
            else
                % For 1D data
                map(rowIndex) = resultImages;
            end
        end
        if optimisation == false
            % Display computation progress
            progress = (100*(rowIndex/size(sortedDates,1)));
            fprintf(1,'\b\b\b\b%3.0f%%',progress);
        end
    end
    if outputType == 2 && bootstrap == false && targetDim ~= 1
        % Close the main netCDF file after the loop
        netcdf.close(ncid);
    elseif bootstrap == true
        netcdf.close(ncid_min);
        netcdf.close(ncid_det);
        netcdf.close(ncid_max);
    end

    % Save text file
    if targetDim == 1
        % Write the resulting image to a GeoTIFF file
        outputBaseName = strcat(targetVarL(i), '.txt');
        fullDestinationFileName = fullfile(outputDir, outputBaseName);
        Dates = cell2mat(sortedDates(:,1));
        Discharge = map;
        T = table(Dates, Discharge);
        writetable(T,fullDestinationFileName,'Delimiter','\t');
    end

    % Data outputs
    if i == 1
        synImages.date = cell2mat(sortedDates(:,1));
        fprintf('\n')
    end
    synImages.(targetVarL(i)) = map;
    varDist = strcat(targetVarL(i), "_Distances");
    minDist = single(nan(size(sortedDates,1),1));
    for c = 1:size(sortedDates, 1)
        values = sortedDates{c,3};
        minDist(c) = min(values);
    end
    %synImages.(varDist) = minDist;
    synImages.(varDist) = sortedDates(:,3);
    varPix = strcat(targetVarL(i), "_AvailablePixels");
    synImages.(varPix) = (availablePix./nbImages).*100;
    varName = strcat(targetVarL(i), "_Variance");
    synImages.(varName) = varMap;
    if bootstrap == true
        varBS = strcat(targetVarL(i), "_Bootstrap");
        BSvar = strcat(varBS, "Variance");
        synImages.(varBS) = imagesSynAll;
        synImages.(BSvar) = varianceBS;
    end
end

if optimisation == false && validation == true
    %fprintf('\n')
    disp('Saving synValidation.mat file...')
    save(fullfile(outputDir,'synValidation.mat'),'synImages', '-v7.3','-nocompression');
end

end
