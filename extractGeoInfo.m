function geoRef = extractGeoInfo(targetVar, coordRefSysCode, rawDir, inputDir)

for i = 1:numel(targetVar)
    % Check file extension and use appropriate function
    tarVarL = lower(targetVar(i));
    filePath = fullfile(rawDir,[convertStringsToChars(tarVarL) '.nc']);
    
    if exist(filePath, 'file') == 2 %strcmpi(ext, '.nc')
        ncid = netcdf.open(filePath, 'NOWRITE');
        latVarID = netcdf.inqVarID(ncid, 'lat');
        lonVarID = netcdf.inqVarID(ncid, 'lon');
        lat = netcdf.getVar(ncid, latVarID);
        lon = netcdf.getVar(ncid, lonVarID);
        netcdf.close(ncid);

        lat = unique(lat);
        lat_space = (max(lat) - min(lat)) / (numel(lat) - 1);
        lat = lat - lat_space / 2;
        lon = unique(lon);
        lon_space = (max(lon) - min(lon)) / (numel(lon) - 1);
        lon = lon - lon_space / 2;

        [minLat, maxLat] = bounds(lat);
        [minLon, maxLon] = bounds(lon);
        ylimits = double([minLat maxLat + lat_space]);
        xlimits = double([minLon maxLon + lon_space]);
        rasterSize = double([length(lat) length(lon)]);
        GeoRef = georefcells(ylimits, xlimits, rasterSize, ...
            'ColumnsStartFrom', 'north', 'RowsStartFrom', 'west');

        if i == 1
            GeoRef.CellExtentInLatitude = double((maxLat - minLat) / (numel(lat) - 1));
            GeoRef.CellExtentInLongitude = double((maxLon - minLon) / (numel(lon) - 1));
        end
        GeoRef.GeographicCRS = geocrs(coordRefSysCode);
    else
        tiffDir = fullfile(rawDir, targetVar(i));
        tiffFiles = dir(fullfile(tiffDir, '*.tif'));
        
        if isempty(tiffFiles)
            warning('No GeoTIFF files found in the directory.');
        end
        
        % Select the first GeoTIFF file in the directory
        file_path = fullfile(tiffDir, tiffFiles(1).name);
        info = geotiffinfo(file_path);
        GeoRef = info.SpatialRef;
        GeoRef.GeographicCRS = geocrs(coordRefSysCode);
    end
    geoRef.(targetVar(i)) = GeoRef;
end

disp('  Saving geoRef data...')
save(fullfile(inputDir, 'geoRef.mat'), 'geoRef');

end
