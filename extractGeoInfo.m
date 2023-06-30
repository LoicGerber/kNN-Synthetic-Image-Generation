function geoRef = extractGeoInfo(var,coordRefSysCode,rawDir,inputDir)

%
%
%
% REDO DOCUMENTATION
%
%
%

% Create R matrix containing InputVar product informations
% Open the netCDF file
for i = 1:numel(var)
    ncid = netcdf.open(fullfile(rawDir,var(i)+'.nc'),'NOWRITE');
    % Get the latitude and longitude variables
    try
        latVarID = netcdf.inqVarID(ncid,'lat');
        lonVarID = netcdf.inqVarID(ncid,'lon');
    catch
        latVarID = netcdf.inqVarID(ncid,'latitude');
        lonVarID = netcdf.inqVarID(ncid,'longitude');
    end
    % Read the latitude and longitude data
    lat       = netcdf.getVar(ncid,latVarID);
    lat       = unique(lat);
    lat_space = (max(lat)-min(lat))/(numel(lat)-1);
    lat       = lat - lat_space/2;
    lon       = netcdf.getVar(ncid,lonVarID);
    lon       = unique(lon);
    lon_space = (max(lon)-min(lon))/(numel(lon)-1);
    lon       = lon - lon_space/2;
    [minLat,maxLat] = bounds(lat);
    [minLon,maxLon] = bounds(lon);
    ylimits         = double([minLat maxLat+lat_space]);
    xlimits         = double([minLon maxLon+lon_space]);
    rasterSize      = double([length(lat) length(lon)]);
    GeoRef = georefcells(ylimits,xlimits,rasterSize, ...
        'ColumnsStartFrom','north','RowsStartFrom','west');
    if i == 1
        GeoRef.CellExtentInLatitude  = double((maxLat-minLat)/(numel(lat)-1));
        GeoRef.CellExtentInLongitude = double((maxLon-minLon)/(numel(lon)-1));
    end
    GeoRef.GeographicCRS = geocrs(coordRefSysCode);
    geoRef.(var(i)) = GeoRef;
end

disp(strcat('Saving geoRef data...'))
save(fullfile(inputDir,strcat('geoRef.mat')), 'geoRef');

end
