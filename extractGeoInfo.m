function GeoRef = extractGeoInfo(var,coordRefSysCode,rawDir,inputDir)

%
%
%
% REDO DOCUMENTATION
%
%
%

% Create R matrix containing InputVar product informations
% Open the netCDF file
ncid = netcdf.open(fullfile(rawDir,var+'.nc'),'NOWRITE');
% Get the latitude and longitude variables
latVarID = netcdf.inqVarID(ncid,'lat'); 
lonVarID = netcdf.inqVarID(ncid,'lon'); 
% Read the latitude and longitude data
lat       = netcdf.getVar(ncid,latVarID);
lat       = unique(lat);
lat_space = (max(lat)-min(lat))/(numel(lat)-1);
lat       = lat - lat_space/2;
lon       = netcdf.getVar(ncid,lonVarID);
lon       = unique(lon);
lon_space = (max(lon)-min(lon))/(numel(lon)-1);
lon       = lon -lon_space/2;
[minLat,maxLat] = bounds(lat);
[minLon,maxLon] = bounds(lon);
ylimits         = [minLat maxLat+lat_space];
xlimits         = [minLon maxLon+lon_space];
rasterSize      = [length(lat) length(lon)];
GeoRef = georefcells(ylimits,xlimits,rasterSize, ...
    'ColumnsStartFrom','north','RowsStartFrom','west');
GeoRef.CellExtentInLatitude  = (maxLat-minLat)/(numel(lat)-1);
GeoRef.CellExtentInLongitude = (maxLon-minLon)/(numel(lon)-1);
GeoRef.GeographicCRS = geocrs(coordRefSysCode);
disp('Saving GeoRef data...')
save(fullfile(inputDir,'GeoRef.mat'), 'GeoRef');

end
