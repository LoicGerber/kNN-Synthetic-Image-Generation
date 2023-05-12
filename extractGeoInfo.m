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
lat = netcdf.getVar(ncid,latVarID);
lat = unique(lat);
lon = netcdf.getVar(ncid,lonVarID);
lon = unique(lon);
[minLat,maxLat] = bounds(lat);
[minLon,maxLon] = bounds(lon);
ylimits         = [minLat maxLat];
xlimits         = [minLon maxLon];
rasterSize      = [length(lat) length(lon)];
GeoRef = georefcells(ylimits,xlimits,rasterSize, ...
    'ColumnsStartFrom','north','RowsStartFrom','west');
GeoRef.GeographicCRS = geocrs(coordRefSysCode);
save(fullfile(inputDir,'GeoRef.mat'), 'GeoRef');

end
