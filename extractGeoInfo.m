function R = extractGeoInfo(var,coordRefSysCode,rawDir,inputDir)

% Create R matrix containing InputVar product informations
% Open the netCDF file
ncid = netcdf.open(fullfile(rawDir,var+'.nc'),'NOWRITE');
% Get the latitude and longitude variables
latvarid = netcdf.inqVarID(ncid,'lat');
lonvarid = netcdf.inqVarID(ncid,'lon');
% Read the latitude and longitude data
lat = netcdf.getVar(ncid,latvarid);
lat = unique(lat);
lon = netcdf.getVar(ncid,lonvarid);
lon = unique(lon);
[minLat,maxLat] = bounds(lat);
[minLon,maxLon] = bounds(lon);
ylimits         = [minLat maxLat];
xlimits         = [minLon maxLon];
rasterSize      = [length(lat) length(lon)];
R = georefcells(xlimits,ylimits,rasterSize, ...
    'ColumnsStartFrom','north','RowsStartFrom','west');
R.GeographicCRS = geocrs(coordRefSysCode);
save(fullfile(inputDir,'R.mat'), 'R');

end
