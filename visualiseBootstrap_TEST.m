Et_3D = synImages.Et_BootstrapVariance;
num_slices = size(Et_3D, 3);

for i = 1:num_slices
    figure('WindowState', 'maximized');
    hold on
    slice = Et_3D(:,:,i);
    fig = imagesc(slice);
    colormap(gca, jet(256));
    set(fig, 'AlphaData', ~isnan(slice))
    set(gcf, 'color', 'white');
    colorbar;
    caxis([0 max(slice(:))])
    axis off
    title(datestr(datetime(synImages.date(i), 'ConvertFrom', 'yyyymmdd'),'dd.mm.yyyy'));
    hold off
end
