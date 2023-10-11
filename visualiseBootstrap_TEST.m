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

meanVar = mean(Et_3D,3,'omitnan');
sumVar = sum(Et_3D,3,'omitnan');

colormap jet

figure
figMean = imagesc(meanVar);
set(figMean, 'AlphaData', ~isnan(meanVar))
title('Mean variance')
set(gcf, 'color', 'white');
colorbar
axis equal
axis off

figure
figSum = imagesc(sumVar);
set(figSum, 'AlphaData', ~(sumVar==0))
title('Variance sum')
set(gcf, 'color', 'white');
colorbar
axis equal
axis off
