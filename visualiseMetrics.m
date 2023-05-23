function visualiseMetrics(validationMetric,metric,LdateStart,LdateEnd,outputDir)

%
%
%
% REDO DOCUMENTATION
%
%
%

figure(1)
plot(datetime(validationMetric(:,1),'ConvertFrom','yyyyMMdd','Format','dd/MM/yyyy'), ...
    validationMetric(:,2));
yline(mean(validationMetric(:,2)),'-',['Mean: ' num2str(mean(validationMetric(:,2)))],'Color','r')
ylim([0 1.4])
if metric == 1
    title('RMSE')
elseif metric == 2
    title('SPEM')
elseif metric == 3
    title('SPAEF')
elseif metric == 4
    title('SPOMF')
end
subtitle(['Learning periode: ' char(datetime(LdateStart,'ConvertFrom','yyyyMMdd','Format','dd/MM/yyyy')) ...
    ' - ' char(datetime(LdateEnd,'ConvertFrom','yyyyMMdd','Format','dd/MM/yyyy'))])
xlabel('Date')
if metric == 1
    ylabel('RMSE')
elseif metric == 2
    ylabel('SPEM')
elseif metric == 3
    ylabel('SPAEF')
elseif metric == 4
    ylabel('SPOMF')
end
set(gcf, 'color', 'white');
grid on
saveas(gcf,strcat(outputDir,'validation.png'))


% Set the output GIF file name
gifName = fullfile(outputDir,'validation.gif');

% Specify the directories containing the images
syntheticDir = fullfile(outputDir,'syntheticImages/');
referenceDir = fullfile(outputDir,'referenceImages/');

% Get a list of all the TIFF files in each directory
syntheticFiles = dir(fullfile(syntheticDir, '*.tif'));
referenceFiles = dir(fullfile(referenceDir, '*.tif'));

% Create an empty figure
figure(2);

% Loop over each file in the synthetic directory and find the corresponding
% file in the reference directory
for i = 1:numel(syntheticFiles)
    % Extract the filename and extension of the current synthetic file
    [~, syntheticName, syntheticExt] = fileparts(syntheticFiles(i).name);

    % Find the corresponding file in the reference directory with the same name
    referenceIndex = find(strcmp({referenceFiles.name}, [syntheticName syntheticExt]));

    % If a matching file is found, display the two images side by side
    if ~isempty(referenceIndex)
        % Load the two images
        synthetic = imread(fullfile(syntheticDir, syntheticFiles(i).name));
        synthetic(synthetic==-999) = NaN;
        reference = imread(fullfile(referenceDir, referenceFiles(referenceIndex).name));
        reference(reference==-999) = NaN;
        %minColor = min(min(reference(:),min(synthetic(:))));
        %maxColor = max(max(reference(:),max(synthetic(:))));

        % Create a figure with two subplots
        subplot(1,3,1);
        img1 = imshow(synthetic);
        colormap(gca, jet(256));
        set(img1, 'AlphaData', ~isnan(synthetic))
        %caxis([0 maxColor])
        caxis([0 6])
        title('Synthetic');
        %colorbar(gca,'southoutside')

        subplot(1,3,2);
        img2 = imshow(reference);
        colormap(gca, jet(256));
        set(img2, 'AlphaData', ~isnan(synthetic))
        %caxis([0 maxColor])
        caxis([0 6])
        title('Reference');
        h = colorbar(gca,'southoutside');
        
        % !!! ADD ERROR MAP !!!
        subplot(1,3,3);
        img3 = imshow(abs(reference-synthetic));
        set(img3, 'AlphaData', ~isnan(synthetic))
        %caxis([0 6])
        title('Absolute error');
        h_err = colorbar(gca,'eastoutside');
        set(get(h_err,'label'),'string','mm/day');
        % !!! ADD ERROR MAP !!!
        
        % Add a colorbar to the reference image subplot
        %h = colorbar('southoutside');
        set(h, 'Position', [0.25 0.1 0.5 0.04]);
        set(get(h,'label'),'string','mm/day');

        % Set the title of the figure to the name of the images
        if metric == 1
            sgtitle({['{\bf\fontsize{14}' syntheticName '}'], ...
                ['{\fontsize{13}' 'RMSE: ' num2str(validationMetric(i,2),'%1.5f') '}']})
        elseif metric == 2
            sgtitle({['{\bf\fontsize{14}' syntheticName '}'], ...
                ['{\fontsize{13}' 'SPEM: ' num2str(validationMetric(i,2),'%1.5f') '}']})
        elseif metric == 3
            sgtitle({['{\bf\fontsize{14}' syntheticName '}'], ...
                ['{\fontsize{13}' 'SPAEF: ' num2str(validationMetric(i,2),'%1.5f') '}']})
        elseif metric == 4
            sgtitle({['{\bf\fontsize{14}' syntheticName '}'], ...
                ['{\fontsize{13}' 'SPOMF: ' num2str(validationMetric(i,2),'%1.5f') '}']})
        end

        % Save the current frame as a GIF
        set(gcf, 'color', 'white');
        frame = getframe(gcf);
        im = frame2im(frame);
        [imind, cm] = rgb2ind(im, 256);
        if i == 1
            imwrite(imind, cm, gifName, 'gif', 'Loopcount', numel(syntheticFiles), 'DelayTime', 0.1);
        else
            imwrite(imind, cm, gifName, 'gif', 'WriteMode', 'append', 'DelayTime', 0.1);
        end

        % Pause for a short time to create animation effect
        %pause(0.05);
    end
end

end


