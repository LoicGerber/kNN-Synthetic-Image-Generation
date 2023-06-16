function visualiseMetrics(refValidation,synImages,validationMetric,metric,LdateStart,LdateEnd,outputDir)

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
gifVal = fullfile(outputDir,'validation.gif');
gifErr = fullfile(outputDir,'error.gif');

refDates = refValidation.date;
refData  = refValidation.maps;
synDates = synImages.date;
synData  = synImages.maps;

% Create an empty figure
figure(2);

% Loop over each file in the synthetic directory and find the corresponding
% file in the reference directory
for i = 1:size(synData,3)
    % Find the corresponding file in the reference directory with the same name
    referenceIndex = find(refDates == synDates(i));
    
    % If a matching file is found, display the two images side by side
    if ~isempty(referenceIndex)
        % Load the two images
        synthetic = synData(:,:,i);
        synthetic(synthetic==-999) = NaN;
        reference = refData(:,:,referenceIndex);
        reference(reference==-999) = NaN;

        % Create a figure with two subplots
        subplot(1,2,1);
        img1 = imshow(synthetic);
        colormap(gca, jet(256));
        set(img1, 'AlphaData', ~isnan(synthetic))
        %caxis([0 maxColor])
        caxis([0 6])
        title('Synthetic');
        %colorbar(gca,'southoutside')

        subplot(1,2,2);
        img2 = imshow(reference);
        colormap(gca, jet(256));
        set(img2, 'AlphaData', ~isnan(synthetic))
        %caxis([0 maxColor])
        caxis([0 6])
        title('Reference');
        h = colorbar(gca,'southoutside');

        % Add a colorbar to the reference image subplot
        %h = colorbar('southoutside');
        set(h, 'Position', [0.25 0.1 0.5 0.04]);
        set(get(h,'label'),'string','mm/day');

        % Set the title of the figure to the name of the images
        if metric == 1
            sgtitle({['{\bf\fontsize{14}' num2str(synDates(i)) '}'], ...
                ['{\fontsize{13}' 'RMSE: ' num2str(validationMetric(i,2),'%1.5f') '}']})
        elseif metric == 2
            sgtitle({['{\bf\fontsize{14}' num2str(synDates(i)) '}'], ...
                ['{\fontsize{13}' 'SPEM: ' num2str(validationMetric(i,2),'%1.5f') '}']})
        elseif metric == 3
            sgtitle({['{\bf\fontsize{14}' num2str(synDates(i)) '}'], ...
                ['{\fontsize{13}' 'SPAEF: ' num2str(validationMetric(i,2),'%1.5f') '}']})
        elseif metric == 4
            sgtitle({['{\bf\fontsize{14}' num2str(synDates(i)) '}'], ...
                ['{\fontsize{13}' 'SPOMF: ' num2str(validationMetric(i,2),'%1.5f') '}']})
        end

        % Save the current frame as a GIF
        set(gcf, 'color', 'white');
        frame = getframe(gcf);
        im = frame2im(frame);
        [imind, cm] = rgb2ind(im, 256);
        if i == 1
            imwrite(imind, cm, gifVal, 'gif', 'Loopcount', size(synData,3), 'DelayTime', 0.1);
        else
            imwrite(imind, cm, gifVal, 'gif', 'WriteMode', 'append', 'DelayTime', 0.1);
        end
    end
end

% Create an empty figure
figure(3);
% Loop over each file in the synthetic directory and find the corresponding
% file in the reference directory
for i = 1:size(synData,3)
    % Find the corresponding file in the reference directory with the same name
    referenceIndex = find(refDates == synDates(i));
    
    % If a matching file is found, display the two images side by side
    if ~isempty(referenceIndex)
        % Load the two images
        synthetic = synData(:,:,i);
        synthetic(synthetic==-999) = NaN;
        reference = refData(:,:,referenceIndex);
        reference(reference==-999) = NaN;

        % ERROR MAP
        subplot(1,1,1);
        errMap = imshow(reference-synthetic);
        set(errMap, 'AlphaData', ~isnan(synthetic))
        %caxis([0 6])
        title('Absolute error');
        h_err = colorbar(gca,'eastoutside');
        set(get(h_err,'label'),'string','mm/day');

        % Set the title of the figure to the name of the images
        if metric == 1
            sgtitle({['{\bf\fontsize{14}' num2str(synDates(i)) '}'], ...
                ['{\fontsize{13}' 'RMSE: ' num2str(validationMetric(i,2),'%1.5f') '}']})
        elseif metric == 2
            sgtitle({['{\bf\fontsize{14}' num2str(synDates(i)) '}'], ...
                ['{\fontsize{13}' 'SPEM: ' num2str(validationMetric(i,2),'%1.5f') '}']})
        elseif metric == 3
            sgtitle({['{\bf\fontsize{14}' num2str(synDates(i)) '}'], ...
                ['{\fontsize{13}' 'SPAEF: ' num2str(validationMetric(i,2),'%1.5f') '}']})
        elseif metric == 4
            sgtitle({['{\bf\fontsize{14}' num2str(synDates(i)) '}'], ...
                ['{\fontsize{13}' 'SPOMF: ' num2str(validationMetric(i,2),'%1.5f') '}']})
        end

        % Save the current frame as a GIF
        set(gcf, 'color', 'white');
        frame = getframe(gcf);
        im = frame2im(frame);
        [imind, cm] = rgb2ind(im, 256);
        if i == 1
            imwrite(imind, cm, gifErr, 'gif', 'Loopcount', size(synData,3), 'DelayTime', 0.1);
        else
            imwrite(imind, cm, gifErr, 'gif', 'WriteMode', 'append', 'DelayTime', 0.1);
        end
    end
end

end
