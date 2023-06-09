function visualiseMetrics(var,refValidation,synImages,validationMetric,metric,LdateStart,LdateEnd,QdateStart,QdateEnd,bootstrap,outputDir)

%
%
%
% REDO DOCUMENTATION
%
%
%

startLdate = char(datetime(LdateStart,'ConvertFrom','yyyyMMdd','Format','dd/MM/yyyy'));
startQdate = char(datetime(QdateStart,'ConvertFrom','yyyyMMdd','Format','dd/MM/yyyy'));
endLdate   = char(datetime(LdateEnd,'ConvertFrom','yyyyMMdd','Format','dd/MM/yyyy'));
endQdate   = char(datetime(QdateEnd,'ConvertFrom','yyyyMMdd','Format','dd/MM/yyyy'));

for k = 1:numel(var)
    if bootstrap == true
        dates = datetime(cell2mat(validationMetric.(var(k))(:,1)),'ConvertFrom','yyyyMMdd','Format','dd/MM/yyyy');
        ensembleData = validationMetric.(var(k))(:,2);
        singleData = cell2mat(validationMetric.(var(k))(:,3));
        maxValues = zeros(size(validationMetric.(var(k)), 1), 1);
        meanValues = zeros(size(validationMetric.(var(k)), 1), 1);
        minValues = zeros(size(validationMetric.(var(k)), 1), 1);
        for i = 1:size(validationMetric.(var(k)), 1)
            values = sort(validationMetric.(var(k)){i, 2});  % Extract values from the second column of the cell
            maxValues(i) = values(end);  % Compute the maximum value
            meanValues(i) = mean(values);
            minValues(i) = values(1);  % Compute the minimum value
        end
        figure(k)
        hold on
        inBetweenRegionX = [dates', fliplr(dates')];
        inBetweenRegionY = [maxValues', fliplr(minValues')];
        patch(inBetweenRegionX, inBetweenRegionY, 'r', 'LineStyle', 'none', 'FaceAlpha', 0.5)
        plot(dates, singleData, 'k-', 'LineWidth', 2)
        for i = 1:numel(dates)
            ensemble = ensembleData{i};
            for j = 1:numel(ensemble)
                plot(dates(i), ensemble(j), 'k.'); % Plot the value against the date
            end
        end
        hold off
        xtickangle(45); % Rotate x-axis labels for better readability
        %yline(mean(validationMetric.(var(k))(:,2)),'-',['Mean: ' num2str(mean(validationMetric.(var(k))(:,2)))],'Color','r')
        %ylim([0 1.4])
        if metric == 1
            title([convertStringsToChars(var(k)) ' - RMSE'])
        elseif metric == 2
            title([convertStringsToChars(var(k)) ' - SPEM'])
        elseif metric == 3
            title([convertStringsToChars(var(k)) ' - SPAEF'])
        elseif metric == 4
            title([convertStringsToChars(var(k)) ' - SPOMF'])
        end
        subtitle(['Learning periode: ' startLdate '-' startQdate ' - ' endQdate '-' endLdate])
        str = {['Mean RMSE: ' num2str(mean(singleData))], ['Mean ensemble RMSE: ' num2str(mean(meanValues))]};
        text(15,.8,str);
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
        saveas(gcf,strcat(outputDir,'bsValidation.png'))
    else
        figure(k)
        plot(datetime(validationMetric.(var(k))(:,1),'ConvertFrom','yyyyMMdd','Format','dd/MM/yyyy'), ...
            validationMetric.(var(k))(:,2));
        yline(mean(validationMetric.(var(k))(:,2)),'-',['Mean: ' num2str(mean(validationMetric.(var(k))(:,2)))],'Color','r')
        %ylim([0 1.4])
        if metric == 1
            title([convertStringsToChars(var(k)) ' - RMSE (mean: ' num2str(mean(validationMetric.(var(k))(:,2))) ')'])
        elseif metric == 2
            title([convertStringsToChars(var(k)) ' - SPEM (mean: ' num2str(mean(validationMetric.(var(k))(:,2))) ')'])
        elseif metric == 3
            title([convertStringsToChars(var(k)) ' - SPAEF (mean: ' num2str(mean(validationMetric.(var(k))(:,2))) ')'])
        elseif metric == 4
            title([convertStringsToChars(var(k)) ' - SPOMF (mean: ' num2str(mean(validationMetric.(var(k))(:,2))) ')'])
        end
        subtitle(['Learning periode: ' startLdate '-' startQdate ' - ' endQdate '-' endLdate])
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
        refData  = refValidation.(var(k));
        synDates = synImages.date;
        synData  = synImages.(var(k));

        % Create an empty figure
        figure(k+2);

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

                sgtitle(var(k))

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
                set(get(h,'label'),'string','Evapotranspiration [mm/day]');

                % Set the title of the figure to the name of the images
                if metric == 1
                    sgtitle({['{\bf\fontsize{14}' num2str(synDates(i)) '}'], ...
                        ['{\fontsize{13}' 'RMSE: ' num2str(validationMetric.(var(k))(i,2),'%1.5f') '}']})
                elseif metric == 2
                    sgtitle({['{\bf\fontsize{14}' num2str(synDates(i)) '}'], ...
                        ['{\fontsize{13}' 'SPEM: ' num2str(validationMetric.(var(k))(i,2),'%1.5f') '}']})
                elseif metric == 3
                    sgtitle({['{\bf\fontsize{14}' num2str(synDates(i)) '}'], ...
                        ['{\fontsize{13}' 'SPAEF: ' num2str(validationMetric.(var(k))(i,2),'%1.5f') '}']})
                elseif metric == 4
                    sgtitle({['{\bf\fontsize{14}' num2str(synDates(i)) '}'], ...
                        ['{\fontsize{13}' 'SPOMF: ' num2str(validationMetric.(var(k))(i,2),'%1.5f') '}']})
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
        figure(k+3);
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
                title([convertStringsToChars(var(k)) ' - Absolute error']);
                h_err = colorbar(gca,'eastoutside');
                set(get(h_err,'label'),'string','Evapotranspiration [mm/day]');

                % Set the title of the figure to the name of the images
                if metric == 1
                    sgtitle({['{\bf\fontsize{14}' num2str(synDates(i)) '}'], ...
                        ['{\fontsize{13}' 'RMSE: ' num2str(validationMetric.(var(k))(i,2),'%1.5f') '}']})
                elseif metric == 2
                    sgtitle({['{\bf\fontsize{14}' num2str(synDates(i)) '}'], ...
                        ['{\fontsize{13}' 'SPEM: ' num2str(validationMetric.(var(k))(i,2),'%1.5f') '}']})
                elseif metric == 3
                    sgtitle({['{\bf\fontsize{14}' num2str(synDates(i)) '}'], ...
                        ['{\fontsize{13}' 'SPAEF: ' num2str(validationMetric.(var(k))(i,2),'%1.5f') '}']})
                elseif metric == 4
                    sgtitle({['{\bf\fontsize{14}' num2str(synDates(i)) '}'], ...
                        ['{\fontsize{13}' 'SPOMF: ' num2str(validationMetric.(var(k))(i,2),'%1.5f') '}']})
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
end

end
