function visualiseMetrics(targetVar,refValidation,synImages,validationMetric,metricV,metricKNN,LdateStart,LdateEnd,QdateStart,QdateEnd,bootstrap,outputDir)

%
%
%
% REDO DOCUMENTATION
%
%
%

startLdate = char(datetime(LdateStart,'ConvertFrom','yyyyMMdd','Format','dd/MM/yyyy'));
startQdate = char(datetime(QdateStart,'ConvertFrom','yyyyMMdd','Format','dd/MM/yyyy')-days(1));
endLdate   = char(datetime(LdateEnd,'ConvertFrom','yyyyMMdd','Format','dd/MM/yyyy'));
endQdate   = char(datetime(QdateEnd,'ConvertFrom','yyyyMMdd','Format','dd/MM/yyyy')+days(1));

for k = 1:numel(targetVar)
    if bootstrap == true
        dates = datetime(cell2mat(validationMetric.(targetVar(k))(:,1)),'ConvertFrom','yyyyMMdd','Format','dd/MM/yyyy');
        % Validation data
        singleDataVal = cell2mat(validationMetric.(targetVar(k))(:,3));
        maxValuesVal  = zeros(size(validationMetric.(targetVar(k)), 1), 1);
        meanValuesVal = zeros(size(validationMetric.(targetVar(k)), 1), 1);
        minValuesVal  = zeros(size(validationMetric.(targetVar(k)), 1), 1);
        for i = 1:size(validationMetric.(targetVar(k)), 1)
            valuesVal = sort(validationMetric.(targetVar(k)){i, 2});  % Extract values from the second column of the cell
            maxValuesVal(i)  = valuesVal(end);  % Compute the maximum value
            meanValuesVal(i) = mean(valuesVal);
            minValuesVal(i)  = valuesVal(1);  % Compute the minimum value
        end
        % Synthetic data
        varBs = strcat(targetVar(k),'_Bootstrap');
        varBD = strcat(targetVar(k),'_BestDistance');
        refData    = squeeze(mean(mean(refValidation.(targetVar(k)),1,'omitnan'),2,'omitnan'));  % Extract mean of ref variable
        synData    = squeeze(mean(mean(synImages.(targetVar(k)),1,'omitnan'),2,'omitnan'));
        maxValues  = zeros(size(validationMetric.(targetVar(k)), 1), 1);
        meanValues = zeros(size(validationMetric.(targetVar(k)), 1), 1);
        minValues  = zeros(size(validationMetric.(targetVar(k)), 1), 1);
        for i = 1:size(validationMetric.(targetVar(k)), 1)
            synValues     = sort(squeeze(mean(mean(synImages.(varBs(k)){i},'omitnan'),'omitnan')));  % Extract mean of variable and sort
            maxValues(i)  = synValues(end);  % Compute the maximum value
            meanValues(i) = mean(synValues);
            minValues(i)  = synValues(1);  % Compute the minimum value
        end
        bestDist   = synImages.(varBD);
        
        % -----------------------------------------------------------------------------------------------------

        figure('WindowState', 'maximized');
        hold on
        inBetweenRegionX = [dates', fliplr(dates')];
        inBetweenRegionY = [maxValues', fliplr(minValues')];
        patch(inBetweenRegionX, inBetweenRegionY, 'k', 'LineStyle', 'none', 'FaceAlpha', 0.25)
        %plot(dates, meanValues, 'k--', 'LineWidth', 1)
        plot(dates, synData, 'k-', 'LineWidth', 1)
        plot(dates, refData, 'r-', 'LineWidth', 1)
        hold off
        title([convertStringsToChars(targetVar(k)) ' - MEAN'])
        r = corr(synData,refData);
        nseSynRef = 1-(sum((synData-refData).^2)/sum((synData-mean(synData)).^2));
        alpha = std(synData)/std(refData);
        beta  = mean(synData)/mean(refData);
        kgeSynRef = 1-(sqrt((r-1)^2 + (alpha-1)^2 + (beta-1)^2));
        str = {['Corr: ' num2str(r,'%.5f')] ['NSE: ' num2str(nseSynRef,'%.5f')] ['KGE: ' num2str(kgeSynRef,'%.5f')]};
        if strcmp(startLdate, startQdate)
            subtitle([['Learning periode: ' endQdate '-' endLdate] str])
        elseif strcmp(endQdate, endLdate)
            subtitle([['Learning periode: ' startLdate '-' startQdate] str])
        else
            subtitle([['Learning periode: ' startLdate '-' startQdate ' - ' endQdate '-' endLdate] str])
        end
        xlabel('Date')
        ylabel(strcat("Mean ", targetVar(k)))
        %legend('Synthetic data spread','Synthetic data mean','Deterministic mean','Reference data mean')
        legend('Synthetic data spread','Deterministic mean','Reference data mean')
        set(gcf, 'color', 'white');
        grid on
        saveas(gcf,strcat(outputDir,['bsValidation_AVG_' convertStringsToChars(targetVar(k)) '.png']))
        
        % -----------------------------------------------------------------------------------------------------------------------------

        figure('WindowState', 'maximized');
        hold on
        inBetweenRegionX = [dates', fliplr(dates')];
        inBetweenRegionY = [maxValuesVal', fliplr(minValuesVal')];
        patch(inBetweenRegionX, inBetweenRegionY, 'k', 'LineStyle', 'none', 'FaceAlpha', 0.25)
        %plot(dates, meanValuesVal, 'k--', 'LineWidth', 1)
        plot(dates, singleDataVal, 'k-', 'LineWidth', 1)
        yyaxis right
        plot(dates, bestDist, 'r', 'LineWidth', 1)
        if metricKNN == 1
            ylabel('RMSE')
        elseif metricKNN == 2
            ylabel('MAE')
        elseif metricKNN == 3
            ylabel('Manhattan distance')
        elseif metricKNN == 4
            ylabel('Euclidean distance')
        end
%         for i = 1:numel(dates)
%             ensemble = ensembleData{i};
%             for j = 1:numel(ensemble)
%                 plot(dates(i), ensemble(j), 'k.'); % Plot the value against the date
%             end
%         end
        hold off
        %xtickangle(45); % Rotate x-axis labels for better readability
        %yline(mean(validationMetric.(var(k))(:,2)),'-',['Mean: ' num2str(mean(validationMetric.(var(k))(:,2)))],'Color','r')
        %ylim([0 1.4])
        yyaxis left
        if metricV == 1
            title([convertStringsToChars(targetVar(k)) ' - RMSE'])
        elseif metricV == 2
            title([convertStringsToChars(targetVar(k)) ' - SPEM'])
        elseif metricV == 3
            title([convertStringsToChars(targetVar(k)) ' - SPAEF'])
        elseif metricV == 4
            title([convertStringsToChars(targetVar(k)) ' - SPOMF'])
        end
        %str = {['Mean RMSE: ' num2str(mean(singleDataVal),'%.5f')], ['Mean ensemble RMSE: ' num2str(mean(meanValuesVal),'%.5f')]};
        str = {['Mean RMSE: ' num2str(mean(singleDataVal),'%.5f')], ['RMSE - MAE correlation: ' num2str(corr(singleDataVal,bestDist),'%.5f')]};
        if strcmp(startLdate, startQdate)
            subtitle([['Learning periode: ' endQdate '-' endLdate] str])
        elseif strcmp(endQdate, endLdate)
            subtitle([['Learning periode: ' startLdate '-' startQdate] str])
        else
            subtitle([['Learning periode: ' startLdate '-' startQdate ' - ' endQdate '-' endLdate] str])
        end        
        xlabel('Date')
        if metricV == 1
            ylabel('RMSE')
        elseif metricV == 2
            ylabel('SPEM')
        elseif metricV == 3
            ylabel('SPAEF')
        elseif metricV == 4
            ylabel('SPOMF')
        end
        ax = gca;
        ax.YAxis(1).Color = 'k';
        ax.YAxis(2).Color = 'r';
        linkprop(ax.YAxis, 'Limits');
        %legend('Stochastic ensembles','Stochastic mean','Deterministic','Best distance')
        legend('Stochastic ensembles','Deterministic','Best distance')
        set(gcf, 'color', 'white');
        grid on
        saveas(gcf,strcat(outputDir,['bsValidation_RMSE_' convertStringsToChars(targetVar(k)) '.png']))
    else
        figure('WindowState', 'maximized');
        plot(datetime(validationMetric.(targetVar(k))(:,1),'ConvertFrom','yyyyMMdd','Format','dd/MM/yyyy'), ...
            validationMetric.(targetVar(k))(:,2));
        yline(mean(validationMetric.(targetVar(k))(:,2)),'-',['Mean: ' num2str(mean(validationMetric.(targetVar(k))(:,2)))],'Color','r')
        %ylim([0 1.4])
        if metricV == 1
            title([convertStringsToChars(targetVar(k)) ' - RMSE (mean: ' num2str(mean(validationMetric.(targetVar(k))(:,2))) ')'])
        elseif metricV == 2
            title([convertStringsToChars(targetVar(k)) ' - SPEM (mean: ' num2str(mean(validationMetric.(targetVar(k))(:,2))) ')'])
        elseif metricV == 3
            title([convertStringsToChars(targetVar(k)) ' - SPAEF (mean: ' num2str(mean(validationMetric.(targetVar(k))(:,2))) ')'])
        elseif metricV == 4
            title([convertStringsToChars(targetVar(k)) ' - SPOMF (mean: ' num2str(mean(validationMetric.(targetVar(k))(:,2))) ')'])
        end
        if strcmp(startLdate, startQdate)
            subtitle(['Learning periode: ' endQdate '-' endLdate])
        elseif strcmp(endQdate, endLdate)
            subtitle(['Learning periode: ' startLdate '-' startQdate])
        else
            subtitle(['Learning periode: ' startLdate '-' startQdate ' - ' endQdate '-' endLdate])
        end        
        xlabel('Date')
        if metricV == 1
            ylabel('RMSE')
        elseif metricV == 2
            ylabel('SPEM')
        elseif metricV == 3
            ylabel('SPAEF')
        elseif metricV == 4
            ylabel('SPOMF')
        end
        set(gcf, 'color', 'white');
        grid on
        saveas(gcf,strcat(outputDir,['validation_' convertStringsToChars(targetVar(k)) '.png']))

        % --------------------------------------------------------------------

        refData    = squeeze(mean(mean(refValidation.(targetVar(k)),1,'omitnan'),2,'omitnan'));  % Extract mean of ref variable
        synData    = squeeze(mean(mean(synImages.(targetVar(k)),1,'omitnan'),2,'omitnan'));
        figure('WindowState', 'maximized');
        date = datetime(validationMetric.(targetVar(k))(:,1),'ConvertFrom','yyyyMMdd','Format','dd/MM/yyyy');
        plot(date, refData, 'r-', date, synData, 'k-');
        legend('Reference','Synthetic','Location','southeast')
        xlabel('Date')
        ylabel('Evaporation [mm/day]')
        title(['Mean ' convertStringsToChars(targetVar(k))])
        r = corr(synData,refData);
        nseSynRef = 1-(sum((synData-refData).^2)/sum((synData-mean(synData)).^2));
        alpha = std(synData)/std(refData);
        beta  = mean(synData)/mean(refData);
        kgeSynRef = 1-(sqrt((r-1)^2 + (alpha-1)^2 + (beta-1)^2));
        str = {['Corr: ' num2str(r,'%.5f')] ['NSE: ' num2str(nseSynRef,'%.5f')] ['KGE: ' num2str(kgeSynRef,'%.5f')]};
        subtitle(str)
        grid on
        box off
        %legend boxoff 
        set(gcf, 'color', 'white');
        saveas(gcf,strcat(outputDir,['correlation_' convertStringsToChars(targetVar(k)) '.png']))

        % --------------------------------------------------------------------

        refData    = squeeze(var(refValidation.(targetVar(k)),0,[1 2],'omitnan'));
        synData    = squeeze(var(synImages.(targetVar(k)),0,[1 2],'omitnan'));
        figure('WindowState', 'maximized');
        date = datetime(validationMetric.(targetVar(k))(:,1),'ConvertFrom','yyyyMMdd','Format','dd/MM/yyyy');
        plot(date, refData, 'r-', date, synData, 'k-');
        legend('Reference','Synthetic','Location','northeast')
        xlabel('Date')
        ylabel('Evaporation [mm/day]')
        title([convertStringsToChars(targetVar(k)) ' variance'])
        r = corr(synData,refData);
        nseSynRef = 1-(sum((synData-refData).^2)/sum((synData-mean(synData)).^2));
        alpha = std(synData)/std(refData);
        beta  = mean(synData)/mean(refData);
        kgeSynRef = 1-(sqrt((r-1)^2 + (alpha-1)^2 + (beta-1)^2));
        str = {['Corr: ' num2str(r,'%.5f')] ['NSE: ' num2str(nseSynRef,'%.5f')] ['KGE: ' num2str(kgeSynRef,'%.5f')]};
        subtitle(str)
        grid on
        box off
        %legend boxoff 
        set(gcf, 'color', 'white');
        saveas(gcf,strcat(outputDir,['variance_' convertStringsToChars(targetVar(k)) '.png']))

        % -------------------------------------------------------------------------

        % Set the output GIF file name
        gifVal = fullfile(outputDir,['validation_' convertStringsToChars(targetVar(k)) '.gif']);

        refDates = refValidation.date;
        refData  = refValidation.(targetVar(k));
        synDates = synImages.date;
        dates    = datetime(synDates,'ConvertFrom','yyyyMMdd','format','dd/MM/yyyy');
        synData  = synImages.(targetVar(k));
        bdName   = [convertStringsToChars(targetVar(k)) '_BestDistance'];
        bestDist = synImages.(bdName);
        currentBest = nan(size(synDates));

        % Create an empty figure
        figure('WindowState', 'maximized');
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

                sgtitle(targetVar(k))

                % Create a figure with three subplots
                subplot(3,3,[1,4]);
                img1 = imshow(synthetic);
                colormap(gca, jet(256));
                set(img1, 'AlphaData', ~isnan(synthetic))
                %caxis([0 maxColor])
                caxis([0 6])
                axis equal
                title('Synthetic');
                %colorbar(gca,'southoutside')

                subplot(3,3,[2,5]);
                img2 = imshow(reference);
                colormap(gca, jet(256));
                set(img2, 'AlphaData', ~isnan(synthetic))
                %caxis([0 maxColor])
                caxis([0 6])
                axis equal
                title('Reference');
                h = colorbar(gca,'southoutside');

                % ERROR MAP
                error  = synthetic - reference;
                subplot(3,3,[3,6])
                errMap = imshow(error);
                set(errMap, 'AlphaData', ~isnan(synthetic))
                colormap(gca, coolwarm(256));
                caxis([-2 2])
                title('Error');
                axis equal
                h_err = colorbar(gca,'southoutside');

                % Add a colorbar to the reference image subplot
                %h = colorbar('southoutside');
                set(h, 'Position', [0.13 0.4 0.5 0.03]);
                set(get(h,'label'),'string','Evaporation [mm/day]');
                set(h_err, 'Position', [0.7 0.4 0.205 0.03])
                set(get(h_err,'label'),'string','Evaporation [mm/day]');

                % Best candidate MAE
                subplot(3,3,[7,8,9])
                currentBest(1:i) = bestDist(1:i);
                plot(dates,currentBest);
                hold on
                plot(dates(i),bestDist(i),"Marker","o","Color",'red')
                hold off
                xlim([min(dates) max(dates)])
                ylim([0.05 0.35])
                title(['Best candidate MAE: ' num2str(bestDist(i),'%1.5f')]);
                ylabel('MAE')
                xlabel('Date')
                grid on

                % Set the title of the figure to the name of the images
                if metricV == 1
                    sgtitle({['{\bf\fontsize{14}' char(dates(i)) '}'], ...
                        ['{\fontsize{13}' 'RMSE: ' num2str(validationMetric.(targetVar(k))(i,2),'%1.5f') '}']})
                elseif metricV == 2
                    sgtitle({['{\bf\fontsize{14}' char(dates(i)) '}'], ...
                        ['{\fontsize{13}' 'SPEM: ' num2str(validationMetric.(targetVar(k))(i,2),'%1.5f') '}']})
                elseif metricV == 3
                    sgtitle({['{\bf\fontsize{14}' char(dates(i)) '}'], ...
                        ['{\fontsize{13}' 'SPAEF: ' num2str(validationMetric.(targetVar(k))(i,2),'%1.5f') '}']})
                elseif metricV == 4
                    sgtitle({['{\bf\fontsize{14}' char(dates(i)) '}'], ...
                        ['{\fontsize{13}' 'SPOMF: ' num2str(validationMetric.(targetVar(k))(i,2),'%1.5f') '}']})
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
    end
end

end
