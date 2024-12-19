function visualiseMetrics(nbImages,pixelWise,targetVar,targetDim,refValidation,synImages,validationMetric,sortedDates,metricV,nanValue,varLegend,varRange,errRange,metricKNN,LdateStart,LdateEnd,QdateStart,QdateEnd,daysRange,bootstrap,outDir)

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

targetVarL = lower(targetVar);

for k = 1:numel(targetVar)
    if bootstrap == true
        dates = datetime(cell2mat(validationMetric.(targetVarL(k))(:,1)),'ConvertFrom','yyyyMMdd','Format','dd/MM/yyyy');
        % Validation data
        singleDataVal = cell2mat(validationMetric.(targetVarL(k))(:,3));
        maxValuesVal  = zeros(size(validationMetric.(targetVarL(k)), 1), 1);
        meanValuesVal = zeros(size(validationMetric.(targetVarL(k)), 1), 1);
        minValuesVal  = zeros(size(validationMetric.(targetVarL(k)), 1), 1);
        for i = 1:size(validationMetric.(targetVarL(k)), 1)
            valuesVal = sort(validationMetric.(targetVarL(k)){i, 2});  % Extract values from the second column of the cell
            maxValuesVal(i)  = valuesVal(end);  % Compute the maximum value
            meanValuesVal(i) = mean(valuesVal);
            minValuesVal(i)  = valuesVal(1);  % Compute the minimum value
        end
        % Synthetic data
        varBs = strcat(targetVar(k),'_Bootstrap');
        varBD = strcat(targetVar(k),'_BestDistance');
        refData    = squeeze(mean(mean(refValidation.(targetVarL(k)),1,'omitnan'),2,'omitnan'));  % Extract mean of ref variable
        synData    = squeeze(mean(mean(synImages.(targetVarL(k)),1,'omitnan'),2,'omitnan'));
        maxValues  = zeros(size(validationMetric.(targetVarL(k)), 1), 1);
        meanValues = zeros(size(validationMetric.(targetVarL(k)), 1), 1);
        minValues  = zeros(size(validationMetric.(targetVarL(k)), 1), 1);
        for i = 1:size(validationMetric.(targetVarL(k)), 1)
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
        %nseSynRef = 1-(sum((synData-refData).^2)/sum((synData-mean(synData)).^2));
        alpha = std(synData)/std(refData);
        beta  = mean(synData)/mean(refData);
        kgeSynRef = 1-(sqrt((r-1)^2 + (alpha-1)^2 + (beta-1)^2));
        str = {['KGE: ' num2str(kgeSynRef,'%.5f')] ['r: ' num2str(r,'%.5f') ', \alpha: ' num2str(alpha,'%.5f') ', \beta: ' num2str(beta,'%.5f')]};
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
        saveas(gcf,strcat(outDir,['\bsValidation_AVG_' convertStringsToChars(targetVar(k)) '.png']))

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
            ylabel('1-bSPEM')
        elseif metricKNN == 4
            ylabel('Hellinger distance')
        elseif metricKNN == 5
            ylabel('0.5*(1-bSPEM) + 0.5*Hellinger')
        elseif metricKNN == 6
            ylabel('SPAEF')
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
            title([convertStringsToChars(targetVar(k)) ' - MAE'])
        elseif metricV == 2
            title([convertStringsToChars(targetVar(k)) ' - RMSE'])
        elseif metricV == 3
            title([convertStringsToChars(targetVar(k)) ' - SPEM'])
        elseif metricV == 4
            title([convertStringsToChars(targetVar(k)) ' - SPAEF'])
        elseif metricV == 5
            title([convertStringsToChars(targetVar(k)) ' - KGE'])
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
            ylabel('MAE')
        elseif metricV == 2
            ylabel('RMSE')
        elseif metricV == 3
            ylabel('SPEM')
        elseif metricV == 4
            ylabel('SPAEF')
        elseif metricV == 5
            ylabel('KGE')
        end
        ax = gca;
        ax.YAxis(1).Color = 'k';
        ax.YAxis(2).Color = 'r';
        linkprop(ax.YAxis, 'Limits');
        %legend('Stochastic ensembles','Stochastic mean','Deterministic','Best distance')
        legend('Stochastic ensembles','Deterministic','Best distance')
        set(gcf, 'color', 'white');
        grid on
        saveas(gcf,strcat(outDir,['\bsValidation_RMSE_' convertStringsToChars(targetVar(k)) '.png']))
    else
        if targetDim == 1
            figure('WindowState', 'maximized');
            plot(datetime(validationMetric.(targetVarL(k))(:,1),'ConvertFrom','yyyyMMdd','Format','dd/MM/yyyy'), ...
                refValidation.(targetVarL(k)),'Color','r','DisplayName','Reference');
            hold on
            plot(datetime(validationMetric.(targetVarL(k))(:,1),'ConvertFrom','yyyyMMdd','Format','dd/MM/yyyy'), ...
                synImages.(targetVarL(k)),'Color','b','DisplayName','Synthetic')
            hold off
            if strcmp(startLdate, startQdate)
                subtitle(['Learning periode: ' endQdate '-' endLdate])
            elseif strcmp(endQdate, endLdate)
                subtitle(['Learning periode: ' startLdate '-' startQdate])
            else
                subtitle(['Learning periode: ' startLdate '-' startQdate ' - ' endQdate '-' endLdate])
            end
            xlabel('Date')
            ylabel('Discharge [m^{3}/s]')
            title(targetVar)
            legend()
            set(gcf, 'color', 'white');
            grid on
            saveas(gcf,strcat(outDir,['\' convertStringsToChars(targetVar(k)) '.png']))
        end
        figure('WindowState', 'maximized');
        plot(datetime(validationMetric.(targetVarL(k))(:,1),'ConvertFrom','yyyyMMdd','Format','dd/MM/yyyy'), ...
            validationMetric.(targetVarL(k))(:,2));
        yline(mean(validationMetric.(targetVarL(k))(:,2)),'-',['Mean: ' num2str(mean(validationMetric.(targetVarL(k))(:,2)))],'Color','r')
        %ylim([0 1.4])
        if metricV == 1
            title([convertStringsToChars(targetVar(k)) ' - MAE (mean: ' num2str(mean(validationMetric.(targetVarL(k))(:,2))) ')'])
        elseif metricV == 2
            title([convertStringsToChars(targetVar(k)) ' - RMSE (mean: ' num2str(mean(validationMetric.(targetVarL(k))(:,2))) ')'])
        elseif metricV == 3
            title([convertStringsToChars(targetVar(k)) ' - SPEM (mean: ' num2str(mean(validationMetric.(targetVarL(k))(:,2))) ')'])
        elseif metricV == 4
            title([convertStringsToChars(targetVar(k)) ' - SPAEF (mean: ' num2str(mean(validationMetric.(targetVarL(k))(:,2))) ')'])
        elseif metricV == 5
            title([convertStringsToChars(targetVar(k)) ' - KGE (mean: ' num2str(mean(validationMetric.(targetVarL(k))(:,2))) ')'])
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
            ylabel('MAE')
        elseif metricV == 2
            ylabel('RMSE')
        elseif metricV == 3
            ylabel('SPEM')
        elseif metricV == 4
            ylabel('SPAEF')
        elseif metricV == 5
            ylabel('KGE')
        end
        set(gcf, 'color', 'white');
        grid on
        saveas(gcf,strcat(outDir,['\validation_' convertStringsToChars(targetVar(k)) '.png']))

        % --------------------------------------------------------------------
        if targetDim ~= 1

            refData  = refValidation.(targetVarL(k));
            synData = synImages.(targetVarL(k));
            if ~isnan(nanValue)
                refData(refData == nanValue) = nan;
                synData(isnan(refData)) = nan;
            end

            meanRefData = squeeze(mean(mean(refData,1,'omitnan'),2,'omitnan'));  % Extract mean of ref variable
            meanSynData = squeeze(mean(mean(synData,1,'omitnan'),2,'omitnan'));
            figure('WindowState', 'maximized');
            date = datetime(validationMetric.(targetVarL(k))(:,1),'ConvertFrom','yyyyMMdd','Format','dd/MM/yyyy');
            plot(date, meanRefData, 'r-', date, meanSynData, 'k-');
            legend('Reference','Synthetic','Location','southeast')
            xlabel('Date')
            ylabel(varLegend)
            title(['Mean ' convertStringsToChars(targetVarL(k))])
            r = corr(meanSynData,meanRefData);
            %nseSynRef = 1-(sum((synData-refData).^2)/sum((synData-mean(synData)).^2));
            alpha = std(meanSynData)/std(meanRefData);
            beta  = mean(meanSynData)/mean(meanRefData);
            kgeSynRef = 1-(sqrt((r-1)^2 + (alpha-1)^2 + (beta-1)^2));
            str = {['KGE: ' num2str(kgeSynRef,'%.5f')] ['r: ' num2str(r,'%.5f') ', \alpha: ' num2str(alpha,'%.5f') ', \beta: ' num2str(beta,'%.5f')]};
            subtitle(str)
            grid on
            box off
            %legend boxoff
            set(gcf, 'color', 'white');
            saveas(gcf,strcat(outDir,['\correlation_' convertStringsToChars(targetVar(k)) '.png']))

            % --------------------------------------------------------------------

%             % Set the filter order (adjust as needed)
%             filterOrder = 3;
%             % Set the cutoff frequency for the high-pass filter (adjust as needed)
%             cutoffFrequency = 0.01; % Adjust this value based on your data characteristics
%             % Design a Butterworth high-pass filter
%             [b, a] = butter(filterOrder, cutoffFrequency, 'high');
% 
%             % Apply the filter to both reference and synthetic datasets
%             refDataHighPass = filtfilt(b, a, double(meanRefData));
%             synDataHighPass = filtfilt(b, a, double(meanSynData));
% 
%             % Plot the original and high-pass filtered data
%             figure('WindowState', 'maximized');
%             date = datetime(validationMetric.(targetVarL(k))(:, 1), 'ConvertFrom', 'yyyyMMdd', 'Format', 'dd/MM/yyyy');
%             plot(date, refDataHighPass, 'r-', date, synDataHighPass, 'k-');
%             legend('Reference', 'Synthetic', 'Location', 'southeast');
%             title('High-Pass Filtered Data');
%             xlabel('Date');
%             ylabel(strcat('Detrended ', varLegend));
%             grid on;
%             r = corr(synDataHighPass,refDataHighPass);
%             %nseSynRef = 1-(sum((synData-refData).^2)/sum((synData-mean(synData)).^2));
%             alpha = std(synDataHighPass)/std(refDataHighPass);
%             beta  = mean(synDataHighPass)/mean(refDataHighPass);
%             kgeSynRef = 1-(sqrt((r-1)^2 + (alpha-1)^2 + (beta-1)^2));
%             str = {['KGE: ' num2str(kgeSynRef,'%.5f')] ['r: ' num2str(r,'%.5f') ', \alpha: ' num2str(alpha,'%.5f') ', \beta: ' num2str(beta,'%.5f')]};
%             subtitle(str)
%             grid on
%             box off
%             %legend boxoff
%             set(gcf, 'color', 'white');
%             saveas(gcf,strcat(outDir,['\correlation_' convertStringsToChars(targetVar(k)) '_highpass.png']))

            % --------------------------------------------------------------------
            
            varRefData = squeeze(var(refData,0,[1 2],'omitnan'));
            varSynData = squeeze(var(synData,0,[1 2],'omitnan'));
            figure('WindowState', 'maximized');
            date = datetime(validationMetric.(targetVarL(k))(:,1),'ConvertFrom','yyyyMMdd','Format','dd/MM/yyyy');
            plot(date, varRefData, 'r-', date, varSynData, 'k-');
            legend('Reference','Synthetic','Location','northeast')
            xlabel('Date')
            ylabel(strcat(varLegend,'^2'))
            title([convertStringsToChars(targetVar(k)) ' variance'])
            r = corr(varSynData,varRefData);
            %nseSynRef = 1-(sum((synData-refData).^2)/sum((synData-mean(synData)).^2));
            alpha = std(varSynData)/std(varRefData);
            beta  = mean(varSynData)/mean(varRefData);
            kgeSynRef = 1-(sqrt((r-1)^2 + (alpha-1)^2 + (beta-1)^2));
            str = {['KGE: ' num2str(kgeSynRef,'%.5f')] ['r: ' num2str(r,'%.5f') ', \alpha: ' num2str(alpha,'%.5f') ', \beta: ' num2str(beta,'%.5f')]};
            subtitle(str)
            grid on
            box off
            %legend boxoff
            set(gcf, 'color', 'white');
            saveas(gcf,strcat(outDir,['\variance_' convertStringsToChars(targetVar(k)) '.png']))

            % -------------------------------------------------------------------------
            
            refData   = refValidation.(targetVarL(k));
            synData   = synImages.(targetVarL(k));
            absDayErr = zeros(size(synData));
            dayErr    = absDayErr;
            relErr    = absDayErr;
            for i = 1:size(synData,3)
                absDayErr(:,:,i) = abs(synData(:,:,i) - refData(:,:,i));
                dayErr(:,:,i)    = synData(:,:,i) - refData(:,:,i);
                relErr(:,:,i)    = abs((synData(:,:,i) - refData(:,:,i))./refData(:,:,i));
            end
            meanError  = mean(absDayErr,3);
            meanBias   = mean(dayErr,3);
            meanRelErr = mean(relErr,3);
            
            figure;
            figMean = imagesc(meanError);
            set(figMean, 'AlphaData', ~isnan(meanError))
            colormap(gca, turbo(256));
            caxis([0.2 0.6])
            title('Mean absolute error')
            subtitle(['Mean MAE: ' num2str(mean(mean(meanError,'omitnan'),'omitnan'),'%1.5f')])
            set(gcf, 'color', 'white');
            hcb=colorbar;
            set(get(hcb,'label'),'string','Mean absolute error [mm/day]','Rotation',90);
            axis equal off
            saveas(gcf,strcat(outDir,['\mae_' convertStringsToChars(targetVar(k)) '.png']))
            
            figure;
            figMean = imagesc(meanBias);
            set(figMean, 'AlphaData', ~isnan(meanBias))
            colormap(gca, turbo(256));
            caxis([-0.2 0.2])
            title('Bias')
            subtitle(['Mean bias: ' num2str(mean(mean(meanBias,'omitnan'),'omitnan'),'%1.5f')])
            set(gcf, 'color', 'white');
            hcb=colorbar;
            set(get(hcb,'label'),'string','Bias [mm/day]','Rotation',90);
            axis equal off
            saveas(gcf,strcat(outDir,['\bias_' convertStringsToChars(targetVar(k)) '.png']))
            
            figure;
            figMean = imagesc(meanRelErr);
            set(figMean, 'AlphaData', ~isnan(meanRelErr))
            colormap(gca, turbo(256));
            caxis([0 1])
            title('Mean relative error')
            subtitle(['Mean: ' num2str(mean(mean(meanRelErr,'omitnan'),'omitnan'),'%1.5f')])
            set(gcf, 'color', 'white');
            hcb=colorbar;
            set(get(hcb,'label'),'string','Mean relative error','Rotation',90);
            axis equal off
            saveas(gcf,strcat(outDir,['\mre_' convertStringsToChars(targetVar(k)) '.png']))

            % -------------------------------------------------------------------------
            
            synDates = synImages.date;
            dates = datetime(synDates, 'ConvertFrom', 'yyyyMMdd', 'format', 'dd/MM/yyyy');
            cellData = synImages.(strcat(targetVarL(k), "_Distances")); % Assuming this is a cell array with 25x1 arrays
            minValues = cellfun(@min, cellData);
            maxValues = cellfun(@max, cellData);
            medianValues = cellfun(@median, cellData);

            figure;
            hold on;
            fill([dates; flipud(dates)], [minValues; flipud(maxValues)], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
            for i = 1:numel(cellData)
                scatter(repmat(dates(i), 25, 1), cellData{i}, 20, 'k', 'filled', 'MarkerFaceAlpha', 0.5);
            end
            plot(dates, medianValues, 'r--', 'LineWidth', 1.5);
            xlabel('Date');
            switch metricKNN
                case 1
                    ylabel('RMSE');
                case 2
                    ylabel('MAE');
                case 3
                    ylabel('1-bSPEM');
                case 4
                    ylabel('Hellinger Distance');
                case 5
                    ylabel('0.5*(1-bSPEM) + 0.5*Hellinger');
                case 6
                    ylabel('SPAEF');
            end
            title('Daily distance of the k candidates');
            grid on;
            hold off;

            set(gcf, 'color', 'white');
            saveas(gcf,strcat(outDir,['\distance_' convertStringsToChars(targetVar(k)) '.png']))

            % -------------------------------------------------------------------------

            % Set the output GIF file name
            if pixelWise == false
                gifVal = fullfile(outDir,['\validation_' convertStringsToChars(targetVar(k)) '.gif']);

                refDates = refValidation.date;
                synDates = synImages.date;
                dates    = datetime(synDates,'ConvertFrom','yyyyMMdd','format','dd/MM/yyyy');
                bdName   = [convertStringsToChars(targetVar(k)) '_BestDistance'];
                analogs  = sortedDates(:,2);
                currentDOY  = nan(numel(analogs{1}),numel(dates));
                meanDOY     = nan(numel(dates),1);
                diffDOY     = nan(numel(analogs{1}),1);
                diffBest    = nan(numel(dates),1);
                bestAnalog  = nan(numel(dates),1);
                bestDist    = synImages.(bdName);
                currentBest = nan(size(synDates));
                q1DOY     = [];
                q3DOY     = [];
                minOut      = [];
                maxOut      = [];
                minCurDOY   = nan(numel(dates),1);
                maxCurDOY   = nan(numel(dates),1);

                % Create an empty figure
                figure('WindowState', 'maximized');
                % Loop over each file in the synthetic directory and find the corresponding
                % file in the reference directory
                for i = 1:size(synData,3)
                    % Find the corresponding file in the reference directory with the same name
                    referenceIndex = find(refDates == synDates(i));

                    analogDOY  = day(datetime(sort(analogs{i},'descend'),'ConvertFrom','yyyyMMdd'),'dayofyear');
                    bestAnalog = day(datetime(analogs{i}(1),'ConvertFrom','yyyyMMdd'),'dayofyear');
                    refDOY     = day(datetime(refDates(i),'ConvertFrom','yyyyMMdd'),'dayofyear');
                    % Check for circular transition
                    % Compute the minimum and maximum range boundaries
                    minRange = refDOY - daysRange;
                    maxRange = refDOY + daysRange;
                    % Handle circular transition for minRangeQ
                    if minRange <= 0
                        minRange = 365 + minRange;
                    end
                    % Handle circular transition for maxRangeQ
                    if maxRange > 365
                        maxRange = maxRange - 365;
                    end
                    % Construct the rangeTot array
                    if minRange < maxRange
                        rangeTot = minRange:maxRange;
                    else
                        rangeQmin = minRange:365;
                        rangeQmax = 1:maxRange;
                        rangeTot    = [rangeQmin rangeQmax];
                    end
                    if refDOY == 366
                        refDOY = 1;
                    end
                    rangeTot = rangeTot';
                    % Find the index of refDOY and analogDOY(i) in rangeTot
                    indexRef = find(rangeTot == refDOY);
                    for di = 1:numel(analogDOY)
                        if analogDOY(di) == 366
                            analogDOY(di) = 1;
                        end
                        indexAnalog   = find(rangeTot == analogDOY(di));
                        diffDOY(di)   = indexRef - indexAnalog;
                    end
                    if bestAnalog == 366
                        bestAnalog = 1;
                    end
                    idxBestAnalog = find(rangeTot == bestAnalog);
                    diffBest(i)   = indexRef - idxBestAnalog;

                    % If a matching file is found, display the two images side by side
                    if ~isempty(referenceIndex)
                        % Load the two images
                        synthetic = synData(:,:,i);
                        %synthetic(synthetic==min(min(synthetic))) = NaN;
                        %synthetic(synthetic==-999) = NaN;
                        reference = refData(:,:,referenceIndex);
                        %reference(reference==min(min(reference))) = NaN;
                        %reference(reference==-999) = NaN;

                        sgtitle(targetVar(k))

                        % Create a figure with three subplots
                        subplot(3,3,[1,4]);
                        img1 = imshow(synthetic);
                        colormap(gca, turbo(256));
                        set(img1, 'AlphaData', ~isnan(synthetic))
                        %caxis([0 maxColor])
                        caxis(varRange)
                        axis equal
                        title('Synthetic');
                        %colorbar(gca,'southoutside')

                        subplot(3,3,[2,5]);
                        img2 = imshow(reference);
                        colormap(gca, turbo(256));
                        set(img2, 'AlphaData', ~isnan(synthetic))
                        %caxis([0 maxColor])
                        caxis(varRange)
                        axis equal
                        title('Reference');
                        h = colorbar(gca,'southoutside');

                        % ERROR MAP
                        error  = synthetic - reference;
                        subplot(3,3,[3,6])
                        errMap = imshow(error);
                        set(errMap, 'AlphaData', ~isnan(synthetic))
                        colormap(gca, coolwarm(256));
                        caxis(errRange)
                        title('Error');
                        axis equal
                        h_err = colorbar(gca,'southoutside');

                        % Add a colorbar to the reference image subplot
                        %h = colorbar('southoutside');
                        set(h, 'Position', [0.13 0.4 0.5 0.03]);
                        set(get(h,'label'),'string',varLegend);
                        set(h_err, 'Position', [0.7 0.4 0.205 0.03])
                        set(get(h_err,'label'),'string',varLegend);

                        % DOY difference
                        subplot(3,3,[7,8,9])
                        currentDOY(:,i) = diffDOY;
                        meanDOY(i)      = mean(currentDOY(:,i));
                        %if numel(dates)<60
                        %    boxchart(currentDOY);%,synDates(i));
                        %else
                        if nbImages > 1
                            inBetweenRegionX = [(1:i), fliplr(1:i)];
                            %inBetweenRegionX = [1:numel(dates), fliplr(1:numel(dates))];
                            quartileDOY = quantile(currentDOY(:,i),4);
                            q1DOY = [q1DOY quartileDOY(1)];
                            q3DOY = [q3DOY quartileDOY(3)];
                            inBetweenRegionY = [q3DOY, fliplr(q1DOY)];
                            itqDOY(i) = iqr(currentDOY(:,i));
                            minOut = [minOut, min(currentDOY(currentDOY(:,i) > (q1DOY(i) - itqDOY(i)),i))];
                            maxOut = [maxOut, max(currentDOY(currentDOY(:,i) < (q3DOY(i) + itqDOY(i)),i))];
                            inBetweenRegionY2 = [minOut, fliplr(maxOut)];
                            minCurDOY(i) = min(currentDOY(:,i));
                            maxCurDOY(i) = max(currentDOY(:,i));
                            patch(inBetweenRegionX, inBetweenRegionY, 'k', 'LineStyle', 'none', 'FaceAlpha', 0.15)
                            hold on
                            patch(inBetweenRegionX, inBetweenRegionY2, 'k', 'LineStyle', 'none', 'FaceAlpha', 0.1)
                            hold on
                            plot(1:i,minCurDOY(1:i),'LineStyle','--','Color','k')
                            hold on
                            plot(1:i,maxCurDOY(1:i),'LineStyle','--','Color','k')
                        end
                        %end
                        hold on
                        if nbImages > 1
                            plot(1:i,meanDOY(1:i),"Color",'red')
                        end
                        scatter(1:i,diffBest(1:i),"red")
                        %plot(dates(1:i),meanDOY(1:i),"Color",'red')
                        hold off
                        xlim([1 numel(dates)])
                        %xlim([min(dates) max(dates)])
                        ylim([-daysRange-10 daysRange+10])
                        title(['Mean DOY difference: ' num2str(mean(diffDOY),'%2.0f')]);
                        ylabel('DOY difference')
                        xlabel('Date')
                        grid on
                        ax = gca();
                        %ax.XTick = categorical(1:0.5:numel(dates));
                        %ax.XTickLabels = char(dates);

                        %                 % Best candidate MAE
                        %                 subplot(3,3,[7,8,9])
                        %                 currentBest(1:i) = bestDist(1:i);
                        %                 plot(dates,currentBest);
                        %                 hold on
                        %                 plot(dates(i),bestDist(i),"Marker","o","Color",'red')
                        %                 hold off
                        %                 xlim([min(dates) max(dates)])
                        %                 ylim([0.5 2])
                        %                 title(['Best candidate MAE: ' num2str(bestDist(i),'%1.5f')]);
                        %                 ylabel('MAE')
                        %                 xlabel('Date')
                        %                 grid on

                        % Set the title of the figure to the name of the images
                        if metricV == 1
                            sgtitle({['{\bf\fontsize{14}' char(dates(i)) '}'], ...
                                ['{\fontsize{13}' 'MAE: ' num2str(validationMetric.(targetVarL(k))(i,2),'%1.5f') '}']})
                        elseif metricV == 2
                            sgtitle({['{\bf\fontsize{14}' char(dates(i)) '}'], ...
                                ['{\fontsize{13}' 'RMSE: ' num2str(validationMetric.(targetVarL(k))(i,2),'%1.5f') '}']})
                        elseif metricV == 3
                            sgtitle({['{\bf\fontsize{14}' char(dates(i)) '}'], ...
                                ['{\fontsize{13}' 'SPEM: ' num2str(validationMetric.(targetVarL(k))(i,2),'%1.5f') '}']})
                        elseif metricV == 4
                            sgtitle({['{\bf\fontsize{14}' char(dates(i)) '}'], ...
                                ['{\fontsize{13}' 'SPAEF: ' num2str(validationMetric.(targetVarL(k))(i,2),'%1.5f') '}']})
                        elseif metricV == 5
                            sgtitle({['{\bf\fontsize{14}' char(dates(i)) '}'], ...
                                ['{\fontsize{13}' 'KGE: ' num2str(validationMetric.(targetVarL(k))(i,2),'%1.5f') '}']})
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
            else

                % -------------------------------------------------------------------------

                % Set the output GIF file name
                gifVal = fullfile(outDir,['\validation_' convertStringsToChars(targetVar(k)) '.gif']);

                refDates = refValidation.date;
                synDates = synImages.date;
                dates    = datetime(synDates,'ConvertFrom','yyyyMMdd','format','dd/MM/yyyy');

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
                        %synthetic(synthetic==min(min(synthetic))) = NaN;
                        %synthetic(synthetic==-999) = NaN;
                        reference = refData(:,:,referenceIndex);
                        %reference(reference==min(min(reference))) = NaN;
                        %reference(reference==-999) = NaN;

                        sgtitle(targetVar(k))

                        % Create a figure with three subplots
                        subplot(1,3,1);
                        img1 = imshow(synthetic);
                        colormap(gca, turbo(256));
                        set(img1, 'AlphaData', ~isnan(synthetic))
                        %caxis([0 maxColor])
                        caxis(varRange)
                        axis equal
                        title('Synthetic');
                        %colorbar(gca,'southoutside')

                        subplot(1,3,2);
                        img2 = imshow(reference);
                        colormap(gca, turbo(256));
                        set(img2, 'AlphaData', ~isnan(synthetic))
                        %caxis([0 maxColor])
                        caxis(varRange)
                        axis equal
                        title('Reference');
                        h = colorbar(gca,'southoutside');

                        % ERROR MAP
                        error  = synthetic - reference;
                        subplot(1,3,3)
                        errMap = imshow(error);
                        set(errMap, 'AlphaData', ~isnan(synthetic))
                        colormap(gca, coolwarm(256));
                        caxis(errRange)
                        title('Error');
                        axis equal
                        h_err = colorbar(gca,'southoutside');

                        % Add a colorbar to the reference image subplot
                        %h = colorbar('southoutside');
                        set(h, 'Position', [0.13 0.2 0.5 0.03]);
                        set(get(h,'label'),'string',varLegend);
                        set(h_err, 'Position', [0.7 0.2 0.205 0.03])
                        set(get(h_err,'label'),'string',varLegend);

                        % Set the title of the figure to the name of the images
                        if metricV == 1
                            sgtitle({['{\bf\fontsize{14}' char(dates(i)) '}'], ...
                                ['{\fontsize{13}' 'MAE: ' num2str(validationMetric.(targetVarL(k))(i,2),'%1.5f') '}']})
                        elseif metricV == 2
                            sgtitle({['{\bf\fontsize{14}' char(dates(i)) '}'], ...
                                ['{\fontsize{13}' 'RMSE: ' num2str(validationMetric.(targetVarL(k))(i,2),'%1.5f') '}']})
                        elseif metricV == 3
                            sgtitle({['{\bf\fontsize{14}' char(dates(i)) '}'], ...
                                ['{\fontsize{13}' 'SPEM: ' num2str(validationMetric.(targetVarL(k))(i,2),'%1.5f') '}']})
                        elseif metricV == 4
                            sgtitle({['{\bf\fontsize{14}' char(dates(i)) '}'], ...
                                ['{\fontsize{13}' 'SPAEF: ' num2str(validationMetric.(targetVarL(k))(i,2),'%1.5f') '}']})
                        elseif metricV == 5
                            sgtitle({['{\bf\fontsize{14}' char(dates(i)) '}'], ...
                                ['{\fontsize{13}' 'KGE: ' num2str(validationMetric.(targetVarL(k))(i,2),'%1.5f') '}']})
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
end

end
